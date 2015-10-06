/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaTelluricWavelengthCorrection
 Version: 1.0
 Description: Apply wavelength correction based on telluric lines
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2015
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralTools.h"
#include "libraries/operaStats.h"							// for MAXORDERS
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFit.h"							// for operaFitSplineDouble
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include "core-espadons/operaTelluricWavelengthCorrection.h"

/*! \file operaTelluricWavelengthCorrection.cpp */

using namespace std;

operaArgumentHandler args;

/*! 
 * operaTelluricWavelengthCorrection
 * \author Eder Martioli
 * \brief Calculate and apply wavelength correction based on telluric lines.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

// A simple data structure to work with (wavelength, intensity) pairs
struct Spectrum {
	vector <double> wavelength;
	vector <double> intensity;
	Spectrum() { }
	Spectrum(unsigned presize) : wavelength(presize), intensity(presize) { }
	unsigned int size() const { return wavelength.size(); }
	bool empty() const { return wavelength.empty(); }
	void add(double wavelength, double intensity) { this->wavelength.push_back(wavelength); this->intensity.push_back(intensity); }
	void reverse() { std::reverse(wavelength.begin(), wavelength.end()); std::reverse(intensity.begin(), intensity.end()); }
	void resize(unsigned newsize) { wavelength.resize(newsize); intensity.resize(newsize); }
	void shrink_to_fit() { vector <double>(wavelength.begin(), wavelength.end()).swap(wavelength); vector <double>(intensity.begin(), intensity.end()).swap(intensity); }
	double* wl_ptr(unsigned offset = 0) { return &(*(wavelength.begin()+offset)); }
	double* intensity_ptr(unsigned offset = 0) { return &(*(intensity.begin()+offset)); }
	const double* wl_ptr(unsigned offset = 0) const { return &(*(wavelength.begin()+offset)); }
	const double* intensity_ptr(unsigned offset = 0) const { return &(*(intensity.begin()+offset)); }
	double firstwl() const { return wavelength.front(); }
	double midwl() const { return wavelength[wavelength.size()/2]; }
	double lastwl() const { return wavelength.back(); }
};

int main(int argc, char *argv[])
{    
	string inputWaveFile;
    string inputObjectSpectrum;
    string inputFlatFluxCalibration;
	string outputWaveFile;
    string telluric_lines; // HITRAN Library
    string telluric_spectrum; // A telluric spectrum
    string inputWavelengthMaskForTelluric;

	int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    bool StarPlusSky = false;
    
    // The parameters below we don't know yet whether they would be useful if used as input
    double spectralResolution = 80000;
    double radialVelocityRange = 10;
    double radialVelocityStep = 0.05;
    double XCorrelationThreshold = 0.05;
    unsigned normalizationBinsize = 110;
    bool useFitToFindMaximum = false;

    unsigned RVCorrectionMethod = 1; // 1. line matching; 2. x-correlation;
    
    string rvcorrsplotfilename;
    string specplotfilename;
    string rvcorrscriptfilename;
    string specscriptfilename;
    string rvcorrdatafilename;
    string specdatafilename;
    string rvcorrfitdatafilename;
        
    args.AddRequiredArgument("inputWaveFile", inputWaveFile, "input wavelength calibration file (.wcal)");
	args.AddRequiredArgument("inputObjectSpectrum", inputObjectSpectrum, "input object spectrum file (.e or .p)");
    args.AddOptionalArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "", "flat field spectrum ff_");
    args.AddRequiredArgument("outputWaveFile", outputWaveFile, "output telluric wavelength calibration file (.tell)");
    args.AddRequiredArgument("telluric_lines", telluric_lines, "atlas of telluric lines (HITRAN)");
    args.AddOptionalArgument("telluric_spectrum", telluric_spectrum, "", "spectrum of telluric lines");
    args.AddRequiredArgument("inputWavelengthMaskForTelluric", inputWavelengthMaskForTelluric, "telluric wavelength mask");
    
    args.AddRequiredArgument("spectralResolution", spectralResolution, "input spectral resolution (wl/dwl) as reference for line detection");
    args.AddRequiredArgument("radialVelocityRange", radialVelocityRange, "radial Velocity Range (in km/s) to scan for first order correction");
    args.AddRequiredArgument("radialVelocityStep", radialVelocityStep, "radial Velocity Step step (in km/s) to scan for first order correction");
    args.AddRequiredArgument("XCorrelationThreshold", XCorrelationThreshold, "X-correlation lower threshold to consider a match between telluric and object spectra");
    args.AddRequiredArgument("normalizationBinsize", normalizationBinsize, "binsize to normalize input object spectrum");
    args.AddSwitch("useFitToFindMaximum", useFitToFindMaximum, "use gaussian fit to find RV correction");
    
    args.AddRequiredArgument("RVCorrectionMethod", RVCorrectionMethod, "Method for measuring RV correction: 1 = Line matching (default), 2 = Spectral cross-correlation");
    
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddSwitch("StarPlusSky", StarPlusSky, "star plus sky mode");
    
    args.AddOptionalArgument("rvcorrsplotfilename", rvcorrsplotfilename, "", "Output RV correction plot eps file name");
    args.AddOptionalArgument("specplotfilename", specplotfilename, "", "Output spectrum plot eps file name");
    args.AddOptionalArgument("rvcorrscriptfilename", rvcorrscriptfilename, "", "Output cross-correlation gnuplot script file name");
    args.AddOptionalArgument("specscriptfilename", specscriptfilename, "", "Output spectrum gnuplot script file name");
    args.AddOptionalArgument("rvcorrdatafilename", rvcorrdatafilename, "", "Output RV correction data file name");
    args.AddOptionalArgument("specdatafilename", specdatafilename, "", "Output spectrum data file name");
	args.AddOptionalArgument("rvcorrfitdatafilename", rvcorrfitdatafilename, "", "Output RV correction fit data file name");
	
	try {
		args.Parse(argc, argv);
		
		// we need an input wavelength calibration file ...
		if (inputWaveFile.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input object spectrum file ...        
		if (inputObjectSpectrum.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output wavelength calibration file ...
		if (outputWaveFile.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input atlas of telluric lines ...
		if (telluric_lines.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        if (inputWavelengthMaskForTelluric.empty()) {
            throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
		if (args.verbose) {
			cout << "operaTelluricWavelengthCorrection: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaTelluricWavelengthCorrection: inputObjectSpectrum = " << inputObjectSpectrum << endl;
            cout << "operaTelluricWavelengthCorrection: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaTelluricWavelengthCorrection: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaTelluricWavelengthCorrection: telluric_lines =" << telluric_lines << endl;
			cout << "operaTelluricWavelengthCorrection: telluric_spectrum =" << telluric_spectrum << endl;
			cout << "operaTelluricWavelengthCorrection: spectralResolution =" << spectralResolution << endl;
			cout << "operaTelluricWavelengthCorrection: radialVelocityRange =" << radialVelocityRange << endl;
			cout << "operaTelluricWavelengthCorrection: radialVelocityStep =" << radialVelocityStep << endl;
			cout << "operaTelluricWavelengthCorrection: XCorrelationThreshold =" << XCorrelationThreshold << endl;
            cout << "operaTelluricWavelengthCorrection: normalizationBinsize =" << normalizationBinsize << endl;
            cout << "operaTelluricWavelengthCorrection: StarPlusSky = " << StarPlusSky << endl;
            cout << "operaTelluricWavelengthCorrection: useFitToFindMaximum = " << useFitToFindMaximum << endl;
            cout << "operaTelluricWavelengthCorrection: inputWavelengthMaskForTelluric = " << inputWavelengthMaskForTelluric << endl;
            cout << "operaTelluricWavelengthCorrection: RVCorrectionMethod = " << RVCorrectionMethod << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaTelluricWavelengthCorrection: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaTelluricWavelengthCorrection: rvcorrsplotfilename = " << rvcorrsplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specplotfilename = " << specplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: rvcorrscriptfilename = " << rvcorrscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specscriptfilename = " << specscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: rvcorrdatafilename = " << rvcorrdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: rvcorrfitdatafilename = " << rvcorrfitdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: specdatafilename = " << specdatafilename << endl;
            }
		}
		
		ofstream frvcorrdata;
        ofstream frvcorrfitdata;
        ofstream fspecdata;
        if (!rvcorrdatafilename.empty()) frvcorrdata.open(rvcorrdatafilename.c_str());
        if (!rvcorrfitdatafilename.empty()) frvcorrfitdata.open(rvcorrfitdatafilename.c_str());
        if (!specdatafilename.empty()) fspecdata.open(specdatafilename.c_str());
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputObjectSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile); // This merges in the wavelength calibration information

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaTelluricWavelengthCorrection: minorder ="<< minorder << " maxorder=" << maxorder << endl;
    
        // Correct for flat-field
        if (!inputFlatFluxCalibration.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputFlatFluxCalibration);
            bool starplusskyInvertSkyFiber = false;
            spectralOrders.correctFlatField(minorder, maxorder, StarPlusSky, starplusskyInvertSkyFiber);
        }
        
        bool emissionSpectrum = false;
        
        double LocalMaxFilterWidth=4.0;
        double MinPeakDepth=3;
        double DetectionThreshold=0.05;
        double nsigclip=3.0;
       
        // Below it detects absorption lines in the observed spectrum within telluric regions given by the inputWavelengthMaskForTelluric
        Spectrum telluricLinesFromObject(MAXORDERS*MAXREFWAVELENGTHSPERORDER);
        vector <double> telluricLineWidthsFromObject(MAXORDERS*MAXREFWAVELENGTHSPERORDER);
        unsigned nlines = spectralOrders.detectSpectralLinesWithinWavelengthMask(inputWavelengthMaskForTelluric, minorder, maxorder,true, normalizationBinsize, telluricLinesFromObject.wl_ptr(), telluricLinesFromObject.intensity_ptr(), &(*telluricLineWidthsFromObject.begin()),spectralResolution,emissionSpectrum,LocalMaxFilterWidth,MinPeakDepth,DetectionThreshold,nsigclip);
        telluricLinesFromObject.resize(nlines);
        telluricLinesFromObject.shrink_to_fit();

        if(args.debug){
            for (unsigned l=0; l<nlines; l++) {
                cout << telluricLinesFromObject.wavelength[l] << " " << telluricLinesFromObject.intensity[l] << " " << telluricLineWidthsFromObject[l] << endl;
            }
        }
        // ---
        
        // Get object spectrum within wavelength mask input ranges defined in "inputWavelengthMaskForTelluric"
        Spectrum objectSpectrum(MAXORDERS*MAXSPECTRALELEMENTSPERORDER);
        vector <double> objectSpectrumVariance(MAXORDERS*MAXSPECTRALELEMENTSPERORDER); //not really used for anything...
        unsigned nelem = spectralOrders.getSpectrumWithinTelluricMask(inputWavelengthMaskForTelluric, minorder, maxorder, true, normalizationBinsize, objectSpectrum.wl_ptr(), objectSpectrum.intensity_ptr(), &(*objectSpectrumVariance.begin()));
        objectSpectrum.resize(nelem);
        objectSpectrum.shrink_to_fit();

		// Read telluric lines database lambda vs. intensity
		Spectrum telluricLines;
        if (args.debug) cout << "operaTelluricWavelengthCorrection: reading telluric lines database " << telluric_lines << endl;
		readTelluricLines(telluric_lines, telluricLines);
    
        
        // Initialize rvshift to zero, so if telluric SNR is low the RV shift will not do anything.
        float rvshift = 0;
        float rvshifterror = 0;
        
        /* 
         *   Below one can use one of the following methods to measure the RV shift :
         * (RVCorrectionMethod=1): detect lines in observed spectrum and match them with telluric HITRAN lines
         * (RVCorrectionMethod=2): cross-correlation between a synthetic telluric and the observed spectra
         */
        
        // Spectrum plot: plot observed and reference telluric spectra.
        if (fspecdata.is_open()) {
            telluricSpectrumComparison(telluricLines, objectSpectrum, telluric_spectrum, radialVelocityRange, spectralResolution, fspecdata);
            fspecdata.close();
        }
        
        switch (RVCorrectionMethod) {
                
            case 1: {
                if (args.verbose) cout << "operaTelluricWavelengthCorrection: calculating RV shift by matching lines for radialVelocityRange=" << radialVelocityRange << " km/s and radialVelocityStep=" << radialVelocityStep << " km/s" << endl;
                
                // Note: spectralResolution can be updated with measurements from Object spectrum! -- could use median of all orders
                double *telluricMatchedLines = new double[MAXORDERS*MAXREFWAVELENGTHSPERORDER];
                double *objectMatchedLines = new double[MAXORDERS*MAXREFWAVELENGTHSPERORDER];
                double *objectMatchedLineIntensities = new double[MAXORDERS*MAXREFWAVELENGTHSPERORDER];
                float *radialVelocities = new float[MAXORDERS*MAXREFWAVELENGTHSPERORDER];
                
                unsigned nmatchedLines = matchTelluricLines(telluricLines, telluricLinesFromObject, telluricMatchedLines, objectMatchedLines, objectMatchedLineIntensities, radialVelocities, spectralResolution);
                
                if (frvcorrdata.is_open()) {
                    for (unsigned l=0; l<nmatchedLines; l++) {
                        double linewidth = telluricMatchedLines[l]/spectralResolution;
                        frvcorrdata << l << " " << telluricMatchedLines[l] << " " <<  objectMatchedLines[l] << " " << objectMatchedLineIntensities[l] << " " << linewidth << " " << radialVelocities[l] << endl;
                    }
                    frvcorrdata.close();
                }
                
                unsigned histogramVectorSize = (unsigned)ceil(radialVelocityRange/radialVelocityStep);
                
                float *rvVector = new float[histogramVectorSize];
                float *probDensity = new float[histogramVectorSize];
                float *histWaveVector = new float[histogramVectorSize];
                
                unsigned nPointsInHistogram = generateHistogramData(nmatchedLines, telluricMatchedLines, radialVelocities, radialVelocityRange, radialVelocityStep, rvVector, probDensity, histWaveVector);
                
                if (frvcorrfitdata.is_open()) {
                    for (unsigned i=0; i<nPointsInHistogram; i++) {
                        frvcorrfitdata << i << " " << histWaveVector[i] << " " <<  rvVector[i] << " " << probDensity[i] << endl;
                    }
                    frvcorrfitdata.close();
                }
                
                unsigned minNumberOfMatchedLines = 10; // arbitrary threshold to avoid small number statistics.
                
                if(nmatchedLines > minNumberOfMatchedLines) {
                    rvshift = operaArrayMedian(nmatchedLines,radialVelocities);
                    rvshifterror = operaArrayMedianSigma(nmatchedLines,radialVelocities,rvshift);
                }
                
                if (args.verbose) cout << "operaTelluricWavelengthCorrection: (Line Match Method) Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s" << endl;
                
                // Telluric spectrum plot:
                if (!specdatafilename.empty() && !rvcorrdatafilename.empty() && !specscriptfilename.empty()) {
                    GenerateTelluricSpecAndLinesPlot(specscriptfilename, specplotfilename, specdatafilename, rvcorrdatafilename);
                }
                
                if (!rvcorrdatafilename.empty() && !rvcorrfitdatafilename.empty() && !rvcorrscriptfilename.empty()) {
                    GenerateTelluricRVCorrPlot(rvcorrscriptfilename, rvcorrsplotfilename, rvcorrdatafilename, rvcorrfitdatafilename, rvshift, rvshifterror);
                }
                
                break;
            }
            case 2: {
                
                /* -- E. Martioli Aug 17 2015 -- below is the old way of calculating RVshift for telluric wavelength
                 correction, where it uses a cross-correlation between observed and telluric spectra rather than
                 individual line positions. This method seems to be biased by~300m/s because the line profiles are
                 not symmetric.
                 */
                
                if (args.verbose) cout << "operaTelluricWavelengthCorrection: calculating cross-correlation for radialVelocityRange=" << radialVelocityRange << " km/s and radialVelocityStep=" << radialVelocityStep << " km/s" << endl;
                double maxcorr=-BIG, chisqr=0, rvshift_tmp=0, rvshifterror_tmp=0;
                
                bool validXCorrelation = calculateRVShiftByXCorr(telluricLines, objectSpectrum, radialVelocityRange, radialVelocityStep, XCorrelationThreshold, rvshift_tmp, rvshifterror_tmp, maxcorr, frvcorrdata, frvcorrfitdata, spectralResolution, useFitToFindMaximum, chisqr);
                
                if(validXCorrelation) {
                    rvshift = (float)rvshift_tmp;
                    rvshifterror = (float)rvshifterror_tmp;
                }
                
                if(args.verbose) cout << "operaTelluricWavelengthCorrection: (X-Corr Method) Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s, maxXCorr=" << maxcorr << ", chisqr=" << chisqr << "\n" << endl;
                
                // Telluric spectrum plot:
                if (!specdatafilename.empty() && !specscriptfilename.empty()) {
                    GenerateTelluricSpecPlot(specscriptfilename, specplotfilename, specdatafilename);
                }
                
                // Telluric wavelength correction plot:
                if (frvcorrdata.is_open() && frvcorrfitdata.is_open()) {
                    frvcorrdata.close();
                    frvcorrfitdata.close();
                    if (!rvcorrscriptfilename.empty()) {
                        GenerateTelluricXCorrelationPlot(rvcorrscriptfilename, rvcorrsplotfilename, rvcorrdatafilename, rvcorrfitdatafilename);
                    }
                }
                
                break;
            }
                
            default:
                break;
        }
        
        //Apr 15, 2015 CU - Output rvel shift instead of wcal file
        spectralOrders.setTelluricRadialVelocityCorrection((double)rvshift);
        operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputWaveFile, Tell);

    }
	catch (operaException e) {
		cerr << "operaTelluricWavelengthCorrection: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaTelluricWavelengthCorrection: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

// Generate multiple plots containing statistical info about telluric wavelength correction
void GenerateTelluricXCorrelationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string cleanDataFileName)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
	
    fgnu << "reset" << endl;

    fgnu << "\nset xlabel \"Radial Velocity (km/s)\"" << endl;
    fgnu << "set ylabel \"cross-correlation\"" << endl;
    fgnu << "set pointsize 1.5" << endl;

    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "plot \"" << dataFileName << "\" u 1:2 t \"gaussian fit\" w l lt 3 lw 2, ";
        fgnu << "\"" << dataFileName << "\" u 1:3:4 t \"XCorr data\" w yerr pt 6, ";
        fgnu << "\"" << cleanDataFileName << "\" u 1:4:2:5 t \"fit data\" w xyerr pt 7" << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    } else {
        fgnu << endl;
        
        fgnu << "plot \"" << dataFileName << "\" u 1:2 t \"gaussian fit\" w l lt 3 lw 2, ";
        fgnu << "\"" << dataFileName << "\" u 1:3:4 t \"XCorr data\" w yerr pt 6, ";
        fgnu << "\"" << cleanDataFileName << "\" u 1:4:2:5 t \"fit data\" w xyerr pt 7" << endl;
        
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();

    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate multiple plots containing statistical info about telluric wavelength correction
void GenerateTelluricRVCorrPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string histDataFileName, float rvshift, float rvshifterror)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "\nreset" << endl;

    if(!outputPlotEPSFileName.empty()) {
        
        fgnu << "\nNX=1; NY=2" << endl;
        fgnu << "\nDX=0.1; DY=0.1; SX=0.98; SY=0.42" << endl;
        fgnu << "\nset bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
        
        fgnu << "\nset size SX*NX+DX*2,SY*NY+DY*2" << endl;
        
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nmedianRV = " << rvshift << endl;
        fgnu << "\nRVerr = " << rvshifterror << endl;
        
        fgnu << "\nset multiplot" << endl;
        
        fgnu << "\nset size 0.87*SX,1.35*SY" << endl;
        fgnu << "set origin DX,DY+2*SY/3" << endl;
        fgnu << "set xlabel \"Radial velocity (km/s)\"" << endl;
        fgnu << "set ylabel \"Probability Density\"" << endl;
        
        fgnu << "\nset label \"{/Symbol D}RV = " << fixed << setprecision(3) <<  rvshift << "+/-" << rvshifterror << "\" at " << rvshift + 2*rvshifterror << ",0.125" << endl;
        
        fgnu << "\nset arrow from medianRV,0.1 to medianRV,0 lt 3 lw 2" << endl;
        fgnu << "set arrow from medianRV-RVerr,0 to medianRV-RVerr,0.1 nohead lt 2 lw 1" << endl;
        fgnu << "set arrow from medianRV+RVerr,0 to medianRV+RVerr,0.1 nohead lt 2 lw 1" << endl;
        
        fgnu << "\nplot \"" << histDataFileName << "\" u 3:4 t \"RV shift PDF\" w boxes" << endl;
        
        fgnu << "\nset size 0.87*SX,0.45*SY" << endl;
        fgnu << "set origin DX,DY" << endl;
        fgnu << "set xlabel \"{/Symbol l} (nm)\"" << endl;
        fgnu << "set ylabel \"Radial velocity (km/s)\"" << endl;
        
        fgnu << "\nf(x) = medianRV" << endl;
        fgnu << "fl(x) = medianRV - RVerr" << endl;
        fgnu << "fh(x) = medianRV + RVerr" << endl;
        
        fgnu << "\nset yrange[medianRV - 8*RVerr:medianRV + 8*RVerr]" << endl;
        
        fgnu << "\nplot \"" << dataFileName << "\" u 2:6 notitle w p pt 7, f(x) notitle w l lt 3 lw 2, fl(x) notitle w l lt 2 lw 1, fh(x) notitle w l lt 2 lw 1" << endl;
        
        fgnu << "\nunset multiplot" << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    
    } else {
        
        fgnu << "\nNX=1; NY=2" << endl;
        fgnu << "\nDX=0.1; DY=0.1; SX=0.98; SY=0.42" << endl;
        fgnu << "\nset bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
        
        fgnu << "\nset size SX*NX+DX*2,SY*NY+DY*2" << endl;
        
        fgnu << "\nmedianRV = " << rvshift << endl;
        fgnu << "\nRVerr = " << rvshifterror << endl;
        
        fgnu << "\nset multiplot" << endl;
        
        fgnu << "\nset size 0.87*SX,1.35*SY" << endl;
        fgnu << "set origin DX,DY+2*SY/3" << endl;
        fgnu << "set xlabel \"Radial velocity (km/s)\"" << endl;
        fgnu << "set ylabel \"Probability Density\"" << endl;
        
        fgnu << "\nset label \"{/Symbol D}RV = " << fixed << setprecision(3) << rvshift << "+/-" << rvshifterror << "\" at " << rvshift + 2*rvshifterror << ",0.125" << endl;
        
        fgnu << "\nset arrow from medianRV,0.1 to medianRV,0 lt 3 lw 2" << endl;
        fgnu << "set arrow from medianRV-RVerr,0 to medianRV-RVerr,0.1 nohead lt 2 lw 1" << endl;
        fgnu << "set arrow from medianRV+RVerr,0 to medianRV+RVerr,0.1 nohead lt 2 lw 1" << endl;
        
        fgnu << "\nplot \"" << histDataFileName << "\" u 3:4 t \"RV shift PDF\" w boxes" << endl;
        
        fgnu << "\nset size 0.87*SX,0.45*SY" << endl;
        fgnu << "set origin DX,DY" << endl;
        fgnu << "set xlabel \"{/Symbol l} (nm)\"" << endl;
        fgnu << "set ylabel \"Radial velocity (km/s)\"" << endl;
        
        fgnu << "\nf(x) = medianRV" << endl;
        fgnu << "fl(x) = medianRV - RVerr" << endl;
        fgnu << "fh(x) = medianRV + RVerr" << endl;
        
        fgnu << "\nset yrange[medianRV - 8*RVerr:medianRV + 8*RVerr]" << endl;
        
        fgnu << "\nplot \"" << dataFileName << "\" u 2:6 notitle w p pt 7, f(x) notitle w l lt 3 lw 2, fl(x) notitle w l lt 2 lw 1, fh(x) notitle w l lt 2 lw 1" << endl;
        
        fgnu << "\nunset multiplot" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    }
    
    fgnu.close();
    
    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate 2D plot for spectra of atlas + comparison + identified lines
void GenerateTelluricSpecAndLinesPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename, string linesdatafilename)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "#unset key" << endl;
    
    fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
    
    fgnu << "set ylabel \"norm flux\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 1, ";
        fgnu << "\"" << specdatafilename << "\" u 1:3 t \"Telluric Reference (HITRAN)\" w l lt 3, ";
        fgnu << "\"" << linesdatafilename << "\" u 3:(1 - (1-$4)/2):5 t \"Matched lines\" w xerr pt 7 lt 2" << endl;

        fgnu << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    } else {
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 1, ";
        fgnu << "\"" << specdatafilename << "\" u 1:3 t \"Telluric Reference (HITRAN)\" w l lt 3, ";
        fgnu << "\"" << linesdatafilename << "\" u 3:(1 - (1-$4)/2):5 t \"Matched lines\" w xerr pt 7 lt 2" << endl;
        
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate 2D plot for spectra of atlas + comparison
void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename)
{    
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "#unset key" << endl;

    fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;

    fgnu << "set ylabel \"norm flux\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 4, ";
        fgnu << "\"" << specdatafilename << "\" u 1:3 t \"Telluric Reference (HITRAN)\" w l lt 3";
        fgnu << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    } else {
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 4, ";
        fgnu << "\"" << specdatafilename << "\" u 1:3 t \"Telluric Reference (HITRAN)\" w l lt 3";
        fgnu << endl;
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate telluric spectra from HITRAN and KPNO and write both to filestream
void telluricSpectrumComparison(const Spectrum& telluricLines, const Spectrum& objectSpectrum, string telluric_spectrum, double radialVelocityRange, double spectralResolution, ofstream& fspecdata)
{
	// Find line with maximum absorption for normalization
	double maxatlasflux = -BIG;
	if (!telluricLines.empty()) {
		for (unsigned i=0; i<telluricLines.size(); i++) {
			if(1 - telluricLines.intensity[i] > maxatlasflux) maxatlasflux = 1 - telluricLines.intensity[i];
		}
	} else {
		maxatlasflux = 1;
	}
    
	// Define wavelength full range
	double wl0 = objectSpectrum.firstwl() - radialVelocityRange*objectSpectrum.firstwl()/SPEED_OF_LIGHT_KMS;
	double wlf = objectSpectrum.lastwl() + radialVelocityRange*objectSpectrum.lastwl()/SPEED_OF_LIGHT_KMS;
	if(args.verbose) cout << "operaTelluricWavelengthCorrection: wl0=" <<  wl0 << " wlf=" << wlf << endl;
	
	// First read telluric spectrum from input spectrum file.
	Spectrum telluricSpectrum;
	if (!telluric_spectrum.empty()) readTelluricSpectrum(telluric_spectrum, telluricSpectrum);
	
	// Then use hitran lines to generate synthetic spectrum
	vector <double> hitranTelluricSpectrum;
	generateSyntheticTelluricSpectrumUsingGaussianProfile(telluricLines, objectSpectrum.wavelength, hitranTelluricSpectrum, spectralResolution);
	//generateSyntheticTelluricSpectrumUsingVoigtProfile(telluricLines, objectSpectrum.wavelength, hitranTelluricSpectrum, spectralResolution);
    
	// Resample the telluricSpectrum to the wavelengths in objectSpectrum
	vector<double> KPNOTelluricSpectrum;
	if (!telluricSpectrum.empty()) {
		unsigned startindex, endindex;
		getWavelengthSubrange(telluricSpectrum.wavelength, wl0, wlf, startindex, endindex);
		
		// Normalize telluric reference
		vector<double> transmission;
		for (unsigned i=startindex; i<endindex; i++) {
			transmission.push_back(1 - (1-telluricSpectrum.intensity[i])/maxatlasflux);
		}
		
		KPNOTelluricSpectrum.resize(objectSpectrum.size());
		
		// Spline to match same sampling as in observed spectrum
		operaFitSplineDouble(transmission.size(),telluricSpectrum.wl_ptr(startindex),&(*transmission.begin()),objectSpectrum.size(),objectSpectrum.wl_ptr(),&(*KPNOTelluricSpectrum.begin()));
	}
	
	// Print spectral data to file
	for(unsigned i=0; i<objectSpectrum.size(); i++) {
		fspecdata << objectSpectrum.wavelength[i] << " " << objectSpectrum.intensity[i] << " " << hitranTelluricSpectrum[i];
		if (!KPNOTelluricSpectrum.empty()) fspecdata << " " << KPNOTelluricSpectrum[i];
		fspecdata << endl;
	}
}

// Read the the full atmospheric transmission spectrum
void readTelluricSpectrum(string telluric_spectrum, Spectrum& telluricSpectrum)
{
	igzstream astream(telluric_spectrum.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') { // skip blank lines and comments
				double tmpwl;
				double tmpi;
				if(!sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi)) continue; //skip over bad line
				telluricSpectrum.add(tmpwl, tmpi);
            }
		}
		astream.close();
	}
	if (args.verbose) {
		if (telluricSpectrum.empty()) printf("          [Telluric] no points found in telluric spectrum.\n");
		else printf("          [Telluric] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", telluricSpectrum.size(), telluricSpectrum.firstwl(), telluricSpectrum.midwl(), telluricSpectrum.lastwl());
	}
}

// Read the entire set of telluric lines in HITRAN database
void readTelluricLines(string telluric_database_file, Spectrum& telluricLines)
{
	const double N_OVER_V = TYPICAL_PRESSURE_AT_MAUNAKEA/(TYPICAL_TEMPERATURE_AT_MAUNAKEA*k_BOLTZMANN_CONSTANT);
	
    igzstream astream(telluric_database_file.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') { // skip blank lines and comments
				double wave_number;
				float intensity;
				if(!sscanf(dataline.c_str(), "%*d %lf %G %*[^\n]", &wave_number, &intensity)) continue; //skip over bad line
				double wavelength_in_nm = 1e7/wave_number;
				telluricLines.add(convertVacuumToAirWavelength(wavelength_in_nm*10)/10, ((double)intensity/(N_OVER_V*1e-6))/TYPICAL_ATMOSPHERE_PATH_LENGTH);
			}
		}
		astream.close();
	}
	telluricLines.reverse();
    if (args.verbose) {
		if (telluricLines.empty()) printf("          [Telluric] no lines found in telluric database.\n");
		else printf("          [Telluric] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", telluricLines.size(), telluricLines.firstwl(), telluricLines.midwl(), telluricLines.lastwl());
	}
}

// Finds the first and last index between wl0 and wlf. Assumes wavelength vector is in increasing order.
void getWavelengthSubrange(const vector<double>& wavelength, double wl0, double wlf, unsigned& startindex, unsigned& endindex)
{
    endindex = startindex = 0;
	for (unsigned i=0; i < wavelength.size() && wavelength[i] <= wlf; i++) {
		if (wavelength[i] >= wl0) {
			if (endindex == 0) endindex = startindex = i;
			endindex++;
		}
	}
}

// Generates a spectrum in outputSpectrum along the points in wavelengthVector by using a Gaussian profile to fit telluricLines.
void generateSyntheticTelluricSpectrumUsingGaussianProfile(const Spectrum& telluricLines, const vector <double>& wavelengthVector, vector <double>& outputSpectrum, double resolution)
{
	outputSpectrum.resize(wavelengthVector.size());
	fill(outputSpectrum.begin(), outputSpectrum.end(), 1.0); //Initialize outputSpectrum to uniform 1.0
    for(unsigned i=0; i<wavelengthVector.size(); i++) {
		double gaussianWidth = (wavelengthVector[i]/resolution);
        double wl0 = wavelengthVector[i] - 5*gaussianWidth;
        double wlf = wavelengthVector[i] + 5*gaussianWidth;
        unsigned startindex, endindex;
        getWavelengthSubrange(telluricLines.wavelength, wl0, wlf, startindex, endindex);
        if(args.debug) cout << "operaTelluricWavelengthCorrection: " << i << " gaussianwidth=" << gaussianWidth << " wl0=" << wl0 << " wlf=" << wlf << " nlinesInRange=" << endindex-startindex << endl;

        //Uncomment the following line for old functionality (recreate bug?)
        if(endindex > 0) endindex--;
        for(unsigned j=startindex; j<endindex; j++) {
            double opticaldepth = telluricLines.intensity[j]*exp(-((telluricLines.wavelength[j] - wavelengthVector[i])*(telluricLines.wavelength[j] - wavelengthVector[i])/(2*gaussianWidth*gaussianWidth)))/(sqrt(2*M_PI)*gaussianWidth);
            outputSpectrum[i] *= exp(-opticaldepth);
        }
    }
}

// Generates a spectrum in outputSpectrum along the points in wavelengthVector by using a Voigt profile for telluric lines.
void generateSyntheticTelluricSpectrumUsingVoigtProfile(const Spectrum& telluricLines, const vector <double>& wavelengthVector, vector <double>& outputSpectrum, double resolution)
{
    outputSpectrum.resize(wavelengthVector.size());
    fill(outputSpectrum.begin(), outputSpectrum.end(), 1.0); //Initialize outputSpectrum to uniform 1.0
    for(unsigned i=0; i<wavelengthVector.size(); i++) {
        double gaussianWidth = (wavelengthVector[i]/resolution);
        double gamma = 1.0; // set to one for testing
        double wl0 = wavelengthVector[i] - 5*gaussianWidth;
        double wlf = wavelengthVector[i] + 5*gaussianWidth;
        unsigned startindex, endindex;
        getWavelengthSubrange(telluricLines.wavelength, wl0, wlf, startindex, endindex);
        if(args.debug) cout << "operaTelluricWavelengthCorrection: " << i << " gaussianwidth=" << gaussianWidth << " wl0=" << wl0 << " wlf=" << wlf << " nlinesInRange=" << endindex-startindex << endl;
        
        //Uncomment the following line for old functionality (recreate bug?)
        if(endindex > 0) endindex--;
        for(unsigned j=startindex; j<endindex; j++) {
            double gaussian = exp(-((telluricLines.wavelength[j] - wavelengthVector[i])*(telluricLines.wavelength[j] - wavelengthVector[i])/(2*gaussianWidth*gaussianWidth)))/(sqrt(2*M_PI)*gaussianWidth);
            double lorentz = (1/(M_PI*gamma))*(gamma*gamma)/(gamma*gamma + (telluricLines.wavelength[j] - wavelengthVector[i])*(telluricLines.wavelength[j] - wavelengthVector[i]));
            double opticaldepth = telluricLines.intensity[j]*gaussian*lorentz;
            outputSpectrum[i] *= exp(-opticaldepth);
        }
    }
}


bool calculateRVShiftByXCorr(const Spectrum& telluricLines, const Spectrum& objectSpectrum, double radialVelocityRange, double radialVelocityStep, double threshold, double& maxRV, double& sigRV, double& maxcorr, ofstream& frvcorrdata, ofstream& frvcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double& chisqr)
{
    int jmax = -1;
	maxcorr = 0;
	maxRV = 0;
	sigRV = 0;
	chisqr = 0;
    
    vector <double> crosscorrelation;
    vector <double> crosscorrerror;
    vector <double> dRV;
    // We can try to speed up allocation since we already know what size the vectors should be
    unsigned nDataPoints = (unsigned)ceil(radialVelocityRange/radialVelocityStep);
    crosscorrelation.reserve(nDataPoints);
    crosscorrerror.reserve(nDataPoints);
    dRV.reserve(nDataPoints);

    double xcorrerror = 2e-04; //why this value in particular?
    
    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        Spectrum telluricSpectrum;
        // Initalize telluricSpectrum wavelength with wavelength of objectSpectrum shifted by deltaRV
        for (unsigned i=0; i<objectSpectrum.size(); i++) {
            double DWavelength = deltaRV * objectSpectrum.wavelength[i] / SPEED_OF_LIGHT_KMS;
            telluricSpectrum.wavelength.push_back(objectSpectrum.wavelength[i] + DWavelength);
        }
        // Generate a spectrum in telluricSpectrum along points in wavelength vector using the provided telluricLines
        generateSyntheticTelluricSpectrumUsingGaussianProfile(telluricLines, telluricSpectrum.wavelength, telluricSpectrum.intensity, spectralResolution);
        //generateSyntheticTelluricSpectrumUsingVoigtProfile(telluricLines, telluricSpectrum.wavelength, telluricSpectrum.intensity, spectralResolution);
        
        // Calculate the x-corr between the generated shifted telluric spectrum and the object spectrum
        double xcorr = operaCrossCorrelation(telluricSpectrum.size(), objectSpectrum.intensity_ptr(), telluricSpectrum.intensity_ptr());
        if(args.debug) cout << deltaRV << " " << xcorr << endl;
        
        // Check if this is the highest x-corr we have found so far, but filter out values under threshold
        if(xcorr > threshold && (jmax < 0 || xcorr > maxcorr)) {
            maxcorr = xcorr;
            maxRV = deltaRV;
            sigRV = radialVelocityStep;
            jmax = crosscorrelation.size();
        }
        
        crosscorrelation.push_back(xcorr);
        crosscorrerror.push_back(xcorrerror);
        dRV.push_back(deltaRV);
    }
    
    if (jmax < 0) return false; // Didn't find any x-corr values above threshold
    
    if(useFitToFindMaximum) {
		// Set initial values for our Gaussian using the maximum x-corr.
        double a = crosscorrelation[jmax]; //Initial amplitude
        double x0 = dRV[jmax]; //Initial center
        double sig = radialVelocityRange/4.0; //Initial sigma
        double ea;
        double ex0;
        double esig;
        double fitchisqr;
        
        // Update the initial values and get errors for each along with the fit chi-squared.
        operaMPFitGaussian(crosscorrelation.size(), &(*dRV.begin()), &(*crosscorrelation.begin()), &(*crosscorrerror.begin()), &a, &ea, &x0, &ex0, &sig, &esig, &fitchisqr);
        //operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        
        if(args.debug) {
            cout << a << "+/-" << ea << endl;
            cout << x0 << "+/-" << ex0 << endl;
            cout << sig << "+/-" << esig <<  " fitchisqr=" << fitchisqr << endl;
        }
        
        // For plotting
        if(frvcorrdata.is_open()) {
            for(unsigned j=0; j<crosscorrelation.size(); j++) {
                double x = (double)dRV[j];
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                frvcorrdata << dRV[j] << " " <<  gaussfunc << " " <<  crosscorrelation[j] << " " <<  crosscorrerror[j] << " " << crosscorrelation[j] - gaussfunc << endl;
            }
            frvcorrdata << endl;
        }
        if(frvcorrfitdata.is_open()) {
            frvcorrfitdata  << x0 << " " << ex0 << " " <<  a << " " <<  maxcorr  <<  " " <<  crosscorrerror[jmax]  << " " << maxcorr - a << endl;
        }

        maxcorr = a;
        maxRV = x0;
        sigRV = ex0;
        chisqr = fitchisqr;
    } else {
		// For plotting
        if(frvcorrdata.is_open()) {
            for(unsigned j=0; j<crosscorrelation.size(); j++) {
                frvcorrdata  << dRV[j] << " " <<  crosscorrelation[j] << " " <<  crosscorrelation[j] <<  " " << crosscorrerror[j] << " " << 0.0 << endl;
            }
            frvcorrdata << endl;
        }
        if(frvcorrfitdata.is_open()) {
            frvcorrfitdata  << maxRV << " " <<  sigRV << " " << maxcorr << " " <<  maxcorr  <<  " " <<  crosscorrerror[jmax] <<  " " << 0.0 << endl;
        }
    }
    return true;
}

// Function to match telluric lines
unsigned matchTelluricLines(const Spectrum& telluricLinesFromAtlas, const Spectrum& telluricLinesFromObject, double *telluricMatchedLines, double *objectMatchedLines, double *objectMatchedLinesIntensities, float *radialVelocities, double spectralResolution)
{
    unsigned nmatches = 0;
    
    unsigned initiaAtlasIndex = 0;
    
    for (unsigned j=0; j<telluricLinesFromObject.size(); j++) {
        double linewidth = telluricLinesFromObject.wavelength[j]/spectralResolution;
        
        double diffmin = BIG;
        bool thereIsAMatch = false;
        unsigned matched_jindex = 0;
        unsigned matched_iindex = 0;
        
        for (unsigned i=initiaAtlasIndex; i<telluricLinesFromAtlas.size(); i++) {
            double wl_diff = fabs(telluricLinesFromObject.wavelength[j] - telluricLinesFromAtlas.wavelength[i]);
            
            if(wl_diff < linewidth && wl_diff < diffmin) {
                thereIsAMatch = true;
                matched_jindex = j;
                matched_iindex = i;
                diffmin = wl_diff;
            } else {
                if (telluricLinesFromAtlas.wavelength[i] > telluricLinesFromObject.wavelength[j] + 2*linewidth) {
                    //double radialVelocity = (telluricLinesFromObject.wavelength[matched_jindex] - telluricLinesFromAtlas.wavelength[matched_iindex])*SPEED_OF_LIGHT_KMS/telluricLinesFromAtlas.wavelength[matched_iindex];

                    if (thereIsAMatch) {
                        telluricMatchedLines[nmatches] = telluricLinesFromAtlas.wavelength[matched_iindex];
                        objectMatchedLines[nmatches] = telluricLinesFromObject.wavelength[matched_jindex];
                        objectMatchedLinesIntensities[nmatches] = telluricLinesFromObject.intensity[matched_jindex];
                        radialVelocities[nmatches] = (float)calculateDeltaRadialVelocityInKPS(telluricMatchedLines[nmatches],objectMatchedLines[nmatches],spectralResolution);

                        // cout << nmatches << " " << telluricLinesFromObject.wavelength[matched_jindex] << " " <<  telluricLinesFromAtlas.wavelength[matched_iindex] << " " << radialVelocities[nmatches] << endl;
                        nmatches++;
                    }
                    initiaAtlasIndex = i+1;
                    break;
                }
            }
        }
    }
    return nmatches;
}


unsigned generateHistogramData(unsigned nmatchedLines, double *telluricMatchedLines, float *radialVelocities, double radialVelocityRange, double radialVelocityStep, float *rvVector, float *probDensity, float *wavelengthVector) {
    
    unsigned nTotal = 0;
    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        for (unsigned l=0; l<nmatchedLines; l++) {
            if(radialVelocities[l] > (float)(deltaRV - radialVelocityStep/2) &&  radialVelocities[l] <= float(deltaRV + radialVelocityStep/2)) {
                nTotal++;
            }
        }
    }
    
    unsigned nPointsInHistogram = 0;

    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        
        unsigned npbin = 0;
        float meanwl = 0;
        
        for (unsigned l=0; l<nmatchedLines; l++) {
            if(radialVelocities[l] > (float)(deltaRV - radialVelocityStep/2) &&  radialVelocities[l] <= float(deltaRV + radialVelocityStep/2)) {
                meanwl += (float)telluricMatchedLines[l];
                npbin++;
            }
        }
        
        rvVector[nPointsInHistogram] = (float)deltaRV;
        probDensity[nPointsInHistogram] = (float)npbin/(float)nTotal;
        wavelengthVector[nPointsInHistogram] = meanwl/(float)npbin;
        
        nPointsInHistogram++;
    }
    
    return nPointsInHistogram;
}

