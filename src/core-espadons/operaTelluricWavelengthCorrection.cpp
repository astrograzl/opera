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

#include <fstream>
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralTools.h"
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFit.h"							// for operaFitSplineDouble
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include "core-espadons/operaTelluricWavelengthCorrection.h"

/*! \file operaTelluricWavelengthCorrection.cpp */

using namespace std;

operaArgumentHandler args;

/*
 * the reference Telluric spectrum
 */
static unsigned nPointsInTelluricSpectrum = 0;
static double telluricSpectrumWavelength[MAXNUMBEROFPOINTSINTELLURICSPECTRUM];
static double telluricSpectrumIntensity[MAXNUMBEROFPOINTSINTELLURICSPECTRUM];

/*
 * the reference Telluric spectral lines
 */
static unsigned ntelluriclines = 0;
static int telluricMoleculeNumber[MAXNUMBEROFLINESINTELLURICDATABASE];
static double telluricLinesWavelength[MAXNUMBEROFLINESINTELLURICDATABASE];
static double telluricLinesIntensity[MAXNUMBEROFLINESINTELLURICDATABASE];

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
    
    string xcorrsplotfilename;
    string specplotfilename;
    string xcorrscriptfilename;
    string specscriptfilename;
    string xcorrdatafilename;
    string specdatafilename;
    string xcorrfitdatafilename;
        
    args.AddRequiredArgument("inputWaveFile", inputWaveFile, "input wavelength calibration file (.wcal)");
	args.AddRequiredArgument("inputObjectSpectrum", inputObjectSpectrum, "input object spectrum file (.e or .p)");
    args.AddOptionalArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "", "flat field spectrum ff_");
    args.AddRequiredArgument("outputWaveFile", outputWaveFile, "output wavelength calibration file (.auto or .pauto)");
    args.AddOptionalArgument("telluric_lines", telluric_lines, "", "atlas of telluric lines (HITRAN)");
    args.AddOptionalArgument("telluric_spectrum", telluric_spectrum, "", "spectrum of telluric lines");
    args.AddRequiredArgument("inputWavelengthMaskForTelluric", inputWavelengthMaskForTelluric, "telluric wavelength mask");
    
    args.AddOptionalArgument("spectralResolution", spectralResolution, 80000, "input spectral resolution (wl/dwl) as reference for line detection");
    args.AddOptionalArgument("radialVelocityRange", radialVelocityRange, 10, "radial Velocity Range (in km/s) to scan for first order correction");
    args.AddOptionalArgument("radialVelocityStep", radialVelocityStep, 0.05, "radial Velocity Step step (in km/s) to scan for first order correction");
    args.AddOptionalArgument("XCorrelationThreshold", XCorrelationThreshold, 0.05, "X-correlation lower threshold to consider a match between telluric and object spectra");
    args.AddOptionalArgument("normalizationBinsize", normalizationBinsize, 110, "binsize to normalize input object spectrum");
    args.AddSwitch("useFitToFindMaximum", useFitToFindMaximum, "use gaussian fit to find RV with maximum xcorr");
    
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddSwitch("StarPlusSky", StarPlusSky, "star plus sky mode");
    
    args.AddOptionalArgument("xcorrsplotfilename", xcorrsplotfilename, "", "Output cross-correlation plot eps file name");
    args.AddOptionalArgument("specplotfilename", specplotfilename, "", "Output spectrum plot eps file name");
    args.AddOptionalArgument("xcorrscriptfilename", xcorrscriptfilename, "", "Output cross-correlation gnuplot script file name");
    args.AddOptionalArgument("specscriptfilename", specscriptfilename, "", "Output spectrum gnuplot script file name");
    args.AddOptionalArgument("xcorrdatafilename", xcorrdatafilename, "", "Output cross-correlation data file name");
    args.AddOptionalArgument("specdatafilename", specdatafilename, "", "Output spectrum data file name");
	args.AddOptionalArgument("xcorrfitdatafilename", xcorrfitdatafilename, "", "Output cross-correlation fit data file name");
	
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
		if (telluric_lines.empty() && telluric_spectrum.empty()) {
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
            if(ordernumber != NOTPROVIDED) cout << "operaTelluricWavelengthCorrection: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaTelluricWavelengthCorrection: xcorrsplotfilename = " << xcorrsplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specplotfilename = " << specplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: xcorrscriptfilename = " << xcorrscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specscriptfilename = " << specscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: xcorrdatafilename = " << xcorrdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: xcorrfitdatafilename = " << xcorrfitdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: specdatafilename = " << specdatafilename << endl;
            }
		}
		
		ofstream fxcorrdata;
        ofstream fxcorrfitdata;
        ofstream fspecdata;
        if (!xcorrdatafilename.empty()) fxcorrdata.open(xcorrdatafilename.c_str());
        if (!xcorrfitdatafilename.empty()) fxcorrfitdata.open(xcorrfitdatafilename.c_str());
        if (!specdatafilename.empty()) fspecdata.open(specdatafilename.c_str());
        
		operaSpectralOrderVector spectralOrders(inputObjectSpectrum);
        spectralOrders.ReadSpectralOrders(inputWaveFile); // This merges in the wavelength calibration information

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaTelluricWavelengthCorrection: minorder ="<< minorder << " maxorder=" << maxorder << endl;
    
        // Correct for flat-field
        if (!inputFlatFluxCalibration.empty()) {
            bool starplusskyInvertSkyFiber = false;
            spectralOrders.correctFlatField(inputFlatFluxCalibration, minorder, maxorder, StarPlusSky, starplusskyInvertSkyFiber);
        }
        
        // Get object spectrum within wavelength mask input ranges defined in "inputWavelengthMaskForTelluric"
        double *objectSpectrum = new double[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
        double *objectSpectrumVariance = new double[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
        double *wavelength = new double [MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
        unsigned nelem = spectralOrders.getSpectrumWithinTelluricMask(inputWavelengthMaskForTelluric, minorder, maxorder, TRUE, normalizationBinsize, wavelength, objectSpectrum,objectSpectrumVariance);
#ifdef PRINT_DEBUG
        for(unsigned i=0; i<nelem; i++) {
            cout << wavelength[i] << " " << objectSpectrum[i] << " " << objectSpectrumVariance[i] << endl;
        }
#endif

        // Define wavelength full range
        double wl0 = wavelength[0] - radialVelocityRange*wavelength[0]/SPEED_OF_LIGHT_KMS;
        double wlf = wavelength[nelem-1] + radialVelocityRange*wavelength[nelem-1]/SPEED_OF_LIGHT_KMS;
        if(args.verbose) cout << "operaTelluricWavelengthCorrection: wl0=" <<  wl0 << " wlf=" << wlf << endl;

		// Read telluric lines database lambda vs. intensity
        double maxatlasflux = -BIG;
        if (!telluric_lines.empty()) {
            if (args.debug) {
                cout << "operaTelluricWavelengthCorrection: reading telluric lines database " << telluric_lines << endl;
            }
            ntelluriclines = readTelluricLines(telluric_lines,telluricMoleculeNumber,telluricLinesWavelength,telluricLinesIntensity);
#ifdef PRINT_DEBUG
            for(unsigned line=0;line<ntelluriclines;line++) {
                cout << line << " " << telluricMoleculeNumber[line] << " " << telluricLinesWavelength[line] << " " << telluricLinesIntensity[line] << endl;
            }
#endif
            // Find line with maximum absorption for normalization
            for (unsigned i=0; i<ntelluriclines; i++) {
                if(1 - telluricLinesIntensity[i] > maxatlasflux) maxatlasflux = 1 - telluricLinesIntensity[i];
            }
        } else {
            maxatlasflux = 1;
        }
        
        // Read Telluric reference spectrum lambda vs. intensity
        if (!telluric_spectrum.empty()) {
            nPointsInTelluricSpectrum = readTelluricSpectrum(telluric_spectrum, telluricSpectrumWavelength, telluricSpectrumIntensity);
        }
        
        // Calculate radial velocity shift by cross-correlation between observed spectrum and telluric reference.
        if (args.verbose) {
            cout << "operaTelluricWavelengthCorrection: calculating cross-correlation for radialVelocityRange=" << radialVelocityRange << " km/s and radialVelocityStep=" << radialVelocityStep << " km/s" << endl;
        }
        double rvshift = 0;
        double rvshifterror = 0;
        double maxcorr = 0;
        double chisqr = 0;
        bool validXCorrelation = calculateRVShiftByXCorr(nelem, wavelength, objectSpectrum, radialVelocityRange, radialVelocityStep, XCorrelationThreshold, &rvshift, &rvshifterror, &maxcorr, fxcorrdata, fxcorrfitdata, spectralResolution, useFitToFindMaximum, &chisqr);
        if(validXCorrelation) {
            if (args.verbose) cout << "\noperaTelluricWavelengthCorrection: Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s, maxXCorr=" << maxcorr << ", chisqr=" << chisqr << "\n" << endl;
        } else {
            rvshift = 0;
            rvshifterror = 0;
            maxcorr = 0;
        }
        
        //Apr 15, 2015 CU - Output rvel shift instead of wcal file
        spectralOrders.setTelluricRadialVelocityCorrection(rvshift);
        spectralOrders.WriteSpectralOrders(outputWaveFile, Tell);
        
        if(fspecdata.is_open()) {
            /*
             * Generate telluric synthetic spectrum using both HITRAN lines and input telluric spectrum:
             */
            
            // First use hitran lines to generate synthetic spectrum
            double *hitranTelluricSpectrum = new double[nelem];
            if (!telluric_lines.empty()) generateSyntheticTelluricSpectrumUsingGaussianProfile(nelem,wavelength,hitranTelluricSpectrum,spectralResolution);
            
            // Then read telluric spectrum from input spectrum file.
            double *KPNOTelluricSpectrum = new double[nelem];
            if (!telluric_spectrum.empty()) {
                double *wl,*transmission_p;
                unsigned npointsInShortTelluricSpectrum = getTelluricSpectrumRange(wl0,wlf,&wl,&transmission_p);
                double *transmission = new double [npointsInShortTelluricSpectrum];
                // Normalize telluric reference
                for (unsigned i=0; i<npointsInShortTelluricSpectrum; i++) {
                    transmission[i]= 1 - (1-transmission_p[i])/maxatlasflux;
                }
                // Spline to match same sampling as in observed spectrum
                operaFitSplineDouble(npointsInShortTelluricSpectrum,wl,transmission,nelem,wavelength,KPNOTelluricSpectrum);
            }
            
            // Print spectral data to file
            for(unsigned i=0; i<nelem; i++) fspecdata << wavelength[i] << " " << objectSpectrum[i] << " " << hitranTelluricSpectrum[i] << " " << KPNOTelluricSpectrum[i] << endl;
        }
        
        
        // Telluric wavelength correction orders info plot:
        if (fxcorrdata.is_open() && fxcorrfitdata.is_open()) {
            fxcorrdata.close();
            fxcorrfitdata.close();
            if (!xcorrscriptfilename.empty()) {
                GenerateTelluricXCorrelationPlot(xcorrscriptfilename, xcorrsplotfilename, xcorrdatafilename, xcorrfitdatafilename);
            }
        }
        
        // Spectrum plot: plot observed and reference telluric spectra.
        if (fspecdata.is_open()) {
            fspecdata.close();
            if (!specscriptfilename.empty()) {
                GenerateTelluricSpecPlot(specscriptfilename, specplotfilename, specdatafilename);
            }
        }

        //delete[] objectSpectrum;
        //delete[] wavelength;
        //delete[] objectSpectrumVariance;

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

/*
 * Generate multiple plot containing statistical info about telluric wavelength correction
 */
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


/*
 * Generate 2D plot for spectra of atlas + comparison + identified lines
 */
void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename) {
    
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

/*
 * Read the entire set of telluric lines in HITRAN database
 */
unsigned readTelluricLines(string telluric_database_file, int *telluricMoleculeNumber, double *telluricLinesWavelength, double *telluricLinesIntensity) {
   	igzstream astream;
	string dataline;
    int tmpnumber = -1;
	double tmpwn = -1.0;
	float tmpi = -1.0;
	unsigned line = 0;
	char *buff = new char[MAXLENGTHOFLINEINTELLURICDATABASE];
   
    int *tmp_telluricMoleculeNumber = new int[MAXNUMBEROFLINESINTELLURICDATABASE];
    double *tmp_telluricLinesWavelength = new double[MAXNUMBEROFLINESINTELLURICDATABASE];
    double *tmp_telluricLinesIntensity = new double[MAXNUMBEROFLINESINTELLURICDATABASE];
    
	astream.open(telluric_database_file.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
                    sscanf(dataline.c_str(), "%d %lf %G %[^\n]", &tmpnumber, &tmpwn, &tmpi, buff);
                    
                    //tmpi = tmpi*pow(10,tmpexp);
                    double wavelength_in_nm = 1e7/tmpwn;
                    
                    tmp_telluricMoleculeNumber[line] = tmpnumber;
                    
                    tmp_telluricLinesWavelength[line] = convertVacuumToAirWavelength(wavelength_in_nm*10)/10;
                    
                    double NoverV = TYPICAL_PRESSURE_AT_MAUNAKEA/(TYPICAL_TEMPERATURE_AT_MAUNAKEA*k_BOLTZMANN_CONSTANT);
                    tmp_telluricLinesIntensity[line] = ((double)tmpi/(NoverV*1e-6))/TYPICAL_ATMOSPHERE_PATH_LENGTH;

//                    printf("%u %d %lf %lf\n",line,tmpnumber,wavelength_in_nm,(double)tmpi*ONEMOLE);
                    
                    line++;
 				}	// skip comments
			}	// if strlen
		} // while (astream.good())
		line--;
   		if (line > 0) {
			if (args.verbose) {
				printf("          [Telluric] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", line, tmp_telluricLinesWavelength[0], tmp_telluricLinesWavelength[line/2], tmp_telluricLinesWavelength[line-1]);
			}
		} else {
			printf("          [Telluric] no lines found in telluric database.\n");
		}
		astream.close();
	}	// if (astream.open()
    delete[] buff;
    
    unsigned np = line;
    for(unsigned i=0; i<line; i++) {
        np--;
        telluricMoleculeNumber[i] = tmp_telluricMoleculeNumber[np];
        telluricLinesWavelength[i] = tmp_telluricLinesWavelength[np];
        telluricLinesIntensity[i] = tmp_telluricLinesIntensity[np];
//        cout << "operaTelluricWavelengthCorrection: " << telluricMoleculeNumber[i] << " " << telluricLinesWavelength[i] << " " << telluricLinesIntensity[i] << endl;
    }

    delete[] tmp_telluricMoleculeNumber;
    delete[] tmp_telluricLinesWavelength;
    delete[] tmp_telluricLinesIntensity;
    
	return line;
}

/*
 * get a subset of the telluric lines for this order only, between wl0 and wlf
 */
unsigned getTelluricLinesRange(double wl0, double wlf, double **wl, double **intensity) {
    
    unsigned firstline = 0;
	unsigned line = 0;
	
	for (line=0; line<ntelluriclines; line++) {
		if (telluricLinesWavelength[line] >= wl0) {
			if (firstline == 0) {
				*intensity = &telluricLinesIntensity[line];
				*wl = &telluricLinesWavelength[line];
				firstline = line;
			}
			if (telluricLinesWavelength[line] > wlf)
				break;
		}
	}
	if (line) line--;
	if (line > firstline) return (line-firstline);
	return 0;
}


void generateSyntheticTelluricSpectrumUsingGaussianProfile(unsigned np, double *wavelengthVector, double *ouputSpectrum, double resolution) {
    
    for(unsigned i=0;i<np;i++) {    
        ouputSpectrum[i] = 1.0;
    }
    
    double *wl, *intensity;
    
    for(unsigned i=0;i<np;i++) {
        double gaussianWidth = (wavelengthVector[i]/resolution);

        unsigned nlinesInRange = getTelluricLinesRange(wavelengthVector[i] - 5*gaussianWidth,wavelengthVector[i] + 5*gaussianWidth,&wl,&intensity);
        
        if(args.debug) {
            cout << "operaTelluricWavelengthCorrection: " << i <<
            " gaussianwidth=" << gaussianWidth <<
            " wl0=" << wavelengthVector[i] - 5*gaussianWidth <<
            " wlf=" << wavelengthVector[i] + 5*gaussianWidth <<
            " nlinesInRange=" << nlinesInRange << endl;
        }
        for(unsigned j=0; j<nlinesInRange; j++) {
            double opticaldepth = intensity[j]*exp(-((wl[j] - wavelengthVector[i])*(wl[j] - wavelengthVector[i])/(2*gaussianWidth*gaussianWidth)))/(sqrt(2*M_PI)*gaussianWidth);
            ouputSpectrum[i] *= exp(-opticaldepth);
        }
    }
}

/*
 * Read the the full atmospheric transmission spectrum
 */
unsigned readTelluricSpectrum(string telluric_spectrum, double *telluricSpectrumWavelength, double *telluricSpectrumIntensity) {
	igzstream astream;
	string dataline;
    
	double tmpwl = -1.0;
	double tmpi = -1.0;
	unsigned np = 0;
	
	astream.open(telluric_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                    
                    telluricSpectrumWavelength[np] = tmpwl;
                    telluricSpectrumIntensity[np] = tmpi;
                    np++;
                }	// skip comments
            }
		} // while (astream.good())
		
		if (np > 0) {
			if (args.verbose) {
				printf("          [Telluric] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, telluricSpectrumWavelength[0], telluricSpectrumWavelength[np/2], telluricSpectrumWavelength[np-1]);
			}
		} else {
			printf("          [Telluric] no points found in telluric spectrum.\n");
		}
		astream.close();
	}	// if (astream.open()
	return np;
}

/*
 * get a subset of the telluric lines for this order only, between wl0 and wlf
 */
unsigned getTelluricSpectrumRange(double wl0, double wlf, double **wl, double **intensity) {
    
    unsigned firstline = 0;
	unsigned line = 0;
	
	for (line=0; line<nPointsInTelluricSpectrum; line++) {
		if (telluricSpectrumWavelength[line] >= wl0) {
			if (firstline == 0) {
				*intensity = &telluricSpectrumIntensity[line];
				*wl = &telluricSpectrumWavelength[line];
				firstline = line;
			}
			if (telluricSpectrumWavelength[line] > wlf)
				break;
		}
	}
	if (line) line--;
	if (line > firstline) return (line-firstline);
	return 0;
}

unsigned matchTelluricReferencewithObjectLines(double acceptableMismatch,double lineSigma, operaSpectralLines *telluricLines, operaSpectralLines *objectLines, Polynomial *wlcorrection, unsigned order, double wl_central, ofstream *flinesdata) {
    /*
     * vectors for uncalibrated spectral lines
     */
	if (objectLines->getnLines() == 0) {
		throw operaException("operaTelluricWavelengthCorrection: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (telluricLines->getnLines() == 0) {
		throw operaException("operaTelluricWavelengthCorrection: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    double *objectlinecenter = new double[objectLines->getnLines()];
    double *objectlinecenterError = new double[objectLines->getnLines()];
    double *objectlineflux = new double[objectLines->getnLines()];
    double *objectlinesigma = new double[objectLines->getnLines()];
    
    double *objectMatchedData = new double[objectLines->getnLines()];
    unsigned *objectMatchedindex = new unsigned[objectLines->getnLines()];
    double *telluricMatchedData = new double[telluricLines->getnLines()];
    unsigned *telluricMatchedindex = new unsigned[telluricLines->getnLines()];
    
    /*
     * vectors for telluric spectral lines
     */
    double *telluricLineswl = new double[telluricLines->getnLines()];
    double *telluricLineswlError = new double[telluricLines->getnLines()];
    double *telluricLinesflux = new double[telluricLines->getnLines()];
    
    unsigned line = 0;
    
    for(unsigned feature=0;feature<objectLines->getNFeatures();feature++) {
        operaSpectralFeature *currentFeature = objectLines->getSpectralFeature(feature);
        double *center = currentFeature->getGaussianFit()->getCenterVector();
        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
        for(unsigned l=0; l<currentFeature->getnLines(); l++) {
            objectlinecenter[line] = center[l];
            objectlinecenterError[line] = centerError[l];
            objectlineflux[line] = amplitude[l];
            objectlinesigma[line] = sigma[l];
            if(args.debug)
                cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << endl;
            line++;
        }
    }
    line = 0;
    
    for(unsigned feature=0;feature<telluricLines->getNFeatures();feature++) {
        operaSpectralFeature *currentFeature = telluricLines->getSpectralFeature(feature);
        double *center = currentFeature->getGaussianFit()->getCenterVector();
        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
        for(unsigned l=0; l<currentFeature->getnLines(); l++) {
            telluricLineswl[line] = center[l];
            telluricLineswlError[line] = centerError[l];
            telluricLinesflux[line] = amplitude[l];
            if(args.debug)
                cout << center[l] <<  " " << amplitude[l] << " " << sigma[l] << endl;
            line++;
        }
    }
    
    /*
     * Below it identifies and select the set of lines that match both comparison and atlas.
     * The criteria for matching is that the difference between centers must be < acceptableMismatch x sigma
     */
    unsigned nmatch = 0;
    unsigned nextfirstline = 0;
    
    double acceptMismatchInwlUnits = acceptableMismatch*lineSigma;
    
    for (unsigned i=0; i<objectLines->getnLines(); i++) {
        
        unsigned bestAtlasMatchIndex = 0;
        double mindifference = acceptMismatchInwlUnits;
        
        for(unsigned l=nextfirstline;l<telluricLines->getnLines();l++) {
            
            double difference = fabs(objectlinecenter[i] - telluricLineswl[l]);
#ifdef PRINT_DEBUG
            cout << "operaTelluricWavelengthCorrection: " << "mindiff=" << mindifference << " diff=" << difference << " object[" << i << "]=" << objectlinecenter[i] << " telluric[" << l << "]=" << telluricLineswl[l] << endl;
#endif
            if(objectlinecenter[i] > telluricLineswl[l]  && difference < mindifference) {
                mindifference = difference;
                bestAtlasMatchIndex = l;
            } else if (objectlinecenter[i] <= telluricLineswl[l] && difference < mindifference) {
                objectMatchedData[nmatch] = objectlinecenter[i];
                telluricMatchedData[nmatch] = telluricLineswl[l];
                telluricMatchedindex[nmatch] = l;
                objectMatchedindex[nmatch] = i;
                nextfirstline = l+1;
                nmatch++;
                break;
            } else if (objectlinecenter[i] <= telluricLineswl[l]  && difference > mindifference) {
                if(bestAtlasMatchIndex) {
                    objectMatchedData[nmatch] = objectlinecenter[i];
                    telluricMatchedData[nmatch] = telluricLineswl[bestAtlasMatchIndex];
                    telluricMatchedindex[nmatch] = bestAtlasMatchIndex;
                    objectMatchedindex[nmatch] = i;
                    nextfirstline = bestAtlasMatchIndex+1;
                    nmatch++;
                }
                break;
            }
            if(l==telluricLines->getnLines()-1) {
                if(bestAtlasMatchIndex) {
                    objectMatchedData[nmatch] = objectlinecenter[i];
                    telluricMatchedData[nmatch] = telluricLineswl[bestAtlasMatchIndex];
                    telluricMatchedindex[nmatch] = bestAtlasMatchIndex;
                    objectMatchedindex[nmatch] = i;
                    nextfirstline = bestAtlasMatchIndex+1;
                    nmatch++;
                }
            }
        }
    }
    
	if (nmatch > 0) {
        if(flinesdata != NULL){
            for(unsigned index=0; index<nmatch;index++){
                *flinesdata << order << " " << objectMatchedData[index] << " " << telluricMatchedData[index] << " " << objectlineflux[objectMatchedindex[index]] << " " << objectlinesigma[objectMatchedindex[index]] << " " << wl_central << endl;
            }
            *flinesdata << endl;
        }
        if (wlcorrection) {
            int npar = wlcorrection->getOrderOfPolynomial();
            double *par = (double *)wlcorrection->getVector();
            double chisqr;
            
            operaLMFitPolynomial(nmatch, telluricMatchedData, objectMatchedData, npar, par, &chisqr);

            wlcorrection->setChisqr(chisqr);
        }
   	}
    
    delete[] objectMatchedData;
    delete[] objectMatchedindex;
    delete[] telluricMatchedData;
    delete[] telluricMatchedindex;
    
    delete[] objectlinecenter;
    delete[] objectlinecenterError;
    delete[] objectlineflux;
    delete[] objectlinesigma;
    delete[] telluricLineswl;
    delete[] telluricLineswlError;
    delete[] telluricLinesflux;
    
    return nmatch;
}


bool calculateRVShiftByXCorr(unsigned nelem, double *wavelength, double *objectSpectrum, double radialVelocityRange, double radialVelocityStep, double threshold, double *maxRV, double *sigRV, double *maxcorr, ofstream& fxcorrdata, ofstream& fxcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double *chisqr) {
    bool status = true;
    
    double *telluricWavelength = new double[nelem];
    double *telluricSpectrum = new double[nelem];
    
    unsigned nDataPoints = (unsigned)ceil(radialVelocityRange/radialVelocityStep);
    
    double firstRV = -radialVelocityRange/2.0;

    double deltaRV = firstRV;
    
    unsigned jmax = 0;
    *maxcorr = -BIG;
    *maxRV = 0;
    
    double *crosscorrelation = new double [nDataPoints];
    double *crosscorrerror = new double [nDataPoints];
    double *dRV = new double [nDataPoints];
    double xcorrerror = 2e-04;
    
    for(unsigned j=0; j<nDataPoints;j++) {
        
        for (unsigned i=0; i<nelem; i++) {
            double DWavelength = deltaRV * wavelength[i] / SPEED_OF_LIGHT_KMS;
            telluricWavelength[i] = wavelength[i] + DWavelength;
        }
        
        //operaFitSplineDouble(nelem,telluricWavelength,telluricSpectrum,nelem,wavelength,telluricNewSpectrum);
        
        generateSyntheticTelluricSpectrumUsingGaussianProfile(nelem,telluricWavelength,telluricSpectrum,spectralResolution);
        
        crosscorrelation[j] = operaCrossCorrelation(nelem,objectSpectrum,telluricSpectrum);
        crosscorrerror[j] = xcorrerror;
        dRV[j] = deltaRV;
        
        if(args.debug) {
            cout << dRV[j] << " " << crosscorrelation[j] << endl;
        }
        // Test (I) :
        if(crosscorrelation[j] > *maxcorr && crosscorrelation[j] > threshold) {
            *maxcorr = crosscorrelation[j];
            *maxRV = deltaRV;
            *sigRV = radialVelocityStep;
            jmax = j;
        }
        
        deltaRV+=radialVelocityStep;
    }
    
    if (*maxRV == firstRV) {
        *maxcorr = 0;
        *maxRV = 0;
        *sigRV = 0;
        *chisqr = 0;
        status = false;
    }
    
    if(useFitToFindMaximum && status == TRUE) {
        
        double a=(double)crosscorrelation[jmax];
        double x0=(double)dRV[jmax];
        double sig=(double)radialVelocityRange/4;
        double ea;
        double ex0;
        double esig;
        double fitchisqr;
        
        //            operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        operaMPFitGaussian(nDataPoints, dRV, crosscorrelation, crosscorrerror, &a, &ea, &x0, &ex0, &sig, &esig, &fitchisqr);
        
        if(args.debug) {
            cout << a << "+/-" << ea << endl;
            cout << x0 << "+/-" << ex0 << endl;
            cout << sig << "+/-" << esig <<  " fitchisqr=" << fitchisqr << endl;
        }
        
        // Below is for plotting
        if(fxcorrdata.is_open()) {
            for(unsigned j=0; j<nDataPoints;j++) {
                double x = (double)dRV[j];
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                
                fxcorrdata << dRV[j] << " " <<  gaussfunc << " " <<  crosscorrelation[j] << " " <<  crosscorrerror[j] << " " << crosscorrelation[j] - gaussfunc << endl;
            }
            fxcorrdata << endl;
        }
        // Below is for plotting
        if(fxcorrfitdata.is_open()) {
            fxcorrfitdata  << x0 << " " << ex0 << " " <<  a << " " <<  *maxcorr  <<  " " <<  crosscorrerror[jmax]  << " " << *maxcorr - a << endl;
        }

        *maxcorr = a;
        *maxRV = x0;
        *sigRV = ex0;
        *chisqr = fitchisqr;
        
    } else if (!useFitToFindMaximum && status == TRUE) {

        // Below is for plotting
        if(fxcorrdata.is_open()) {
            for(unsigned j=0; j<nDataPoints;j++) {
                fxcorrdata  << dRV[j] << " " <<  crosscorrelation[j] << " " <<  crosscorrelation[j] <<  " " <<  crosscorrerror[j] << " " <<  0.0  << endl;
            }
            fxcorrdata << endl;
        }
        
        // Below is for plotting
        if(fxcorrfitdata.is_open()) {
            fxcorrfitdata  << *maxRV << " " <<  *sigRV << " " << *maxcorr << " " <<  *maxcorr  <<  " " <<  crosscorrerror[jmax] <<  " " << 0.0 << endl;
        }
         *chisqr = 0;
    }
    
    delete[] crosscorrelation;
    delete[] crosscorrerror;
    delete[] dRV;
    
    delete[] telluricWavelength;
    delete[] telluricSpectrum;
    
    return status;
}
