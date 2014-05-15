#ifndef OPERASPECTRALORDER_H
#define OPERASPECTRALORDER_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrder
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
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

#include "operaError.h"

#include "libraries/operaLibCommon.h"           // for operaSpectralOrder_t
#include "libraries/Polynomial.h"				// for Polynomial
#include "libraries/operaGeometry.h"			// for operaGeometry
#include "libraries/operaPolarimetry.h"			// for Polarimetry
#include "libraries/operaWavelength.h"			// for operaWavelength
#include "libraries/operaSpectralElements.h"	// for operaSpectralElements
#include "libraries/operaInstrumentProfile.h"	// for operaInstrumentProfile
#include "libraries/operaSpectralLines.h"       // for operaSpectralLines
#include "libraries/operaFITSImage.h"			// for Matrix
#include "libraries/operaExtractionAperture.h"  // for operaExtractionAperture
#include "libraries/operaSpectralEnergyDistribution.h" // for operaSpectralEnergyDistribution
#include "libraries/GainBiasNoise.h" // for GainBiasNoise

using namespace std;

/*! 
 * \sa class operaSpectralOrder
 * \brief operaSpectralOrder
 * \details A spectral order (SO) consists of a data set containing 
 * \details the information concerned with a full spectral order.
 * \return none
 * \file operaSpectralOrder.h
 * \ingroup libraries
 */
class operaSpectralOrder {
	
private:
	
	unsigned orderNumber;
	
	operaSpectralOrder_t SpectrumType;				// what type of spectral order data is stored
	
	double SNR;										// SNR at wlc for this order
    
	doubleValue_t tiltInDegrees;					// tiltInDegrees    
	
	operaSpectralElements *SpectralElements;		// pointer to the operaSpectralElements
	
    operaSpectralElements *SkyElements;             // pointer to the sky elements        
    
	operaGeometry *Geometry;						// pointer to the operaGeometry class instance
	
	operaWavelength *Wavelength;					// pointer to the operaWavelength class instance
	
	operaInstrumentProfile *InstrumentProfile;		// pointer to the operaInstrumentProfile 
    
	operaSpectralLines *SpectralLines;              // pointer to the operaSpectralLines     
	
    operaPolarimetry *Polarimetry;					// pointer to the polarimetry    
    
	operaSpectralEnergyDistribution *SpectralEnergyDistribution;	// pointer to the operaSpectralEnergyDistribution class instance    
    
    unsigned numberOfBeams;
    
    operaSpectralElements *BeamElements[MAXNUMBEROFBEAMS]; // pointer to spectral elements for beams
    
    operaInstrumentProfile *BeamProfiles[MAXNUMBEROFBEAMS];    // pointer to instrument profiles for beams         
    
    operaSpectralElements *BackgroundElements[LEFTANDRIGHT]; // pointer to spectral elements for background 
    
	operaExtractionAperture *ExtractionApertures[MAXNUMBEROFBEAMS];    // pointer to the apertures for extraction
	
    operaExtractionAperture *BackgroundApertures[LEFTANDRIGHT];    // pointer to the apertures for background    
	
	operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];	// pointer to the operaSpectralEnergyDistribution class instance        
    
	bool hasSpectralElements;						// true if we have information of this type about this order
	bool hasSkyElements;    
	bool hasGeometry;
	bool hasWavelength;
	bool hasInstrumentProfile;
	bool hasSpectralLines;
	bool hasExtractionApertures;
	bool hasPolarimetry;    
	bool hasSNR;
	bool hasCenterSNROnly;
	bool hasSpectralEnergyDistribution;    
	
public:
	
	/*
	 * Constructors
	 */
	operaSpectralOrder();
	
	operaSpectralOrder(unsigned order) ;
	
	operaSpectralOrder(unsigned order, unsigned maxDatapoints, unsigned maxValues, unsigned maxElements, operaSpectralOrder_t format);
	
	/*
	 * Destructor
	 */
	~operaSpectralOrder();
	
	/*
	 * Common Methods
	 */
	
	void deleteAll();
	
	void deleteInstrumentProfile(void);	
	
	void deleteBeamProfiles(void);
    
	void deleteApertures(void);    
    
	void deleteBeamsAndBackgroundElements(void);    
	
	void createBeamsAndBackgrounds(unsigned nElements, unsigned nBeams, operaSpectralOrder_t format);
    
    void deletePolarimetry(void);
    
    void createPolarimetry(unsigned nElements);
    
	unsigned getorder(void);
	
	void sethasSpectralElements(bool HasSpectralElement);
	
	void sethasSkyElements(bool HasSkyElements);
	
	void sethasGeometry(bool HasGeometry);
	
	void sethasWavelength(bool HasWavelength);
	
	void sethasInstrumentProfile(bool hasInstrumentProfile);
	
	void sethasSpectralLines(bool HasSpectralLines);	
    
	void sethasExtractionApertures(bool HasSpectralLines); 
    
	void sethasPolarimetry(bool HasPolarimetry);
	
	void sethasSNR(bool HasSNR);
    
	void sethasCenterSNROnly(bool HasCenterSNROnly);
    
	void sethasSpectralEnergyDistribution(bool HasSpectralEnergyDistribution);     
	
	bool gethasSpectralElements(void);
	
	bool gethasSkyElements(void);	
    
	bool gethasGeometry(void);
	
	bool gethasWavelength(void);
	
	bool gethasInstrumentProfile(void);
	
	bool gethasSpectralLines(void);    
    
	bool gethasExtractionApertures(void);     
	
	bool gethasPolarimetry(void);     
    
	bool gethasSNR(void);     
    
	bool gethasCenterSNROnly(void);     
    
	bool gethasSpectralEnergyDistribution(void);       
	
	void createGeometry(unsigned maxdatapoints, unsigned maxValues);
	
	void createWavelength(unsigned maxnumberofcoefficients);
	
	void createSpectralElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType, bool extended = false);
	
    void createSkyElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType);    
#if 0
    -- Jan 2013 DT deprecated -- bad interface
	void createSpectralEnergyDistribution(void);
#endif	
	operaSpectralElements *getSpectralElements(void);
	
	operaSpectralElements *getSkyElements(void);	
    
	operaGeometry *getGeometry(void);
	
	operaWavelength *getWavelength(void);
	
	operaInstrumentProfile *getInstrumentProfile(void);
	
	void setSpectralLines(operaSpectralLines *spectralLines);
	
	operaSpectralLines *getSpectralLines(void);
	
    operaPolarimetry *getPolarimetry(void);
    
    operaSpectralEnergyDistribution *getSpectralEnergyDistribution(void);
	
    operaSpectralElements *getBeamElements(unsigned beam);
	
    operaInstrumentProfile *getBeamProfiles(unsigned beam);    
    
    operaSpectralEnergyDistribution *getBeamSED(unsigned beam);
	
    operaSpectralElements *getBackgroundElements(unsigned LeftOrRight);
    
	operaExtractionAperture *getExtractionApertures(unsigned beam);
	
	operaExtractionAperture *calculateMainApertureFromExtractionBeams(bool useIP);
    
    operaExtractionAperture *getBackgroundApertures(unsigned LeftOrRight);
    
    void setBeamElements(unsigned beam, operaSpectralElements *beamElements);
    
    void setBeamProfiles(unsigned beam, operaInstrumentProfile *beamProfiles);    
    
    void setBackgroundElements(unsigned LeftOrRight, operaSpectralElements *backgroundElements);
    
	void setExtractionApertures(unsigned beam, operaExtractionAperture *extractionApertures);
	
    void setBackgroundApertures(unsigned LeftOrRight, operaExtractionAperture *backgroundApertures);
	
	double getCenterSNR(void);
	
	void setCenterSNR(double Snr);
    
	unsigned getnumberOfBeams(void);   
    
	void setnumberOfBeams(unsigned NumberOfBeams);    
    
	doubleValue_t getTiltInDegrees(void);   
    
	void setTiltInDegrees(doubleValue_t TiltInDegrees); 
    
    void setTiltInDegrees(double tilt, double error);
    
    double getTiltInDegreesValue(void);
    
    double getTiltInDegreesError(void);
	
    double *getSNRVector(void);    
	
	void calculateSNR(void);
	
	float getCentralSmoothedSNR(int upperlowerbound);
	
	float getPeakSmoothedSNR(int upperlowerbound);

	float getLECompatibleSNR(void);

	void NormalizeFlat(operaFITSImage &flatMatrix, operaFITSImage &outputMatrix, unsigned nx, unsigned ny, unsigned binsize);
	
	void extractRawSum(operaFITSImage &inputImage, ofstream &sout);
	
	void extractRawSum(operaFITSImage &inputImage, operaFITSProduct &outputSpectrum);
	
	void extractRawSum(operaFITSImage &inputImage, float noise, float gain);
	
    void measureInstrumentProfileAlongRows(operaFITSImage &masterFlatImage, unsigned binsize, unsigned sampleElementForPlot, ostream *pout);
	
    void measureInstrumentProfileAlongRowsInto2DWithGaussian(operaFITSImage &masterFlatImage, operaFITSImage &badpix, unsigned binsize, float gaussSig, float tiltInDegrees, bool witherrors, unsigned sampleElementForPlot, ostream *pout);
	
    void measureInstrumentProfileAlongRowsInto2DWithGaussian(operaFITSImage &masterFlatImage, unsigned binsize,float gaussSig, float tiltInDegrees, ostream *pout);
	
	void CalculateWavelengthSolution(void);
	
	void setSpectralElementsByHeight(double Height);

    void setSpectralElementsByStitchingApertures(double effectiveApertureFraction);
    
	void setInstrumentProfileVector(unsigned IPxsize, unsigned IPxsampling, unsigned IPysize, unsigned IPysampling, unsigned NDataPoints);
	
    void setSpectralLines(operaFITSImage &masterCompImage, operaFITSImage &badpix, operaFITSImage &bias, float noise, float gain, float ReferenceLineWidth,float DetectionThreshold, float LocalMaxFilterWidth, float MinPeakDepth);        
	
    void measureInstrumentProfile(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout);
	
    void measureInstrumentProfileWithBinning(operaFITSImage &masterCompImage, operaFITSImage &badpix, double binsize, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout);
    
    void measureInstrumentProfileUsingMedian(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout);
    
    void measureInstrumentProfileUsingWeightedMean(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout);
    
    void recenterOrderPosition(void);
    
    void setApertureElements(operaSpectralOrder_t format);
    
    void extractRawSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, ostream *pout);
	
    void extractStandardSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, ostream *pout);
    
    void extractStandardSpectrumNoBackground(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, ostream *pout);
	
    void extractOptimalSpectrum(operaFITSImage &objectImage, operaFITSImage &flatImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, unsigned sigmaClip, unsigned iterations, bool onTargetProfile, bool usePolynomialFit, bool removeBackground, bool verbose, bool calculateXCorrelation, ostream *pout);
	
    void measureBeamSpatialProfiles(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, bool usePolynomialFit);
    
    void measureOptimalSpectrum(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned sigmaClip, ostream *pout);
        
    void refineBeamSpatialProfiles(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix,GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned NumberofElementsToBin, unsigned sigmaClip, bool usePolynomialFit);

	void calculateXCorrBetweenIPandImage(operaFITSImage &Image, operaFITSImage &badpix, ostream *pout);
	
	operaSpectralOrder_t getSpectrumType(void);
	
	void setSpectrumType(operaSpectralOrder_t format);
	
    void printBeamSpectrum(ostream *pout);
	
	void printBeamSpectrum(string addFirstColumnEntry, ostream *pout);
	
	/*
	 * Normalization/Flux Calibration...
	 */
	void applyNormalization(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, ostream *poutspec, ostream *poutcontinuum, bool overwriteUncalFlux, unsigned numberOfprintouts);
	
    void applyNormalizationForEmissionSpectrum(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, ostream *poutspec, ostream *poutcontinuum, bool overwriteUncalFlux, unsigned numberOfprintouts);

    void applyNormalizationFromExistingContinuum(ostream *poutspec, ostream *poutcontinuum, bool overwriteUncalFlux, bool normalizeBeams, unsigned numberOfprintouts);

    void normalizeSpectrum(operaFluxVector &uncalibratedFlux, operaFluxVector &normalizedFlux, operaFluxVector &outputContinuum, operaSpectralEnergyDistribution &spectralEnergyDistribution, unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial);
	
    void measureContinuum(operaFluxVector &uncalibratedFlux,operaFluxVector &outputContinuum, operaSpectralEnergyDistribution &spectralEnergyDistribution, unsigned binsize, unsigned nsigcut, unsigned orderOfPolynomial, bool usePolynomial);    
	
    void deleteSpectralEnergyDistribution(void);
    
    void createSpectralEnergyDistribution(unsigned binsize);
    
    void createSpectralEnergyDistributionElements(unsigned nElements);
    
    void calculateContinuum(unsigned binsize, unsigned nsigcut, ostream *poutspec, ostream *poutcontinuum);
    
    void calculateFluxCalibrationFromExistingContinuum(unsigned nPointsInReference,double *refwl,double *refflux,double referenceFluxForNormalization,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS],double uncalibratedContinuumFluxForNormalization,double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS], ostream *poutspec, ostream *poutcontinuum);
        
    void calculateFluxCalibration(unsigned nPointsInReference,double *refwl,double *refflux,double referenceFluxForNormalization, unsigned binsize,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS],double uncalibratedContinuumFluxForNormalization,double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS], ostream *poutspec, ostream *poutcontinuum);
    
    void calculateFluxCalibration(unsigned nPointsInReference,double *refwl,double *refflux, unsigned binsize,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS], ostream *poutspec, ostream *poutcontinuum);    
    
    void applyFluxCalibration(double exposureTime, ostream *poutcontinuum);
    
    void applyFluxCalibration(double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS],double uncalibratedContinuumFluxForNormalization,double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS], bool absoluteCalibration, ostream *poutspec);
        
    void divideSpectralElementsBySEDElements(bool useThroughput, ostream *poutspec, bool StarPlusSky, bool starplusskyInvertSkyFiber);

    void multiplySpectralElementsBySEDElements(bool useThroughput,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS], ostream *poutspec);

    void multiplySpectralElementsBySEDElements(bool useThroughput,double spectralBinConstant, ostream *poutspec);
    
    void normalizeSpectralElementsByConstant(double maxFluxForNormalization, double maxBeamFluxForNormalization[MAXNUMBEROFBEAMS]);

    
	/*
	 * LE compatibility
	 */
	void applyFlatResponse(double exposureTime, operaSpectralElements *fluxCalibrationElements, ostream *poutspec);
    /*
     * Star+Sky Mode
     */
    void calculateStarAndSkyElements(bool starplusskyInvertSkyFiber,ostream *poutspec);
    
	/*
	 * Polar Mode, take the mean of the two beams
	 */
	void calculatePolarElements(ostream *poutspec);
	/*
	 * Barycentric Wavelength Correction
	 */
    void applyBarycentricWavelengthCorrection(double RVcorrectionInMetersPerSecond);
	void setExtendedBarycentricWavelengthCorrection(double RVcorrectionInMetersPerSecond);
    void applyBarycentricWavelengthCorrectionFromExtendedRvel(void);
};

#endif

