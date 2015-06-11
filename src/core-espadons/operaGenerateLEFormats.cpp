/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaGenerateLEFormats
 Version: 1.0
 Description:  Module to generate LE compatible formats.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

unsigned readLEorderwavelength(string LEorderwavelength, int *orders, double *wl0, double *wlf);

/*! \file operaGenerateLEFormats.cpp */

using namespace std;

/*! 
 * operaGenerateLEFormats
 * \author Eder Martioli
 * \brief  Module to generate LE compatible formats
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
    
	string inputOperaSpectrum; 
	string outputLEfilename;
    string LEorderwavelength;
    string object = "Nowhere";
    unsigned LibreEspritSpectrumType_val = LibreEspritsp2Spectrum;
    unsigned fluxType_val = RawFluxInElectronsPerElement;
    unsigned wavelengthType_val = ThArCalibratedInNM;
	int ordernumber = NOTPROVIDED;
    int minorder = 0;
    int maxorder = 0;
	
	args.AddRequiredArgument("inputOperaSpectrum", inputOperaSpectrum, "Extended opera spectrum (.spc)");
	args.AddRequiredArgument("outputLEfilename", outputLEfilename, "Libre-Esprit spectrum (.s)");
	args.AddOptionalArgument("LEorderwavelength", LEorderwavelength, "", "Table with LE order wavelength ranges");
	args.AddRequiredArgument("object", object, "Object name");
	args.AddRequiredArgument("LibreEspritSpectrumType", LibreEspritSpectrumType_val, "Spectrum type");
	args.AddRequiredArgument("fluxType", fluxType_val, "Flux type");
	args.AddRequiredArgument("wavelengthType", wavelengthType_val, "Wavelength type");
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	
	try {
		args.Parse(argc, argv);
		
		operaSpectralOrder_t LibreEspritSpectrumType = operaSpectralOrder_t(LibreEspritSpectrumType_val);
		/*  Available LibreEspritSpectrumType options for LE formats are:
			LibreEspritpolarimetry
			LibreEspritpolSpectrum
			LibreEspritsp1Spectrum
			LibreEspritsp2Spectrum
		*/
		operaFluxType_t fluxType = operaFluxType_t(fluxType_val);
		/*  Available operaFluxType_t options are:
			1 = RawFluxInElectronsPerElement
			2 = NormalizedFluxToContinuum
			3 = CalibratedFluxNormalizedToRefWavelength
		*/
		operaWavelengthType_t wavelengthType = operaWavelengthType_t(wavelengthType_val);
		/*  Available operaWavelengthType_t options are:
			1 = ThArCalibratedInNM
			2 = TelluricCorrectedWavelengthInNM
			3 = RVCorrectedWavelengthInNM
			4 = RVAndTelluricCorrectedWavelengthInNM
		*/
		
		// we need an input .e spectrum...
		if (inputOperaSpectrum.empty()) {
			throw operaException("operaGenerateLEFormats: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputLEfilename.empty()) {
			throw operaException("operaGenerateLEFormats: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		if (args.verbose) {
			cout << "operaGenerateLEFormats: inputOperaSpectrum = " << inputOperaSpectrum << endl; 
            cout << "operaGenerateLEFormats: outputLEfilename = " << outputLEfilename << endl;
            cout << "operaGenerateLEFormats: LEorderwavelength = " << LEorderwavelength << endl;
			cout << "operaGenerateLEFormats: object = " << object << endl;
			cout << "operaGenerateLEFormats: LibreEspritSpectrumType = " << LibreEspritSpectrumType << endl;
			cout << "operaGenerateLEFormats: fluxType = " << fluxType << endl;
			cout << "operaGenerateLEFormats: wavelengthType = " << wavelengthType << endl;
		}
        
		/*
		 * Down to business, read in all the source and calibration data.
		 */        
		operaSpectralOrderVector spectralOrders(inputOperaSpectrum);
		
		UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaGenerateLEFormats: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        unsigned LEnp = 0;
        int *LEorders = new int[MAXORDERS];
        double *LEwl0 = new double[MAXORDERS];
        double *LEwlf = new double[MAXORDERS];
        
        if (!LEorderwavelength.empty()) {
            LEnp = readLEorderwavelength(LEorderwavelength, LEorders, LEwl0, LEwlf);
        }

		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
			if (spectralOrder->gethasSpectralElements()) {
                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                unsigned nElements = spectralElements->getnSpectralElements();
                
                double wl0 = spectralElements->getwavelength(0);
                double wlf = spectralElements->getwavelength(nElements-1);
                
                for (unsigned i=0; i<LEnp; i++) {
                    if(order == LEorders[i]) {
                        if (LEwl0[i] > wl0 && LEwl0[i] < wlf) wl0 = LEwl0[i];
                        if (LEwlf[i] < wlf && LEwlf[i] > wl0) wlf = LEwlf[i];
                        break;
                    }
                }
                
                unsigned newIndexElem = 0;
                
                for(unsigned indexElem=0;indexElem<nElements;indexElem++) {
                    
                    double wl = spectralElements->getwavelength(indexElem);
                    double flux = spectralElements->getFlux(indexElem);
                    double fluxVariance = spectralElements->getFluxVariance(indexElem);
                    double tell = spectralElements->gettell(indexElem);
                    double rvel = spectralElements->getrvel(indexElem);
                    double normalizedFlux = spectralElements->getnormalizedFlux(indexElem);
                    double normalizedFluxVariance = spectralElements->getnormalizedFluxVariance(indexElem);
                    double fcalFlux = spectralElements->getfcalFlux(indexElem);
                    double fcalFluxVariance = spectralElements->getfcalFluxVariance(indexElem);
                    double rawFlux = spectralElements->getrawFlux(indexElem);
                    double rawFluxVariance = spectralElements->getrawFluxVariance(indexElem);
                    
                    if (wl > wl0 && wl < wlf) {
                        spectralElements->setwavelength(wl,newIndexElem);
                        spectralElements->setFlux(flux,newIndexElem);
                        spectralElements->setFluxVariance(fluxVariance,newIndexElem);
                        spectralElements->settell(tell,newIndexElem);
                        spectralElements->setrvel(rvel,newIndexElem);
                        spectralElements->setnormalizedFlux(normalizedFlux,newIndexElem);
                        spectralElements->setfcalFlux(fcalFlux,newIndexElem);
                        spectralElements->setrawFlux(rawFlux,newIndexElem);
                        spectralElements->setnormalizedFluxVariance(normalizedFluxVariance,newIndexElem);
                        spectralElements->setfcalFluxVariance(fcalFluxVariance,newIndexElem);
                        spectralElements->setrawFluxVariance(rawFluxVariance,newIndexElem);
                        
                        newIndexElem++;
                    } else if (wl > wl0 && wl > wlf) {
                        break;
                    }
                }
                
                spectralElements->setnSpectralElements(newIndexElem);
                
                if (spectralElements->getHasExtendedBeamFlux()){
                    switch (fluxType) {
                        case RawFluxInElectronsPerElement:
                            //spectralElements->copyFROMrawFlux();
                            break;
                        case NormalizedFluxToContinuum:
                            spectralElements->copyFROMnormalizedFlux();
                            break;
                        case CalibratedFluxNormalizedToRefWavelength:
                            spectralElements->copyFROMfcalFlux();
                            break;
                        default:
                            break;
                    }
                    switch (wavelengthType) {
                        case ThArCalibratedInNM:
                            break;
                        case TelluricCorrectedWavelengthInNM:
                            spectralElements->copyFROMtell();
                            break;
                        case RVCorrectedWavelengthInNM:
                            spectralOrder->applyWavelengthCorrectionFromExtendedRvel();
                            break;
                        case RVAndTelluricCorrectedWavelengthInNM:
                            spectralElements->copyFROMtell();
                            spectralOrder->applyWavelengthCorrectionFromExtendedRvel();
                            break;
                        default:
                            break;
                    }
                }
			}
		}        
 		// output wavelength/flux calibrated spectrum...
		spectralOrders.setObject(object);
		spectralOrders.WriteSpectralOrders(outputLEfilename, LibreEspritSpectrumType);
	}
	catch (operaException e) {
		cerr << "operaGenerateLEFormats: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaGenerateLEFormats: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*
 * Read LE order wavelength ranges
 */
unsigned readLEorderwavelength(string LEorderwavelength, int *orders, double *wl0, double *wlf) {
    ifstream astream;
    string dataline;
    
    int tmporder = 0;
    double tmpwl0 = -1.0;
    double tmpwlf = -1.0;
    unsigned np = 0;
    
    astream.open(LEorderwavelength.c_str());
    if (astream.is_open()) {
        while (astream.good()) {
            getline(astream, dataline);
            if (strlen(dataline.c_str())) {
                if (dataline.c_str()[0] == '#') {
                    // skip comments
                } else {
                    sscanf(dataline.c_str(), "%d %lf %lf", &tmporder, &tmpwl0, &tmpwlf);
                    orders[np] = tmporder,
                    wl0[np] = tmpwl0;
                    wlf[np] = tmpwlf;
                    np++;
                }	// skip comments
            }
        } // while (astream.good())
        astream.close();
    }	// if (astream.open())
    return np;
}
