/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaCreateProduct
 Version: 1.0
 Description: Bundle .s and .sn files into a product
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


#include <stdio.h>
#include <getopt.h>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/gzstream.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaStokesVector.h"
#include "libraries/operaPolarimetry.h"
#include "libraries/Polynomial.h"
#include "libraries/LaurentPolynomial.h"
#include "libraries/operaCCD.h"					// for MAXORDERS
#include "libraries/operastringstream.h"		// for Double, Float
#include "libraries/operaArgumentHandler.h"

/*! \file operaCreateProduct.cpp */

using namespace std;

/*!
 * operaCreateProduct
 * \author Doug Teeple
 * \brief Bundle files into an i.fits, p.fits, m.fits Product.
 * \arg argc
 * \arg argv
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

operaArgumentHandler args;

// Returns the number of rows from a Libre-Esprit spectrum file (.s).
unsigned int GetRowCountFromLESpectrumFile(const string spectrumfile);

// Returns the number of rows from an extended spectrum file (.es).
unsigned int GetRowCountFromExtendedSpectrumFile(const string spectrumfile, const instrumentmode_t instrumentmode);

/* Updates the specified columns of a FITS product using the values from a Libre-Esprit spectrum file (.s).
   colgroup - which group of columns (index from 0 to 3) will be updated 
   groupsize - the number of columns in this group */
void UpdateProductFromLESpectrumFile(operaFITSProduct& Product, const string spectrumfile, const unsigned groupsize, const unsigned int colgroup, const unsigned int rowcount);

/* Updates the specified columns of a FITS product using the values from an extended spectrum file (.es).
   colgroup - which group of columns (index from 0 to 3) will be updated */
void UpdateProductFromExtendedSpectrumFile(operaFITSProduct& Product, const string spectrumfile, const instrumentmode_t instrumentmode, const unsigned int colgroup);

double readRadialVelocityCorrection(string filename) {
	if (!filename.empty()) {
		operaistream fin(filename.c_str());
		if (fin.is_open()) {
			while (fin.good()) {
				string dataline;
				getline (fin,dataline);
				if (!dataline.empty() && dataline[0] != '#') {
					stringstream ss (dataline);
					double rvel;
					ss >> rvel;
					return rvel;
				}									
			}
		}
	}
	return 0.00;
}

// Adds various information to the header of the FITS product.
void AddFITSHeaderToProduct(operaFITSProduct& Product, const string version, const string date, const string spectralOrderType, const string parametersfilename, const string snrfilename, const string rvelfilename) {
	Product.operaFITSDeleteHeaderKey("DATASEC");
	Product.operaFITSDeleteHeaderKey("DETSEC");
	Product.operaFITSAddComment("----------------------------------------------------");
	Product.operaFITSAddComment("| Processed by the CFHT OPERA Open Source Pipeline |");
	Product.operaFITSAddComment("----------------------------------------------------");
	Product.operaFITSAddComment(version);
	Product.operaFITSAddComment("Processing Date");
	Product.operaFITSAddComment("---------------");
	Product.operaFITSAddComment(date);
	Product.operaFITSAddComment("------------------------------------------------------------------------");
	Product.operaFITSAddComment("upena-compatible headers for column names in primary extension:");
	Product.operaFITSAddComment("(1) Spectroscopy Star only mode");
	Product.operaFITSAddComment("    First column = wavelength in nanometres");
	Product.operaFITSAddComment("    Second column = intensity");
	Product.operaFITSAddComment("    Third column = error bar");
	Product.operaFITSAddComment("(2) Polarimetry");
	Product.operaFITSAddComment("    1st col = wavelength in nanometres");
	Product.operaFITSAddComment("    2d  col = intensity");
	Product.operaFITSAddComment("    3rd col = polarisation (Q or U or V or W)");
	Product.operaFITSAddComment("    4th col = Check Spectra #1");
	Product.operaFITSAddComment("    5th col = Check Spectra #2");
	Product.operaFITSAddComment("     6th col = error bar");
	Product.operaFITSAddComment("(3) Spectroscopy Star + Sky");
	Product.operaFITSAddComment("    1st col = wavelength");
	Product.operaFITSAddComment("    2d  col = star spectra (sky subtracted)");
	Product.operaFITSAddComment("    3rd col = star + sky spectra");
	Product.operaFITSAddComment("    4th col = sky spectra");
	Product.operaFITSAddComment("    5, 6, 7 = error bars for each column 2, 3, 4");
	Product.operaFITSAddComment("------------------------------------------------------------------------");
	Product.operaFITSAddComment(spectralOrderType);
	Product.operaFITSAddComment("OPERA Processing Parameters");
	Product.operaFITSAddComment("---------------------------");
	if (!parametersfilename.empty()) {
		if (args.verbose) cout << "operaCreateProduct: adding parameters " << endl;
		ifstream parameters(parametersfilename.c_str());
		while (parameters.good()) {
			string dataline;
			getline(parameters, dataline);
			if (strlen(dataline.c_str())) Product.operaFITSAddComment(dataline);
		}
	}
	//To do: replace this with more useful SNR information (i.e. peak SNR)
	if (!snrfilename.empty()) {
		if (args.verbose) cout << "operaCreateProduct: adding SNR comments " << endl;
		operaSpectralOrderVector spectralOrders(snrfilename);
		unsigned minorder = spectralOrders.getMinorder();
		unsigned maxorder = spectralOrders.getMaxorder();
		Product.operaFITSAddComment("SNR Table");
		Product.operaFITSAddComment("---------");
		Product.operaFITSAddComment("Format: <order number><center SNR><center wavelength><SNR>");
		for (unsigned order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if (spectralOrder->gethasSpectralElements()) {
				operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
				string out = itos(spectralOrder->getorder()) + ' ' + ftos(spectralOrder->getCenterSNR()) + ' ' + ftos(spectralelements->getwavelength(spectralelements->getnSpectralElements()/2)) + ' ' + ftos(spectralelements->getFluxSNR(spectralelements->getnSpectralElements()/2));
				Product.operaFITSAddComment(out.c_str());
			}
		}
	}
	double rvcorr = readRadialVelocityCorrection(rvelfilename);
	Product.operaFITSSetHeaderValue("HRV", rvcorr, "barycentric RV correction");
	//Product.operaFITSSetHeaderValue("TELLRV", value, "telluric RV correction");
}

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
	string version, date;
	args.AddOptionalArgument("version", version, "", "");
	args.AddOptionalArgument("date", date, "", "");
	
	string inputfilename, outputfilename;
	args.AddRequiredArgument("input", inputfilename, "input file (o.fits)");
	args.AddOptionalArgument("output", outputfilename, "", "output file (i.fits/m.fits)");
	string ifilename;
	args.AddOptionalArgument("i", ifilename, "", "i.s");
	string polarfilename;
	args.AddOptionalArgument("polar", polarfilename, "", ".ep");
	string ufile, nfile, uwfile, nwfile;
	args.AddOptionalArgument("ufile", ufile, "", "unnormalized spectrum (iu.s/pu.s)");
	args.AddOptionalArgument("nfile", nfile, "", "normalized spectrum (in.s/pn.s)");
	args.AddOptionalArgument("uwfile", uwfile, "", "unnormalized wavelength corrected spectrum (iuw.s/puw.s)");
	args.AddOptionalArgument("nwfile", nwfile, "", "normalized wavelength corrected spectrum (inw.s/pnw.s)");
	string snrfilename;
	args.AddOptionalArgument("snr", snrfilename, "", ".sn");
	string csvfilename;
	args.AddOptionalArgument("csv", csvfilename, "", ".csv");
	string esfilename;
	args.AddOptionalArgument("es", esfilename, "", ".es.gz");
	string geomfilename;
	args.AddOptionalArgument("geom", geomfilename, "", ".geom");
	string wavefilename;
	args.AddOptionalArgument("wave", wavefilename, "", ".wcal");
	string aperfilename;
	args.AddOptionalArgument("aper", aperfilename, "", ".aper");
	string proffilename;
	args.AddOptionalArgument("prof", proffilename, "", ".prof");
	string ordpfilename;
	args.AddOptionalArgument("ordp", ordpfilename, "", ".ordp");
	string beamfilename;
	args.AddOptionalArgument("beam", beamfilename, "", ".es");
	string gainfilename;
	args.AddOptionalArgument("gain", gainfilename, "", ".gain");
	string biasfilename;
	args.AddOptionalArgument("bias", biasfilename, "", ".bias");
	string fcalfilename;
	args.AddOptionalArgument("fcal", fcalfilename, "", ".fcal");
	string dispfilename;
	args.AddOptionalArgument("disp", dispfilename, "", ".disp");
	string rvelfilename;
	args.AddOptionalArgument("rvel", rvelfilename, "", "i.rvel");
	string tellfilename;
	args.AddOptionalArgument("tell", tellfilename, "", "i.tell");
	string prvelfilename;
	args.AddOptionalArgument("prvel", prvelfilename, "", "p.rvel");
	string ptellfilename;
	args.AddOptionalArgument("ptell", ptellfilename, "", "p.tell");
	string parametersfilename;
	args.AddOptionalArgument("parameters", parametersfilename, "", ".parm");
	string object;
	args.AddOptionalArgument("object", object, "", "object name, needed for Libre-Esprit output");
	
	string spectralOrderType;
	args.AddOptionalArgument("spectrumtype", spectralOrderType, "", "spectral order type");
	int compressionVal;
	args.AddOptionalArgument("compressiontype", compressionVal, cNone, "compression type");
	bool centralsnr;
	args.AddSwitch("centralsnr", centralsnr, "");
	unsigned sequence;
	args.AddOptionalArgument("sequence", sequence, 0, "for polar sequence in case of masterfluxcalibrations");
		
	try {
		args.Parse(argc, argv);
		
		if (inputfilename.empty()) // we need an input...
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		if (outputfilename.empty() && csvfilename.empty()) // we need an output...
			throw operaException("operaCreateProduct: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		
		eCompression compression = (eCompression)compressionVal;
		unsigned extensions = 0;	
		if (!beamfilename.empty()) extensions++;
		if (!polarfilename.empty()) extensions++;
		if (!geomfilename.empty()) extensions++;
		if (!wavefilename.empty()) extensions++;
		if (!aperfilename.empty()) extensions++;
		if (!proffilename.empty()) extensions++;
		if (!ordpfilename.empty()) extensions++;
		if (!gainfilename.empty()) extensions++;
		if (!biasfilename.empty()) extensions++;
		if (!fcalfilename.empty()) extensions++;
		if (!dispfilename.empty()) extensions++;
		if (!rvelfilename.empty()) extensions++;
		if (!tellfilename.empty()) extensions++;
		if (!prvelfilename.empty()) extensions++;
		if (!ptellfilename.empty()) extensions++;
		if (args.verbose) {
			cout << "operaCreateProduct: input= " << inputfilename << endl;
			cout << "operaCreateProduct: output= " << outputfilename << endl;
			cout << "operaCreateProduct: ifilename= " << ifilename << endl;
			cout << "operaCreateProduct: ufile= " << ufile << endl;
			cout << "operaCreateProduct: nfile= " << nfile << endl;
			cout << "operaCreateProduct: uwfile= " << uwfile << endl;
			cout << "operaCreateProduct: nwfile= " << nwfile << endl;
			cout << "operaCreateProduct: snrfilename= " << snrfilename << endl;
			cout << "operaCreateProduct: csvfilename= " << csvfilename << endl;
			cout << "operaCreateProduct: esfilename= " << esfilename << endl;
			cout << "operaCreateProduct: parametersfilename= " << parametersfilename << endl;
			cout << "operaCreateProduct: geomfilename= " << geomfilename << endl;
			cout << "operaCreateProduct: wavefilename= " << wavefilename << endl;
			cout << "operaCreateProduct: aperfilename= " << aperfilename << endl;
			cout << "operaCreateProduct: proffilename= " << proffilename << endl;
			cout << "operaCreateProduct: ordpfilename= " << ordpfilename << endl;
			cout << "operaCreateProduct: beamfilename= " << beamfilename << endl;
			cout << "operaCreateProduct: polarfilename= " << polarfilename << endl;
			cout << "operaCreateProduct: gainfilename= " << gainfilename << endl;
			cout << "operaCreateProduct: biasfilename= " << biasfilename << endl;
			cout << "operaCreateProduct: fcalfilename= " << fcalfilename << endl;
			cout << "operaCreateProduct: dispfilename= " << dispfilename << endl;
			cout << "operaCreateProduct: rvelfilename= " << rvelfilename << endl;
			cout << "operaCreateProduct: tellfilename= " << tellfilename << endl;
			cout << "operaCreateProduct: prvelfilename= " << prvelfilename << endl;
			cout << "operaCreateProduct: ptellfilename= " << ptellfilename << endl;
			cout << "operaCreateProduct: OPERA version= " << version << endl;
			cout << "operaCreateProduct: Reduction date= " << date << endl;
			cout << "operaCreateProduct: compression= " << compression << endl;
			cout << "operaCreateProduct: spectralOrderType= " << spectralOrderType << endl;
			cout << "operaCreateProduct: extensions= " << (extensions+1) << endl;
		}
		
		operaFITSImage input(inputfilename, tfloat, READONLY, cNone, true);
		string mode = input.operaFITSGetHeaderValue("INSTMODE");
		input.operaFITSImageClose();
		if (args.verbose) {
			cout << "operaCreateProduct: " << mode << endl;
		}
		instrumentmode_t instrumentmode;
		unsigned int numcols = 0;
		if (mode.find("Polarimetry") != string::npos) {
			instrumentmode = MODE_POLAR;
			numcols = MODE_POLAR_COLS;
		} else if (mode.find("Spectroscopy, star+sky") != string::npos) {
			instrumentmode = MODE_STAR_PLUS_SKY;
			numcols = MODE_STAR_PLUS_SKY_COLS;
		} else if (mode.find("Spectroscopy, star only") != string::npos) {
			instrumentmode = MODE_STAR_ONLY;
			numcols = MODE_STAR_ONLY_COLS;
		} else {
			throw operaException("operaCreateProduct: "+mode+" ", operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
		}
		unsigned minorder = 0;
		unsigned maxorder = 0;
		/***********************************************************************************************
		 * PART 0 - Do the Intensity libre-esprit compatible products
		 ***********************************************************************************************/
		if (!ifilename.empty()) {
			operaostream fout;
			fout.open(outputfilename.c_str());
			if (args.verbose) {
				cout << "operaCreateProduct: mode=" << instrumentmode << endl;
				cout << "operaCreateProduct: object='" << object << "'"<< endl;
			}
			operaSpectralOrderVector spectralOrders(ifilename);
			minorder = spectralOrders.getMinorder();
			maxorder = spectralOrders.getMaxorder();
			unsigned rows = 0;
			for (unsigned order=MAXORDERS-1; order>=minorder; order--) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
					if (spectralOrder->gethasWavelength()) {
						rows += spectralelements->getnSpectralElements();
					}
				}
			}
			fout << "***Reduced spectrum of '" << object << "'" << endl;
			switch (instrumentmode) {
				case MODE_POLAR:
				case MODE_STAR_ONLY:
					fout << rows << " 2" << endl;
					for (unsigned order=MAXORDERS-1; order>=minorder; order--) {
						operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
						if (spectralOrder->gethasSpectralElements()) {
							operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
							if (spectralOrder->gethasWavelength()) {
								for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
									fout << fixed << setprecision(4)
									<< spectralelements->getwavelength(i) << ' ';
									fout << scientific << spectralelements->getFlux(i) << ' '
									<< sqrt(spectralelements->getFluxVariance(i)) << endl;
								}
								if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
							}
						}
					}
					break;
				case MODE_STAR_PLUS_SKY:
					fout << rows << " 6" << endl;
					for (unsigned order=MAXORDERS-1; order>=minorder; order--) {
						operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
						operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
						operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
						if (spectralOrder->gethasSpectralElements()) {
							operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
							if (spectralOrder->gethasWavelength()) {
								for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
									fout << fixed << setprecision(4)
									<< spectralelements->getwavelength(i) << ' ';
									fout << scientific << beamElements0->getFlux(i) << ' '
									<< (beamElements0->getFlux(i) + beamElements1->getFlux(i)) << ' '
									<< beamElements1->getFlux(i) << ' '
									<< sqrt(beamElements0->getFluxVariance(i)) << ' '
									<< sqrt(beamElements0->getFluxVariance(i)+beamElements1->getFluxVariance(i)) << ' '
									<< sqrt(beamElements1->getFluxVariance(i)) << endl;
								}
								if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
							}
						}
					}
					break;
				default:
					break;
			}
			fout.close();
			exit(0);
		}
		else if (!polarfilename.empty() && outputfilename.find("m.fits") == string::npos) {
			operaostream fout;
			fout.open(outputfilename.c_str());
			if (args.verbose) {
				cout << "operaCreateProduct: mode=" << instrumentmode << endl;
				cout << "operaCreateProduct: object='" << object << "'"<< endl;
			}
			operaSpectralOrderVector spectralOrders(polarfilename);
			minorder = spectralOrders.getMinorder();
			maxorder = spectralOrders.getMaxorder();
			unsigned rows = 0;
			for (unsigned order=minorder; order<=maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasPolarimetry()) {
					if (spectralOrder->gethasWavelength()) {
						rows += spectralOrder->getPolarimetry()->getLength();
					}
				}
			}
			fout << "***Reduced spectrum of '" << object << "'" << endl;
			fout << rows << " 5" << endl;
			double PolarizationVariance, DegreeOfPolarization, DegreeOfPolarizationVariance, NullSpectrum1Variance, NullSpectrum2Variance;
			double Intensity, Variance, Polarization, NullSpectrum1, NullSpectrum2;
			operaPolarimetry *Polarimetry = NULL;
			for (unsigned order=MAXORDERS-1; order>=minorder; order--) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
					if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasWavelength()) {
						Polarimetry = spectralOrder->getPolarimetry();
						unsigned length = spectralOrder->getPolarimetry()->getLength();
						stokes_parameter_t stokesParameter = StokesI;
						if (Polarimetry->getHasStokesV()) {
							stokesParameter = StokesV;
						} else if (Polarimetry->getHasStokesQ()) {
							stokesParameter = StokesQ;
						} else if (Polarimetry->getHasStokesU()) {
							stokesParameter = StokesU;
						} else if (Polarimetry->getHasStokesI()) {
							stokesParameter = StokesI;
						}
						for (unsigned index = 0 ; index < length ; index++) {
							Intensity = Polarimetry->getStokesParameter(StokesI)->getflux(index);
							Variance = Polarimetry->getStokesParameter(StokesI)->getvariance(index);
							Polarization = Polarimetry->getStokesParameter(stokesParameter)->getflux(index);
							PolarizationVariance = Polarimetry->getStokesParameter(stokesParameter)->getvariance(index);
							DegreeOfPolarization = Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(index);
							DegreeOfPolarizationVariance = Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(index);
							NullSpectrum1 = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(index);
							NullSpectrum1Variance = Polarimetry->getFirstNullPolarization(stokesParameter)->getvariance(index);
							NullSpectrum2 = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(index);
							NullSpectrum2Variance = Polarimetry->getSecondNullPolarization(stokesParameter)->getvariance(index);
							fout << fixed << setprecision(4)
							<< SpectralElements->getwavelength(index) << ' ';
							fout << scientific << Intensity << ' '
							<< DegreeOfPolarization << ' '
							<< NullSpectrum1 << ' '
							<< NullSpectrum2 << ' '
							<< sqrt(DegreeOfPolarizationVariance) << endl;
						}
						if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
					}
				}
			}
			fout.close();
			exit(0);
		}

		/***********************************************************************************************
		 * PARTS I & II - Create the Intensity product i.fits or the Polarimetry product p.fits
		 ***********************************************************************************************/
		else if (!ufile.empty() && !nfile.empty() && ! uwfile.empty() && !nwfile.empty()) {
			const bool extended = (ufile.find(".es") != string::npos); // Check if using extended spectrum
			const unsigned int datapoints = (extended ? GetRowCountFromExtendedSpectrumFile(ufile, instrumentmode) : GetRowCountFromLESpectrumFile(ufile));
			operaFITSProduct Product(outputfilename, inputfilename, instrumentmode, datapoints, numcols, compression);
			if(extended) { //Using .es extended spectrum
				UpdateProductFromExtendedSpectrumFile(Product, ufile, instrumentmode, 3);
				UpdateProductFromExtendedSpectrumFile(Product, nfile, instrumentmode, 2);
				UpdateProductFromExtendedSpectrumFile(Product, uwfile, instrumentmode, 1);
				UpdateProductFromExtendedSpectrumFile(Product, nwfile, instrumentmode, 0);
			}
			else { //Using .s Libre-Esprit spectrum (has no order information)
				UpdateProductFromLESpectrumFile(Product, ufile, numcols/4, 3, datapoints);
				UpdateProductFromLESpectrumFile(Product, nfile, numcols/4, 2, datapoints);
				UpdateProductFromLESpectrumFile(Product, uwfile, numcols/4, 1, datapoints);
				UpdateProductFromLESpectrumFile(Product, nwfile, numcols/4, 0, datapoints);
			}
			AddFITSHeaderToProduct(Product, version, date, spectralOrderType, parametersfilename, snrfilename, rvelfilename);
			Product.operaFITSImageSave();
			Product.operaFITSImageClose();
			if (args.verbose && instrumentmode == MODE_POLAR) cout << "operaCreateProduct: done polarimetry " << endl;
			else if (args.verbose) cout << "operaCreateProduct: done intensity " << endl;
        }
		
        /***********************************************************************************************
         * PART III - Do the CSV
         ***********************************************************************************************/
		else if (!csvfilename.empty() && !esfilename.empty()) {
			if (args.verbose) {
				cout << "operaCreateProduct: mode=" << instrumentmode << endl;
				cout << "operaCreateProduct: object='" << object << "'"<< endl;
			}
			operaSpectralOrderVector spectralOrders(esfilename);
			if (!polarfilename.empty()) {
				spectralOrders.ReadSpectralOrders(polarfilename);
			}
			spectralOrders.setInstrumentmode(instrumentmode);
			spectralOrders.setObject(object);
			spectralOrders.WriteSpectralOrders(csvfilename, CSV);
		}
		/***********************************************************************************************
         * PART IV - Do the calibration MEF
         ***********************************************************************************************/
        else if (extensions) {
			unsigned extension = 0;
			compression = cNone;	// cfitsio creates bintables for compressed MEFS, so we can't compress
			operaMEFFITSProduct *product = NULL;
			unsigned Columns = 44;
			unsigned MaxColumns = 44;
			DATASEC_t detsec = {1,Columns,1,1};
			if (!beamfilename.empty()) {
				// <order number><nElements><nBeams><elementindex><SpectralElements photoCenterX><SpectralElements photoCenterY><SpectralElements dist><SpectralElements flux><SpectralElements flux variance><XCorrelation><nBeams><beam><BeamElements[beam] photoCenterX><BeamElements[beam] photoCenterY><BeamElements[beam] flux><BeamElements[beam] flux variance>
				extension++;
				operaSpectralOrderVector spectralOrders(beamfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasSpectralElements()) {
						Rows += spectralOrder->getSpectralElements()->getnSpectralElements();
					}
				}
				product = new operaMEFFITSProduct(outputfilename, inputfilename, Columns, Rows, extension, compression);
				operaFITSProduct Product(*product);
				/*
				 * Now populate headers as comments
				 */
				product->operaFITSDeleteHeaderKey("DATASEC");
				product->operaFITSDeleteHeaderKey("DETSEC");
				product->operaFITSAddComment("----------------------------------------------------");
				product->operaFITSAddComment("| Processed by the CFHT OPERA Open Source Pipeline |");
				product->operaFITSAddComment("----------------------------------------------------");
				product->operaFITSAddComment(version);
				product->operaFITSAddComment("Processing Date");
				product->operaFITSAddComment("---------------");
				product->operaFITSAddComment(date);
				product->operaFITSAddComment("OPERA Processing Parameters");
				product->operaFITSAddComment("---------------------------");
				/*
				 * Add in the reduction parameters
				 */
				if (!parametersfilename.empty()) {
					if (args.verbose) {
						cout << "operaCreateProduct: adding parameters " << endl;
					}
					ifstream parameters(parametersfilename.c_str());
					string dataline;
					while (parameters.good()) {
						getline(parameters, dataline);
						if (strlen(dataline.c_str())) {
							product->operaFITSAddComment(dataline);
						}
					}
				}
				MaxColumns = 25;
				detsec.x2 = MaxColumns;
				detsec.y2 = Rows;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " beam spectra " << spectralOrderType << " " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("BEAMFLUX", MaxColumns, Rows, detsec, true);	// reuse the first extension
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("OPERA Calibration Reduction Data", extension);
				product->operaFITSAddComment("----------------------", extension);
				product->operaFITSAddComment("Beam Flux Variance", extension);
				product->operaFITSAddComment("--------------------", extension);
				product->operaFITSAddComment("ORDERTYPE: " + lowerCase(spectralOrderType), extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><nElements><nBeams><elementindex><SpectralElements photoCenterX><SpectralElements photoCenterY><SpectralElements dist><SpectralElements wl><SpectralElements flux><SpectralElements flux variance><XCorrelation><nBeams>[repeat <beam><BeamElements[beam] photoCenterX><BeamElements[beam] photoCenterY><BeamElements[beam] flux><BeamElements[beam] flux variance>]", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0; indexElem < spectralelements->getnSpectralElements(); indexElem++) {
							Product[Row][Column++] = maxorder - minorder + 1;
							Product[Row][Column++] = order;
							Product[Row][Column++] = spectralelements->getnSpectralElements();
							Product[Row][Column++] = spectralOrder->getnumberOfBeams();
							Product[Row][Column++] = indexElem;
							Product[Row][Column++] = spectralelements->getphotoCenterX(indexElem);
							Product[Row][Column++] = spectralelements->getphotoCenterY(indexElem);
							Product[Row][Column++] = spectralelements->getdistd(indexElem);
							Product[Row][Column++] = spectralelements->getwavelength(indexElem);
							Product[Row][Column++] = spectralelements->getFlux(indexElem);
							Product[Row][Column++] = spectralelements->getFluxVariance(indexElem);
							Product[Row][Column++] = spectralelements->getXCorrelation(indexElem);
							Product[Row][Column++] = spectralOrder->getnumberOfBeams();
							for (unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								Product[Row][Column++] = beam;
								Product[Row][Column++] = spectralOrder->getBeamElements(beam)->getphotoCenterX(indexElem);
								Product[Row][Column++] = spectralOrder->getBeamElements(beam)->getphotoCenterY(indexElem);
								Product[Row][Column++] = spectralOrder->getBeamElements(beam)->getFlux(indexElem);
								Product[Row][Column++] = spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem);
								if (Column >= MaxColumns) {
									cout << "operaCreateProduct: ***Warning: BEAMFLUX increase the column size to at least " << Column << endl;
								}
							}
							Row++;
							Column = 0;
						}
					}
				}
				product->saveExtension(extension);
			}
			if (!snrfilename.empty()) {
				extension++;
				// <order number><Center SNR><wl><SNR>
				operaSpectralOrderVector spectralOrders(snrfilename);
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (centralsnr) {
						Rows++;
					} else if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						Rows += spectralelements->getnSpectralElements();
					}
				}
                if (centralsnr) {
                    MaxColumns = 5;
				} else {
                    MaxColumns = 7;
                }
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " SNR " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("SNR", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
                if (centralsnr) {
                    product->operaFITSAddComment("Center Signal to Noise Ratio", extension);
				} else {
                    product->operaFITSAddComment("Signal to Noise Ratio", extension);
                }
				product->operaFITSAddComment("---------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
                if (centralsnr) {
                    product->operaFITSAddComment("<orders><Columns><order number><wavelength><Center SNR>", extension);
				} else {
                    product->operaFITSAddComment("<orders><Columns><order number><nElements><Center SNR><wavelength><SNR>", extension);
                }
                product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
					if (centralsnr && spectralOrder->gethasSpectralElements()) {
						Product[Row][Column++] = maxorder - minorder + 1;
						Product[Row][Column++] = MaxColumns;
						Product[Row][Column++] = order;
						Product[Row][Column++] = spectralElements->getwavelength(spectralOrder->getSpectralElements()->getnSpectralElements()/2);
						Product[Row][Column++] = spectralOrder->getCenterSNR();
						Column = 0;
						Row++;
					} else if (spectralOrder->gethasSpectralElements()) {
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							Product[Row][Column++] = maxorder - minorder + 1;
							Product[Row][Column++] = MaxColumns;
							Product[Row][Column++] = order;
							Product[Row][Column++] = spectralElements->getnSpectralElements();
							Product[Row][Column++] = spectralOrder->getCenterSNR();
							Product[Row][Column++] = spectralElements->getwavelength(k);
							Product[Row][Column++] = spectralElements->getFluxSNR(k);
							if (Column >= MaxColumns) {
								cout << "operaCreateProduct: ***Warning: SNR increase the column size to at least " << Column << endl;
							}
							Column = 0;
							Row++;
						}
					}
				}
				product->saveExtension(extension);
			}
			double rvcorr = readRadialVelocityCorrection(rvelfilename);
			product->operaFITSSetHeaderValue("HRV", rvcorr, "barycentric RV correction");
			if (!polarfilename.empty()) {
				extension++;
				// <order number> <StokesParameter_t> <length> <index> <wavelength> <Stokes(Q,U,V) flux> <Stokes(Q,U,V) variance> <StokesI flux> <StokesI variance> <degree of polarization flux> <degree of polarization variance> <first null polarization> <first null polarization variance> <second null polarization> <second null polarization variance>
				operaSpectralOrderVector spectralOrders(polarfilename);
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasPolarimetry()) {
						operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
						Rows += Polarimetry->getLength();
					}
				}
				MaxColumns = 20;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " POLAR " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("POLAR", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Polarimetry", extension);
				product->operaFITSAddComment("-----------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><StokesParameter_t><method_t><length><distance><wavelength><Stokes(Q,U,V) flux><Stokes(Q,U,V) variance><StokesI flux><StokesI variance><degree of polarization flux><degree of polarization variance> <first null polarization><first null polarization variance><second null polarization><second null polarization variance>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasSpectralElements()) {
						operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
						operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
						operaStokesVector *StokesVector = Polarimetry->getStokesVector();
						operaStokesVector *DegreeOfPolarization = Polarimetry->getDegreeOfPolarization();
						operaStokesVector *FirstNullPolarization = Polarimetry->getFirstNullPolarization();
						operaStokesVector *SecondNullPolarization = Polarimetry->getSecondNullPolarization();
						SpectralElements->setHasDistance(true);
						SpectralElements->setHasWavelength(true);
						SpectralElements->setHasXCorrelation(true);
						unsigned length = StokesVector->getLength();
						for (unsigned index = 0 ; index < length; index++) {
							Product[Row][Column++] = maxorder - minorder + 1;
							Product[Row][Column++] = order;
							stokes_parameter_t stokes = StokesI;
							if (Polarimetry->getHasStokesQ()) {
								stokes = StokesQ;
							}
							if (Polarimetry->getHasStokesU()) {
								stokes = StokesU;
							}
							if (Polarimetry->getHasStokesV()) {
								stokes = StokesV;
							}
							Product[Row][Column++] = stokes;
							Product[Row][Column++] = Polarimetry->getmethod();
							Product[Row][Column++] = length;
							Product[Row][Column++] = SpectralElements->getdistd(index);
							Product[Row][Column++] = SpectralElements->getwavelength(index);
							Product[Row][Column++] = SpectralElements->getXCorrelation(index);
							Product[Row][Column++] = StokesVector->getStokesParameterFlux(stokes, index);
							Product[Row][Column++] = StokesVector->getStokesParameterVariance(stokes, index);
							Product[Row][Column++] = StokesVector->getStokesParameterFlux(StokesI, index);
							Product[Row][Column++] = StokesVector->getStokesParameterVariance(StokesI, index);
							Product[Row][Column++] = DegreeOfPolarization->getStokesParameterFlux(stokes, index);
							Product[Row][Column++] = DegreeOfPolarization->getStokesParameterVariance(stokes, index);
							if (Polarimetry->getHasFirstNullPolarization()) {
								Product[Row][Column++] = FirstNullPolarization->getStokesParameterFlux(stokes, index);
								Product[Row][Column++] = FirstNullPolarization->getStokesParameterVariance(stokes, index);
							} else {
								Product[Row][Column++] = 0.0;
								Product[Row][Column++] = 0.0;
							}
							if (Polarimetry->getHasSecondNullPolarization()) {
								Product[Row][Column++] = SecondNullPolarization->getStokesParameterFlux(stokes, index);
								Product[Row][Column++] = SecondNullPolarization->getStokesParameterVariance(stokes, index);
								
							} else {
								Product[Row][Column++] = 0.0;
								Product[Row][Column++] = 0.0;
							}
							if (Column >= MaxColumns) {
								cout << "operaCreateProduct: ***Warning: POLAR increase the column size to at least " << Column << endl;
							}
							Column = 0;
							Row++;
						}
					}
				}
				product->saveExtension(extension);
			}
			if (!geomfilename.empty()) {
				extension++;
				// <order number><number of coefficients><ndatapoints> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL <chisqr><YBinning><miny><maxy>
				operaSpectralOrderVector spectralOrders(geomfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasGeometry()) {
						Rows++;
					}
				}
				MaxColumns = 9+MAXPOLYNOMIAL*2;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " geometry " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("GEOMETRY", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Geometry Polynomial", extension);
				product->operaFITSAddComment("-------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><number of coefficients><ndatapoints>[<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL <chisqr><YBinning><miny><maxy>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasGeometry()) {
						operaGeometry *geometry = spectralOrder->getGeometry();
						Polynomial *polynomial = geometry->getCenterPolynomial();
						unsigned npar = polynomial->getOrderOfPolynomial();
						Product[Row][Column++] = maxorder - minorder + 1;
						Product[Row][Column++] = order;
						Product[Row][Column++] = npar;
						Product[Row][Column++] = geometry->getNdatapoints();
						for (unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
							if (coeff < npar) {
								Product[Row][Column++] = polynomial->getCoefficient(coeff);
								Product[Row][Column++] = polynomial->getCoefficientError(coeff);
							} else {
								Product[Row][Column++] = 0.0;
								Product[Row][Column++] = 0.0;
							}
						}
						Product[Row][Column++] = polynomial->getChisqr();
						Product[Row][Column++] = geometry->getNumberofPointsToBinInYDirection();
						Product[Row][Column++] = geometry->getYmin();
						Product[Row][Column++] = geometry->getYmax();
						if (Column >= MaxColumns) {
							cout << "operaCreateProduct: ***Warning: GEOMETRY increase the column size to at least " << Column << endl;
						}
						Row++;
						Column = 0;
					}
				}
				product->saveExtension(extension);
			}
			if (!ordpfilename.empty()) {
				extension++;
				// <number of coefficients> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL ...
				operaSpectralOrderVector spectralOrders(ordpfilename);
				unsigned Rows = 1;
				unsigned Row = 0;
				unsigned Column = 0;
				MaxColumns = 2+MAXPOLYNOMIAL*2;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " order spacing polynomial " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("ORDERPOLY", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Order Polynomial in Dispersion Direction", extension);
				product->operaFITSAddComment("----------------------------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<number of coefficients> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				Polynomial *polynomial = spectralOrders.getOrderSpacingPolynomial();
				unsigned npar = polynomial->getOrderOfPolynomial();
				Product[Row][Column++] = npar;
				for (unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
					if (coeff < npar) {
						Product[Row][Column++] = polynomial->getCoefficient(coeff);
						Product[Row][Column++] = polynomial->getCoefficientError(coeff);
					} else {
						Product[Row][Column++] = 0.0;
						Product[Row][Column++] = 0.0;
					}
					if (Column >= MaxColumns) {
						cout << "operaCreateProduct: ***Warning: ORDERPOLY increase the column size to at least " << Column << endl;
					}
				}
				product->saveExtension(extension);
			}
			if (!dispfilename.empty()) {
				extension++;
				// <NumberOfDispersionPolynomials> <PolynomialIndex> <min coefficient> <max coefficient> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL ...
				operaSpectralOrderVector spectralOrders(dispfilename);
                unsigned Rows = spectralOrders.getnumberOfDispersionPolynomials();
				unsigned Row = 0;
				unsigned Column = 0;
				MaxColumns = 5+MAXPOLYNOMIAL*2;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " dispersion polynomial " << Columns << " x " << Rows << endl;
				}
				product->addExtension("DISPERSION", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Echelle Dispersion LaurentPolynomials", extension);
				product->operaFITSAddComment("-------------------------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<NumberOfDispersionPolynomials> <PolynomialIndex> <min coefficient> <max coefficient> <npar> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				
                for (unsigned PolynomialIndex = 0; PolynomialIndex < Rows; PolynomialIndex++) {
                    
                    LaurentPolynomial *polynomial = spectralOrders.getDispersionPolynomial(PolynomialIndex);
                    int minordcoeff = polynomial->getMinorderOfLaurentPolynomial();
                    int maxordcoeff = polynomial->getMaxorderOfLaurentPolynomial();
                    unsigned npar = polynomial->getNumberOfCoefficients();
                    
					Product[Row][Column++] = Rows;
					Product[Row][Column++] = PolynomialIndex;
                    Product[Row][Column++] = minordcoeff;
                    Product[Row][Column++] = maxordcoeff;
                    for (unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
                        if (coeff < npar) {
                            Product[Row][Column++] = polynomial->getCoefficient(coeff);
                            Product[Row][Column++] = polynomial->getCoefficientError(coeff);
                        } else {
                            Product[Row][Column++] = 0.0;
                            Product[Row][Column++] = 0.0;
                        }
						if (Column >= MaxColumns) {
							cout << "operaCreateProduct: ***Warning: DISPERSION increase the column size to at least " << Column << endl;
						}
                    }
					Column = 0;
					Row++;
                }
				product->saveExtension(extension);
			}
			if (!rvelfilename.empty()) {
				extension++;
				operaSpectralOrderVector spectralOrders(rvelfilename);
                float rvel = spectralOrders.getBarycentricRadialVelocityCorrection();
				unsigned Rows = 1;
				unsigned Row = 0;
				unsigned Column = 0;
				MaxColumns = 1;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " radial velocity " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("RADIALVELOCITY", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Barycentric Radial Velocity Correction (km/s)", extension);
				product->operaFITSAddComment("----------------------------------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<BarycentricRadialVelocityCorrection>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				
				Product[Row][Column] = rvel;
				product->saveExtension(extension);
			}
			if (!prvelfilename.empty()) {
				extension++;
				// <rvel> 
				operaSpectralOrderVector spectralOrders(prvelfilename);
                float rvel = spectralOrders.getBarycentricRadialVelocityCorrection();
				unsigned Rows = 1;
				unsigned Row = 0;
				unsigned Column = 0;
				MaxColumns = 1;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " polarimetry radial velocity " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("PRADIALVELOCITY", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Polarimetry Barycentric Radial Velocity Correction (km/s)", extension);
				product->operaFITSAddComment("----------------------------------------------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<BarycentricRadialVelocityCorrection>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				
				Product[Row][Column] = rvel;
				product->saveExtension(extension);
			}
			if (!tellfilename.empty()) {	// looks just like wcal
				extension++;
				// <orders><order number><number of coefficients><polynomial coefficient><polynomial coefficient error>...
				operaSpectralOrderVector spectralOrders(tellfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasWavelength()) {
						Rows++;
					}
				}
				MaxColumns = 4+MAXPOLYNOMIAL*2;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " telluric wavelength " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("TELLCORR", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Telluric Wavelength Polynomial for Intensity", extension);
				product->operaFITSAddComment("--------------------------------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><number of coefficients><polynomial coefficient><polynomial coefficient error>...", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasWavelength()) {
						Polynomial *polynomial = spectralOrder->getWavelength()->getWavelengthPolynomial() ;
						unsigned npar = polynomial->getOrderOfPolynomial();
						Product[Row][Column++] = maxorder - minorder + 1;
						Product[Row][Column++] = order;
						Product[Row][Column++] = npar;
						for (unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
							if (coeff < npar) {
								Product[Row][Column++] = polynomial->getCoefficient(coeff);
								Product[Row][Column++] = polynomial->getCoefficientError(coeff);
							} else {
								Product[Row][Column++] = 0.0;
								Product[Row][Column++] = 0.0;
							}
						}
						if (Column >= MaxColumns) {
							cout << "operaCreateProduct: ***Warning: TELLCORR increase the column size to at least " << Column << endl;
						}
						Row++;
						Column = 0;
					}
				}
				product->saveExtension(extension);
			}
			if (!ptellfilename.empty()) {	// looks just like wcal
				extension++;
				// <orders><order number><number of coefficients><polynomial coefficient><polynomial coefficient error>...
				operaSpectralOrderVector spectralOrders(ptellfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasWavelength()) {
						Rows++;
					}
				}
				MaxColumns = 4+MAXPOLYNOMIAL*2;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " polarimetry telluric wavelength " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("PTELLCORR", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Telluric Wavelength Polynomial for Polarimetry", extension);
				product->operaFITSAddComment("----------------------------------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><number of coefficients><polynomial coefficients>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasWavelength()) {
						Polynomial *polynomial = spectralOrder->getWavelength()->getWavelengthPolynomial() ;
						unsigned npar = polynomial->getOrderOfPolynomial();
						Product[Row][Column++] = maxorder - minorder + 1;
						Product[Row][Column++] = order;
						Product[Row][Column++] = npar;
						for (unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
							if (coeff < npar) {
								Product[Row][Column++] = polynomial->getCoefficient(coeff);
								Product[Row][Column++] = polynomial->getCoefficientError(coeff);
							} else {
								Product[Row][Column++] = 0.0;
								Product[Row][Column++] = 0.0;
							}
						}
						if (Column >= MaxColumns) {
							cout << "operaCreateProduct: ***Warning: PTELLCORR increase the column size to at least " << Column << endl;
						}
						Row++;
						Column = 0;
					}
				}
				product->saveExtension(extension);
			}
			if (!wavefilename.empty()) {
				extension++;
				// <order number><number of coefficients><polynomial coefficients>
				operaSpectralOrderVector spectralOrders(wavefilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasWavelength()) {
						Rows++;
					}
				}
				MaxColumns = 4+MAXPOLYNOMIAL*2;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " wavelength " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("WAVELENGTH", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Wavelength Polynomial", extension);
				product->operaFITSAddComment("---------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><number of coefficients><polynomial coefficients>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasWavelength()) {
						Polynomial *polynomial = spectralOrder->getWavelength()->getWavelengthPolynomial() ;
						unsigned npar = polynomial->getOrderOfPolynomial();
						Product[Row][Column++] = maxorder - minorder + 1;
						Product[Row][Column++] = order;
						Product[Row][Column++] = npar;
						for (unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
							if (coeff < npar) {
								Product[Row][Column++] = polynomial->getCoefficient(coeff);
								Product[Row][Column++] = polynomial->getCoefficientError(coeff);
							} else {
								Product[Row][Column++] = 0.0;
								Product[Row][Column++] = 0.0;
							}
						}
						if (Column >= MaxColumns) {
							cout << "operaCreateProduct: ***Warning: WAVELENGTH increase the column size to at least " << Column << endl;
						}
						Row++;
						Column = 0;
					}
				}
				product->saveExtension(extension);
			}
			if (!aperfilename.empty()) {
				extension++;
				// <order number><number of beams><measured tilt><tilt error><leftBackgroundIndex><xsampling><ysampling><lb height><lb width><lb slope><lb xcenter><lb ycenter><lb fluxFraction><rightBackgroundIndex><xsampling><ysampling><rb height><rb width><rb slope><rb xcenter><rb ycenter><rb fluxFraction><beam><xsampling><ysampling><beam height><beam width><beam slope><beam xcenter><beam ycenter><beam fluxFraction>
				operaSpectralOrderVector spectralOrders(aperfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasExtractionApertures()) {
						Rows++;
					}
				}
				MaxColumns = 43;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " aperture " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("APERTURE", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Extraction Aperture", extension);
				product->operaFITSAddComment("-------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><number of beams><measured tilt><tilt error><leftBackgroundIndex><xsampling><ysampling><lb height><lb width><lb slope><lb xcenter><lb ycenter><lb fluxFraction><rightBackgroundIndex><xsampling><ysampling><rb height><rb width><rb slope><rb xcenter><rb ycenter><rb fluxFraction><beam><xsampling><ysampling><beam height><beam width><beam slope><beam xcenter><beam ycenter><beam fluxFraction>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasExtractionApertures()) {
						Product[Row][Column++] = maxorder - minorder + 1;
						Product[Row][Column++] = order;
						Product[Row][Column++] = spectralOrder->getnumberOfBeams();
						Product[Row][Column++] = spectralOrder->getTiltInDegreesValue();
						Product[Row][Column++] = spectralOrder->getTiltInDegreesError();
						for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
							operaExtractionAperture *backgroundAperture = spectralOrder->getBackgroundApertures(backgroundIndex);
							Line *backgroundLineAperture = backgroundAperture->getLineAperture();
							Product[Row][Column++] = backgroundIndex;
							Product[Row][Column++] = backgroundAperture->getXsampling();
							Product[Row][Column++] = backgroundAperture->getYsampling();
							Product[Row][Column++] = backgroundLineAperture->getWidth();
							Product[Row][Column++] = backgroundLineAperture->getLength();
							Product[Row][Column++] = backgroundLineAperture->getSlope();
							Product[Row][Column++] = backgroundLineAperture->getMidPoint()->getXcoord();
							Product[Row][Column++] = backgroundLineAperture->getMidPoint()->getYcoord();
							Product[Row][Column++] = backgroundAperture->getFluxFraction();
						}
						for(unsigned beam=0; beam<spectralOrder->getnumberOfBeams(); beam++) {
							operaExtractionAperture *beamAperture = spectralOrder->getExtractionApertures(beam);
							Line *beamLineAperture = beamAperture->getLineAperture();
							Product[Row][Column++] = beam;
							Product[Row][Column++] = beamAperture->getXsampling();
							Product[Row][Column++] = beamAperture->getYsampling();
							Product[Row][Column++] = beamLineAperture->getWidth();
							Product[Row][Column++] = beamLineAperture->getLength();
							Product[Row][Column++] = beamLineAperture->getSlope();
							Product[Row][Column++] = beamLineAperture->getMidPoint()->getXcoord();
							Product[Row][Column++] = beamLineAperture->getMidPoint()->getYcoord();
							Product[Row][Column++] = beamAperture->getFluxFraction();
						}
						if (Column >= MaxColumns) {
							cout << "operaCreateProduct: ***Warning: APERTURE increase the column size to at least " << Column << endl;
						}
						Row++;
						Column = 0;
					}
				}
				product->saveExtension(extension);
			}
			if (!proffilename.empty()) {
				extension++;
				// <number of columns i><number of rows j><xsize><xsampling><ysize><ysampling><i><j><number of coefficients><ndatapoints><polynomial coefficients><chisqr>
				operaSpectralOrderVector spectralOrders(proffilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				unsigned Row = 0;
				unsigned Column = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasInstrumentProfile()) {
						operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
						Rows += instrumentProfile->getNYPoints()*instrumentProfile->getNXPoints();
					}
				}
				MaxColumns = 44;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " instrument profile " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("PROFILE", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Instrument Profile", extension);
				product->operaFITSAddComment("------------------", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("<orders><order><number of columns i><number of rows j><xsize><xsampling><ysize><ysampling><i><j><number of coefficients><ndatapoints><polynomial coefficients><chisqr>", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasInstrumentProfile()) {
						operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
						for (unsigned j=0; j<instrumentProfile->getNYPoints(); j++) {
							for (unsigned i=0; i<instrumentProfile->getNXPoints(); i++) {
								Product[Row][Column++] = maxorder - minorder + 1;
								Product[Row][Column++] = order;
								Product[Row][Column++] = instrumentProfile->getNXPoints();
								Product[Row][Column++] = instrumentProfile->getNYPoints();
								Product[Row][Column++] = instrumentProfile->getxsize();
								Product[Row][Column++] = instrumentProfile->getXsampling();
								Product[Row][Column++] = instrumentProfile->getysize();
								Product[Row][Column++] = instrumentProfile->getYsampling();
								Product[Row][Column++] = i;
								Product[Row][Column++] = j;
								PolynomialCoeffs_t *pp = instrumentProfile->getipPolyModelCoefficients(i,j);
								unsigned npar = pp->orderofPolynomial;
								Product[Row][Column++] = npar;
								Product[Row][Column++] = 0;//instrumentProfile->getnDataPoints();
								for	(unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
									if (coeff < npar) {
										Product[Row][Column++] = pp->p[coeff];
									} else {
										Product[Row][Column++] = 0.0;
									}
								}
								Product[Row][Column++] = instrumentProfile->getchisqrMatrixValue(i, j);
								if (Column >= MaxColumns) {
									cout << "operaCreateProduct: ***Warning: PROFILE increase the column size to at least " << Column << endl;
								}
								Row++;
								Column = 0;
							}
						}
					}
				}
				product->saveExtension(extension);
			}
			if (!fcalfilename.empty()) {
				extension++;
				// <order number><nElements><nBeams><elementindex><wavelength><SpectralElements flux conversion><flux conversion variance><SpectralElements throughput><throughput variance><nElements><nBeams><beam><BeamSED[beam] flux conversion><flux conversion variance><BeamSED[beam] throughput><throughput variance>
				operaSpectralOrderVector spectralOrders(fcalfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Rows = 0;
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasSpectralEnergyDistribution()) {
						operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
						if (spectralEnergyDistribution->getHasFluxCalibration() && spectralEnergyDistribution->getHasInstrumentThroughput()){
							operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();
							Rows += FluxCalibration->getnSpectralElements();
						}
					}
				}
				unsigned Row = 0;
				unsigned Column = 0;
				MaxColumns = 30;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " flux calibration " << MaxColumns << " x " << Rows << endl;
				}
				if (sequence == 0) {
					product->addExtension("FLUXCALIBRATION", MaxColumns, Rows, detsec);
				} else {
					product->addExtension("FLUXCALIBRATION"+itos(sequence), MaxColumns, Rows, detsec);
				}
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Flux Calibration", extension);
				product->operaFITSAddComment("----------------", extension);
				if (sequence > 0) {
					product->operaFITSAddComment("Sequence number "+itos(sequence)+" of 4", extension);
				}
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<orders><order number><nElements><nBeams><wavelengthForNormalization><elementindex><wavelength><SpectralElements flux conversion><flux conversion variance><SpectralElements throughput><throughput variance><nElements><nBeams><beam><BeamSED[beam] flux conversion><flux conversion variance><BeamSED[beam] throughput><throughput variance>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned order=minorder; order<=maxorder; order++) {
					operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
					if (spectralOrder->gethasSpectralEnergyDistribution()) {
						operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
						if (spectralEnergyDistribution->getHasFluxCalibration() && spectralEnergyDistribution->getHasInstrumentThroughput()){
							
							operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();
							operaSpectralElements *InstrumentThroughput = spectralEnergyDistribution->getThroughputElements();
							unsigned nElements = FluxCalibration->getnSpectralElements();
							unsigned nBeams = spectralOrder->getnumberOfBeams();
                            double wavelengthForNormalization = spectralEnergyDistribution->getwavelengthForNormalization();
							for (unsigned indexElem=0;indexElem < nElements; indexElem++) {
								Product[Row][Column++] = maxorder - minorder + 1;
								Product[Row][Column++] = order;
								Product[Row][Column++] = nElements;
								Product[Row][Column++] = nBeams;
								Product[Row][Column++] = wavelengthForNormalization;
								Product[Row][Column++] = indexElem;
								Product[Row][Column++] = FluxCalibration->getwavelength(indexElem);
								Product[Row][Column++] = FluxCalibration->getFlux(indexElem);
								Product[Row][Column++] = FluxCalibration->getFluxVariance(indexElem);
								Product[Row][Column++] = InstrumentThroughput->getFlux(indexElem);
								Product[Row][Column++] = InstrumentThroughput->getFluxVariance(indexElem);
								for (unsigned beam = 0; beam < nBeams; beam++) {
									Product[Row][Column++] = beam;
									Product[Row][Column++] = spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFlux(indexElem);
									Product[Row][Column++] = spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFluxVariance(indexElem);
									Product[Row][Column++] = spectralOrder->getBeamSED(beam)->getThroughputElements()->getFlux(indexElem);
									Product[Row][Column++] = spectralOrder->getBeamSED(beam)->getThroughputElements()->getFluxVariance(indexElem);
									if (Column >= MaxColumns) {
										cout << "operaCreateProduct: ***Warning: FLUXCALIBRATION increase the column size to at least " << Column << endl;
									}
								}
								Row++;
								Column = 0;
							}
						}
					}
				}
				product->saveExtension(extension);
			}
			if (!gainfilename.empty()) {
				// <amps><amp><gain><noise><gainerror><bias><amp><gain><noise><gainerror><bias><detsec x1><detsec x2><detsec y1><detsec y2>
				extension++;
				operaSpectralOrderVector spectralOrders(gainfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Row = 0;
				unsigned Column = 0;
				unsigned Rows = spectralOrders.getGainBiasNoise()->getAmps();
				MaxColumns = 14;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " gain " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("GAIN", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Gain", extension);
				product->operaFITSAddComment("----", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<amps><amp><gain><noise><gainerror><bias><amp><gain><noise><gainerror><bias><detsec x1><detsec x2><detsec y1><detsec y2>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned i=0; i<spectralOrders.getGainBiasNoise()->getAmps(); i++) {
					DATASEC_t datasec;
					spectralOrders.getGainBiasNoise()->getDatasec(i, datasec);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getAmps();
					Product[Row][Column++] = i;
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getGain(i);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getNoise(i);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getGainError(i);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getBias(i);
					Product[Row][Column++] = datasec.x1;
					Product[Row][Column++] = datasec.x2;
					Product[Row][Column++] = datasec.y1;
					Product[Row][Column++] = datasec.y2;
					if (Column >= Columns) {
						cout << "operaCreateProduct: ***Warning: GAIN increase the column size to at least " << Column << endl;
					}
					Row++;
					Column = 0;
				}
				product->saveExtension(extension);
			}
			if (!biasfilename.empty()) {
				// <amps><amp><gain><noise><gainerror><bias><amp><gain><noise><gainerror><bias><detsec x1><detsec x2><detsec y1><detsec y2>
				extension++;
				operaSpectralOrderVector spectralOrders(biasfilename);
				minorder = spectralOrders.getMinorder();
				maxorder = spectralOrders.getMaxorder();
				unsigned Row = 0;
				unsigned Column = 0;
				unsigned Rows = spectralOrders.getGainBiasNoise()->getAmps();
				MaxColumns = 14;
				detsec.x1 = detsec.x2 + 1;
				detsec.y2 = Rows;
                detsec.x2 = detsec.x1 + MaxColumns;
				if (args.verbose) {
					cout << "operaCreateProduct: adding extension " << extension << " bias " << MaxColumns << " x " << Rows << endl;
				}
				product->addExtension("BIAS", MaxColumns, Rows, detsec);
				operaFITSProduct Product(*product);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				product->operaFITSAddComment("Bias", extension);
				product->operaFITSAddComment("----", extension);
				product->operaFITSAddComment("Columns:", extension);
				product->operaFITSAddComment("<amps><amp><gain><noise><gainerror><bias><amp><gain><noise><gainerror><bias><detsec x1><detsec x2><detsec y1><detsec y2>", extension);
				product->operaFITSAddComment("------------------------------------------------------------------------", extension);
				for (unsigned i=0; i<spectralOrders.getGainBiasNoise()->getAmps(); i++) {
					DATASEC_t datasec;
					spectralOrders.getGainBiasNoise()->getDatasec(i, datasec);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getAmps();
					Product[Row][Column++] = i;
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getGain(i);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getNoise(i);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getGainError(i);
					Product[Row][Column++] = spectralOrders.getGainBiasNoise()->getBias(i);
					Product[Row][Column++] = datasec.x1;
					Product[Row][Column++] = datasec.x2;
					Product[Row][Column++] = datasec.y1;
					Product[Row][Column++] = datasec.y2;
					if (Column >= Columns) {
						cout << "operaCreateProduct: ***Warning: BIAS increase the column size to at least " << Column << endl;
					}
					Row++;
					Column = 0;
				}
				product->saveExtension(extension);
			}
			if (args.verbose) {
				cout << "operaCreateProduct: completed adding " << extension << " extensions." << endl;
			}
			product->operaFITSImageClose();
		}
    }
    catch (operaException e) {
        cerr << "operaCreateProduct: " << e.getFormattedMessage() << endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        cerr << "operaCreateProduct: " << operaStrError(errno) << endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

unsigned int GetRowCountFromLESpectrumFile(const string spectrumfile) {
	unsigned int rows = 0;
	operaistream fin(spectrumfile.c_str());
	if (fin.is_open()) {
		string dataline;
		getline(fin, dataline);
		getline(fin, dataline);
		stringstream ss (dataline);
		ss >> rows;
		fin.close();
	}
	return rows;
}

unsigned int GetRowCountFromExtendedSpectrumFile(const string spectrumfile, const instrumentmode_t instrumentmode) {
	operaSpectralOrderVector spectralOrders(spectrumfile);
	const unsigned int minorder = spectralOrders.getMinorder();
	const unsigned int maxorder = spectralOrders.getMaxorder();
	unsigned totaldatapoints = 0;
	for (unsigned int order = minorder; order <= maxorder; order++) {
		operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
		if (instrumentmode == MODE_POLAR && spectralOrder->gethasPolarimetry() && spectralOrder->gethasWavelength()) {
			totaldatapoints += spectralOrder->getPolarimetry()->getLength();
		}
		else if (instrumentmode != MODE_POLAR && spectralOrder->gethasSpectralElements()) {
			operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
			totaldatapoints += spectralelements->getnSpectralElements();
		}
	}
	return totaldatapoints;
}

void UpdateProductFromLESpectrumFile(operaFITSProduct& Product, const string spectrumfile, const unsigned groupsize, const unsigned int colgroup, const unsigned int rowcount) {
	operaistream fin(spectrumfile.c_str());
	if (fin.is_open()) {
		string dataline;
		getline(fin, dataline);
		getline(fin, dataline);
		unsigned row = 0;
		while (fin.good() && row < rowcount) {
			getline(fin, dataline);
			if (!dataline.empty()) {
				stringstream ss (dataline);
				for (unsigned col = colgroup * groupsize; col < (colgroup + 1) * groupsize; col++) {
					Float NanTolerantFloat = 0.0;
					ss >> NanTolerantFloat;
					Product[col][row] = NanTolerantFloat.f;
				}
				row++;
			}
		}
		fin.close();
	}	
}

void UpdateProductFromExtendedSpectrumFile(operaFITSProduct& Product, const string spectrumfile, const instrumentmode_t instrumentmode, const unsigned int colgroup) {
	operaSpectralOrderVector spectralOrders(spectrumfile);
	const unsigned int minorder = spectralOrders.getMinorder();
	const unsigned int maxorder = spectralOrders.getMaxorder();
	unsigned int startcol = 0;
	unsigned int datapoint = 0;
	switch (instrumentmode) {
		case MODE_POLAR:
			startcol = colgroup * MODE_POLAR_COLS / 4;
			for (unsigned order = minorder; order <= maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasPolarimetry() && spectralOrder->gethasWavelength()) {
					operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
					operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
					unsigned length = Polarimetry->getLength();
					stokes_parameter_t stokesParameter = StokesI;
					if (Polarimetry->getHasStokesV()) stokesParameter = StokesV;
					else if (Polarimetry->getHasStokesQ()) stokesParameter = StokesQ;
					else if (Polarimetry->getHasStokesU()) stokesParameter = StokesU;
					else if (Polarimetry->getHasStokesI()) stokesParameter = StokesI; //this is currently the default value anyway
					for (unsigned index = 0; index < length; index++) {
						Product[startcol+0][datapoint] = SpectralElements->getwavelength(index);
						Product[startcol+1][datapoint] = Polarimetry->getStokesParameter(stokesParameter)->getflux(index); //Polarization
						Product[startcol+2][datapoint] = stokesParameter;
						Product[startcol+3][datapoint] = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(index); //Null Spectrum 1
						Product[startcol+4][datapoint] = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(index); //Null Spectrum 2
						Product[startcol+5][datapoint] = sqrt(Polarimetry->getStokesParameter(stokesParameter)->getvariance(index)); //Polarization Error
						datapoint++;
					}
				}
			}
			break;
		case MODE_STAR_ONLY:
			startcol = colgroup * MODE_STAR_ONLY_COLS / 4;
			for (unsigned int order = minorder; order <= maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
					for (unsigned int i = 0; i < spectralelements->getnSpectralElements(); i++) {
						Product[startcol+0][datapoint] = spectralelements->getwavelength(i);
						Product[startcol+1][datapoint] = spectralelements->getFlux(i);
						Product[startcol+2][datapoint] = sqrt(spectralelements->getFluxVariance(i));
						datapoint++;
					}
				}
			}
			break;
		case MODE_STAR_PLUS_SKY:
			startcol = colgroup * MODE_STAR_PLUS_SKY_COLS / 4;
			for (unsigned order = minorder; order <= maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
					operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
					operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
					for (unsigned int i = 0; i < spectralelements->getnSpectralElements(); i++) {
						Product[startcol+0][datapoint] = spectralelements->getwavelength(i);
						Product[startcol+1][datapoint] = beamElements0->getFlux(i);
						Product[startcol+2][datapoint] = beamElements0->getFlux(i) + beamElements1->getFlux(i);
						Product[startcol+3][datapoint] = beamElements1->getFlux(i);
						Product[startcol+4][datapoint] = sqrt(beamElements0->getFluxVariance(i));
						Product[startcol+5][datapoint] = sqrt(beamElements0->getFluxVariance(i)+beamElements1->getFluxVariance(i));
						Product[startcol+6][datapoint] = sqrt(beamElements1->getFluxVariance(i));
						datapoint++;
					}
				}
			}
			break;
		default:
			throw operaException("operaCreateProduct: ", operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
	}
}
