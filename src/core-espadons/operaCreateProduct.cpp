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

#include <iomanip>
#include "libraries/gzstream.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaFITSProduct.h"
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

class WavelengthRange {
public:
	WavelengthRange() : wl0(0.0), wlf(0.0) {}
	WavelengthRange(double wl0, double wlf) : wl0(wl0), wlf(wlf) {}
	bool Contains(double wl) const { return wl0 <= wl && wl <= wlf; }
private:
	double wl0;
	double wlf;
};

class OrderBounds {
public:
	OrderBounds() {}
	OrderBounds(const string filename) { ReadFromFile(filename); }
	void ReadFromFile(const string filename);
	bool empty() const { return orderBounds.empty(); }
	bool OrderContains (unsigned order, double wl) const { return hasOrder(order) && atOrder(order).Contains(wl); }
private:
	typedef std::map<unsigned, WavelengthRange> OrderMap;
	bool hasOrder(unsigned order) const { return orderBounds.find(order) != orderBounds.end(); }
	WavelengthRange atOrder(unsigned order) const { OrderMap::const_iterator i = orderBounds.find(order); if (i != orderBounds.end()) return i->second; return WavelengthRange(); }
	OrderMap orderBounds;
};

// Reads in a table from a file into matrix, skipping the first skiplines of the file.
// orderBounds is an optional argument used when reading in spc files to trim the wavelength ranges of each order.
void GetMatrixFromDataFile(const string filename, vector<vector<float> > &matrix, const unsigned skiplines, const OrderBounds &orderBounds = OrderBounds());

// Updates the FITS product to contain the values in matrix. Starts at coloffset.
// Product must have at least as many rows as matrix and at least coloffset more columns than matrix.
void UpdateProductFromMatrix(operaFITSProduct& Product, const vector<vector<float> > &matrix, const unsigned coloffset = 0);

// Reads in the radial velocity correction from a rvel or tell file.
double readRadialVelocityCorrection(string filename);

// Adds various information to the header of the FITS product.
void AddFITSHeaderToProduct(operaFITSProduct& Product, const string version, const string date, const operaSpectralOrder_t spectralOrderType, const string snrfilename, const string rvelfilename, const string tellfilename);

int main(int argc, char *argv[])
{
	string version, date;
	string inputfilename;
	string outputfilename;
	string ufile, nfile, uwfile, nwfile;
	string spectrumfile;
	string wlrangefilename;
	string snrfilename;
	string rvelfilename;
	string tellfilename;
	string object;
	unsigned spectralOrderType_val = LibreEspritsp2Spectrum;
	int compressionVal;
	args.AddOptionalArgument("version", version, "", "");
	args.AddOptionalArgument("date", date, "", "");
	args.AddRequiredArgument("input", inputfilename, "input file (o.fits)");
	args.AddRequiredArgument("output", outputfilename, "output file (i.fits/m.fits)");
	args.AddOptionalArgument("ufile", ufile, "", "unnormalized spectrum (iu.s/pu.s)");
	args.AddOptionalArgument("nfile", nfile, "", "normalized spectrum (in.s/pn.s)");
	args.AddOptionalArgument("uwfile", uwfile, "", "unnormalized wavelength corrected spectrum (iuw.s/puw.s)");
	args.AddOptionalArgument("nwfile", nwfile, "", "normalized wavelength corrected spectrum (inw.s/pnw.s)");
	args.AddOptionalArgument("spectrumfile", spectrumfile, "", "extended spectrum (.spc)");
	args.AddOptionalArgument("wlrangefile", wlrangefilename, "", "LE order wavelength ranges");
	args.AddRequiredArgument("spectrumtype", spectralOrderType_val, "spectral order type");
	args.AddOptionalArgument("snr", snrfilename, "", ".sn");
	args.AddOptionalArgument("rvel", rvelfilename, "", "i.rvel");
	args.AddOptionalArgument("tell", tellfilename, "", "i.tell");
	args.AddOptionalArgument("object", object, "", "object name, needed for Libre-Esprit output");
	args.AddOptionalArgument("compressiontype", compressionVal, cNone, "compression type");
		
	try {
		args.Parse(argc, argv);
		
		operaSpectralOrder_t spectralOrderType = operaSpectralOrder_t(spectralOrderType_val);
		
		if (inputfilename.empty()) // we need an input...
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		if (outputfilename.empty()) // we need an output...
			throw operaException("operaCreateProduct: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		if ((ufile.empty() || nfile.empty() || uwfile.empty() || nwfile.empty()) && spectrumfile.empty()) // we need either a spectrum or all four LE files...
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		if ((!ufile.empty() || !nfile.empty() || !uwfile.empty() || !nwfile.empty()) && !spectrumfile.empty()) // ...but we don't wont both
			throw operaException("operaCreateProduct: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		
		eCompression compression = (eCompression)compressionVal;
		if (args.verbose) {
			cout << "operaCreateProduct: input= " << inputfilename << endl;
			cout << "operaCreateProduct: output= " << outputfilename << endl;
			cout << "operaCreateProduct: ufile= " << ufile << endl;
			cout << "operaCreateProduct: nfile= " << nfile << endl;
			cout << "operaCreateProduct: uwfile= " << uwfile << endl;
			cout << "operaCreateProduct: nwfile= " << nwfile << endl;
			cout << "operaCreateProduct: wlrangefile= " << wlrangefilename << endl;
			cout << "operaCreateProduct: snrfilename= " << snrfilename << endl;
			cout << "operaCreateProduct: rvelfilename= " << rvelfilename << endl;
			cout << "operaCreateProduct: tellfilename= " << tellfilename << endl;
			cout << "operaCreateProduct: OPERA version= " << version << endl;
			cout << "operaCreateProduct: Reduction date= " << date << endl;
			cout << "operaCreateProduct: compression= " << compression << endl;
			cout << "operaCreateProduct: spectrumtype= " << spectralOrderType << endl;
		}
		
		instrumentmode_t instrumentmode;
		switch(spectralOrderType) {
			case LibreEspritsp1Spectrum:
				instrumentmode = MODE_STAR_PLUS_SKY;
				break;
			case LibreEspritsp2Spectrum:
				instrumentmode = MODE_STAR_ONLY;
				break;
			case LibreEspritpolarimetry:
			case LibreEspritpolSpectrum:
				instrumentmode = MODE_POLAR;
				break;
			default:
				throw operaException("operaCreateProduct: ", operaErrorCodeBadInstrumentModeError, __FILE__, __FUNCTION__, __LINE__);
				break;
		}
		
		operaFITSProduct Product(outputfilename, inputfilename, instrumentmode, 1, 1, compression); //Start off 1x1, resize once we know the actual dimensions
		
		if (!ufile.empty() && !nfile.empty() && ! uwfile.empty() && !nwfile.empty()) {
			string inputfiles[4] = {nfile, ufile, nwfile, uwfile};
			for(int i = 0; i < 4; i++) {
				vector<vector<float> > readdata;
				GetMatrixFromDataFile(inputfiles[i], readdata, 2);
				if(i == 0) Product.resize(readdata.size(), readdata[0].size()*4);
				UpdateProductFromMatrix(Product, readdata, readdata[0].size()*i);
			}
        }
        
        else if (!spectrumfile.empty()) {
			OrderBounds orderBounds(wlrangefilename); //If filename is empty, this will not be used
			vector<vector<float> > readdata;
			GetMatrixFromDataFile(spectrumfile, readdata, 1, orderBounds);
			Product.resize(readdata.size(), readdata[0].size());
			UpdateProductFromMatrix(Product, readdata);
		}
		
		AddFITSHeaderToProduct(Product, version, date, spectralOrderType, snrfilename, rvelfilename, tellfilename);
		Product.operaFITSImageSave();
		Product.operaFITSImageClose();
		
		if (args.verbose && spectralOrderType == LibreEspritpolarimetry) cout << "operaCreateProduct: done polarimetry " << endl;
		else if (args.verbose) cout << "operaCreateProduct: done intensity " << endl;
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

void OrderBounds::ReadFromFile(const string filename) {
    if(filename.empty()) return;
    ifstream fin(filename.c_str());
    if (fin.is_open()) {
		string dataline;
        while (getline(fin, dataline)) {
            if (!dataline.empty() && dataline[0] != '#') { // skip comments
				istringstream ss(dataline);
				unsigned order;
				double wl0, wlf;
				ss >> order >> wl0 >> wlf;
				orderBounds[order] = WavelengthRange(wl0, wlf);
            }
		}
        fin.close();
    }
}

void GetMatrixFromDataFile(const string filename, vector<vector<float> > &matrix, const unsigned skiplines, const OrderBounds &orderBounds) {
	operaistream fin(filename.c_str());
	if (fin.is_open()) {
		string dataline;
		unsigned line = 0;
		while (getline(fin, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') {
				if (line >= skiplines) {
					istringstream ss (dataline);
					vector<float> datarow;
					for (Float NanTolerantFloat = 0.0; ss >> NanTolerantFloat; datarow.push_back(NanTolerantFloat.f));
					if (orderBounds.empty()) matrix.push_back(datarow);
					else if(datarow.size() > 4 && orderBounds.OrderContains(datarow[0], datarow[4])) matrix.push_back(datarow);
				}
				line++;
			}
		}
		fin.close();
	}
	if(matrix.empty()) matrix.push_back(vector<float>()); //to make sure we can call matrix[0].size()
}

void UpdateProductFromMatrix(operaFITSProduct& Product, const vector<vector<float> > &matrix, const unsigned coloffset) {
	for (unsigned row = 0; row < matrix.size(); row++) {
		for (unsigned col = 0; col < matrix[row].size(); col++) {
			Product[col+coloffset][row] = matrix[row][col];
		}
	}
}

double readRadialVelocityCorrection(string filename) {
	if (!filename.empty()) {
		operaistream fin(filename.c_str());
		if (fin.is_open()) {
			string dataline;
			while (getline (fin,dataline)) {
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

void AddFITSHeaderToProduct(operaFITSProduct& Product, const string version, const string date, const operaSpectralOrder_t spectralOrderType, const string snrfilename, const string rvelfilename, const string tellfilename) {
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
	Product.operaFITSAddComment(itos(spectralOrderType));
	//unused
	/*Product.operaFITSAddComment("OPERA Processing Parameters");
	Product.operaFITSAddComment("---------------------------");
	if (!parametersfilename.empty()) {
		if (args.verbose) cout << "operaCreateProduct: adding parameters " << endl;
		ifstream parameters(parametersfilename.c_str());
		while (parameters.good()) {
			string dataline;
			getline(parameters, dataline);
			if (strlen(dataline.c_str())) Product.operaFITSAddComment(dataline);
		}
	}*/
	//To do: replace this with more useful SNR information (i.e. peak SNR per ccd/spec bin)
	/*if (!snrfilename.empty()) {
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
	}*/
	double rvcorr = readRadialVelocityCorrection(rvelfilename);
	double tellcorr = readRadialVelocityCorrection(tellfilename);
	Product.operaFITSSetHeaderValue("HRV", rvcorr, "heliocentric RV correction");
	Product.operaFITSSetHeaderValue("TELLRV", tellcorr, "telluric RV correction");
}
