#include "libraries/operaIOFormats.h"
#include "libraries/gzstream.h"
#include "libraries/operastringstream.h"
#include <algorithm>
#include <iomanip>

class IOFormatFlags {
	public:
		IOFormatFlags(operaSpectralOrder_t format) : libreEsprit(false), singleval(false) {
			switch (format) {
				case GainNoise:
				case PRVel:
				case RVel:
				case Tell:
				case Orderspacing:
				case Disp:
					formattype = NONORDER;
					break;
				case Aperture:
				case Geom:
				case Wave:
				case SNR:
					formattype = OTHER;
					singleval = true;
					break;
				case Fcal:
				case Prof:
					formattype = OTHER;
					break;
				case RawSpectrum:
				case StandardSpectrum:
				case OptimalSpectrum:
				case OperaOptimalSpectrum:
				case RawBeamSpectrum:
				case StandardBeamSpectrum:
				case OptimalBeamSpectrum:
				case OperaOptimalBeamSpectrum:
					formattype = SPECTRUM;
					break;
				case CalibratedRawSpectrum:
				case CalibratedStandardSpectrum:
				case CalibratedOptimalSpectrum:
				case CalibratedOperaOptimalSpectrum:
				case CalibratedRawBeamSpectrum:
				case CalibratedStandardBeamSpectrum:
				case CalibratedOptimalBeamSpectrum:
				case CalibratedOperaOptimalBeamSpectrum:
				case CalibratedExtendedBeamSpectrum:
					formattype = CALIBRATED;
					break;
				case Polarimetry:
				case ExtendedPolarimetry:
					formattype = POLARIMETRY;
					break;
				case LibreEspritsp1Spectrum:
				case LibreEspritpolSpectrum:
				case LibreEspritsp2Spectrum:
					formattype = CALIBRATED;
					libreEsprit = true;
					break;
				case LibreEspritpolarimetry:
					formattype = POLARIMETRY;
					libreEsprit = true;
					break;
				case LibreEspritSNR:
					formattype = OTHER;
					libreEsprit = true;
					break;
				case CSV:
					formattype = OTHER;
					break;
				default:
					throw operaException("operaSpectralOrderVector: unkown IO format ", operaErrorCodeIncorrectFileType, __FILE__, __FUNCTION__, __LINE__);
					break;
			}
		}
		
		bool typeIsNonOrder() { return formattype == NONORDER; }
		bool typeIsSpectrum() { return formattype == SPECTRUM; }
		bool typeIsCalibrated() { return formattype == CALIBRATED; }
		bool typeIsPolarimetry() { return formattype == POLARIMETRY; }
		bool typeIsOther() { return formattype == OTHER; }
		bool isLE() { return libreEsprit; }
		bool isSingle() { return singleval; }
		bool isSNR() {return snr; }
	private:
		enum FormatType { NONORDER, SPECTRUM, CALIBRATED, POLARIMETRY, OTHER };
		FormatType formattype;
		bool libreEsprit;
		bool singleval;
		bool snr;
};

namespace operaIOFormats {
	const vector<string> GetFormatNames();
	const vector<string> formatnames = GetFormatNames();
	operaSpectralOrder_t FormatFromName(string name);
	string NameFromFormat(operaSpectralOrder_t format);
	/* Read/Write subroutines: */
	void writeFormatWithoutOrders(const operaSpectralOrderVector& orders, ostream &fout, operaSpectralOrder_t format);
	void readFormatWithoutOrders(operaSpectralOrderVector& orders, istream &fin, operaSpectralOrder_t format);
	void writeFormatWithOrders(const operaSpectralOrderVector& orders, ostream &fout, operaSpectralOrder_t format);
	void readFormatWithOrders(operaSpectralOrderVector& orders, istream &fin, operaSpectralOrder_t format);
	void writeLibreEsprit(const operaSpectralOrderVector& orders, ostream &fout, operaSpectralOrder_t format);
	void writeLibreEspritCenterSNR(const operaSpectralOrderVector& orders, ostream &fout);
	void writeCSVFromOrders(const operaSpectralOrderVector& orders, ostream& fout);
	void readOrdersFromCSV(operaSpectralOrderVector& orders, istream& fin);
	
	//Format selectors
	string getFormatHeader(operaSpectralOrder_t format);
	string getLinesFromFormatWithoutOrders(const operaSpectralOrderVector& orders, operaSpectralOrder_t format);
	string getLineFromFormatWithOrders(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format, unsigned index);
	void setFromLineWithoutOrders(operaSpectralOrderVector& orders, string dataline, unsigned linenumber, operaSpectralOrder_t format);
	void setFromLineWithOrders(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, operaSpectralOrder_t format, const std::vector<unsigned> &headervals);
	
	//Formats without orders
	string getGainNoiseHeader();
	string getOrderSpacingHeader();
	string getDispHeader();
	string getRvelHeader();
	string getTellHeader();
	//Formats with orders
	string getApertureHeader();
	string getFcalHeader();
	string getGeomHeader();
	string getWaveHeader();
	string getProfHeader();
	string getSNRHeader();
	string getSpectrumHeader();
	string getCalibratedSpectrumHeader();
	string getBeamSpectrumHeader();
	string getCalibratedBeamSpectrumHeader();
	string getCalibratedExtendedBeamSpectrumHeader();
	string getPolarimetryHeader();
	string getExtendedPolarimetryHeader();
	
	//Get all output lines for formats without individual orders
	string getLinesFromGainNoise(const operaSpectralOrderVector& orders);
	string getLinesFromOrderSpacing(const operaSpectralOrderVector& orders);
	string getLinesFromDisp(const operaSpectralOrderVector& orders);
	string getLinesFromRvel(const operaSpectralOrderVector& orders);
	string getLinesFromTell(const operaSpectralOrderVector& orders);
	//Get a single line from a spectral order format
	string getLineFromAperture(const operaSpectralOrder *spectralOrder);
	string getLineFromFcal(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromGeom(const operaSpectralOrder *spectralOrder);
	string getLineFromWave(const operaSpectralOrder *spectralOrder);
	string getLineFromProf(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromSNR(const operaSpectralOrder *spectralOrder);
	string getLineFromSpectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromBeamSpectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromCalibratedSpectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromCalibratedBeamSpectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromCalibratedExtendedBeamSpectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromPolarimetry(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromExtendedPolarimetry(const operaSpectralOrder *spectralOrder, unsigned index);
	//Libre-Esprit formats
	string getLineFromLibreEspritsp2Spectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromLibreEspritsp1Spectrum(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromLibreEspritpolarimetry(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromLibreEspritSNR(const operaSpectralOrder *spectralOrder, unsigned index);
	string getLineFromLibreEspritCenterSNR(const operaSpectralOrder *spectralOrder);
	
	//Formats without orders
	void setGainNoiseFromLine(operaSpectralOrderVector& orders, string dataline, unsigned linenumber);
	void setOrderSpacingFromLine(operaSpectralOrderVector& orders, string dataline);
	void setDispFromLine(operaSpectralOrderVector& orders, string dataline, unsigned linenumber);
	void setRvelFromLine(operaSpectralOrderVector& orders, string dataline);
	void setTellFromLine(operaSpectralOrderVector& orders, string dataline);
	//Formats with orders
	void setApertureFromLine(operaSpectralOrder *spectralOrder, string dataline);
	void setFcalFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder);
	void setGeomFromLine(operaSpectralOrder *spectralOrder, string dataline);
	void setWaveFromLine(operaSpectralOrder *spectralOrder, string dataline);
	void setProfFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned xsize, unsigned xsampling, unsigned ysize, unsigned ysampling);	
	void setSNRFromLine(operaSpectralOrder *spectralOrder, string dataline);
	void setSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, operaSpectralOrder_t format);
	void setCalibratedSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, operaSpectralOrder_t format);
	void setBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, operaSpectralOrder_t format);
	void setCalibratedBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, operaSpectralOrder_t format);
	void setCalibratedExtendedBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, operaSpectralOrder_t format);
	void setPolarFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, unsigned method);
	void setExtendedPolarimetryFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned stokespar, unsigned method);
	//Other order formats
	void setWavelengthRangeFromLine(operaSpectralOrder *spectralOrder, string dataline);
	
	// Helper functions
	stokes_parameter_t getStokesParameter(const operaSpectralOrder *spectralOrder);
	bool validFormatOrder(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format);
	bool validLibreEspritElement(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format, unsigned index);
	unsigned sizeOfFormatOrder(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format);
	bool centerSNROnly(const operaSpectralOrderVector& orders);
}

operaSpectralOrder_t operaIOFormats::FormatFromName(string name){
	vector<string>::const_iterator it = find(formatnames.begin(), formatnames.end(), name);
	if (it != formatnames.end()) return operaSpectralOrder_t(it-formatnames.begin());
	return None;
}
string operaIOFormats::NameFromFormat(operaSpectralOrder_t format) {
	unsigned formatindex = format;
	if(formatindex < formatnames.size()) return formatnames[formatindex];
	return "";
}

const vector<string> operaIOFormats::GetFormatNames() {
	vector<string> temp;
	temp.resize(count_SpectralOrder_t);
	temp[RawSpectrum] = "#!rawspectrum";
	temp[StandardSpectrum] = "#!standardspectrum";
	temp[OptimalSpectrum] = "#!optimalspectrum";
	temp[OperaOptimalSpectrum] = "#!operaoptimalspectrum";
	temp[RawBeamSpectrum] = "#!rawbeamspectrum";
	temp[StandardBeamSpectrum] = "#!standardbeamspectrum";
	temp[OptimalBeamSpectrum] = "#!optimalbeamspectrum";
	temp[OperaOptimalBeamSpectrum] = "#!operaoptimalbeamspectrum";
	temp[CalibratedRawSpectrum] = "#!calibratedrawspectrum";
	temp[CalibratedStandardSpectrum] = "#!calibratedstandardspectrum";
	temp[CalibratedOptimalSpectrum] = "#!calibratedoptimalspectrum";
	temp[CalibratedOperaOptimalSpectrum] = "#!calibratedoperaoptimalspectrum";
	temp[CalibratedRawBeamSpectrum] = "#!calibratedrawbeamspectrum";
	temp[CalibratedStandardBeamSpectrum] = "#!calibratedstandardbeamspectrum";
	temp[CalibratedOptimalBeamSpectrum] = "#!calibratedoptimalbeamspectrum";
	temp[CalibratedOperaOptimalBeamSpectrum] = "#!calibratedoperaoptimalbeamspectrum";
	temp[CalibratedExtendedBeamSpectrum] = "#!calibratedextendedbeamspectrum";
	temp[ExtendedPolarimetry] = "#!extendedpolarimetry";
	temp[Geom] = "#!geom";
	temp[Wave] = "#!wave";
	temp[Prof] = "#!prof";
	temp[Disp] = "#!disp";
	temp[SNR] = "#!SNR";
	temp[Polarimetry] = "#!polar";
	temp[Orderspacing] = "#!ordp";
	temp[GainNoise] = "#!gain";
	temp[Aperture] = "#!aper";
	temp[Fcal] = "#!fcal";
	temp[RVel] = "#!rvel";
	temp[Tell] = "#!tell";
	temp[CSV] = "#!csv";
	temp[OrderWavelengthRange] = "#!orderwlrange";
	return temp;
}

/* 
* void WriteSpectralOrders(string Filename, operaSpectralOrder_t Format), unsigned order=0, unsigned min=0;
* \brief Writes a SpectralOrder to a File
* \brief in the right place such that the vector remains ordered.
* \brief the optional order argument permits incremental addition to the output file, where zero means write all.
*/
void operaIOFormats::WriteFromSpectralOrders(const operaSpectralOrderVector& orders, string filename, operaSpectralOrder_t format) {
	operaostream fout(filename.c_str());
	if(fout.is_open()) {
		IOFormatFlags formatflags (format);
		if (format == CSV) writeCSVFromOrders(orders, fout);
		else if (format == LibreEspritSNR && centerSNROnly(orders)) writeLibreEspritCenterSNR(orders, fout);
		else if (formatflags.isLE()) writeLibreEsprit(orders, fout, format);
		else if (formatflags.typeIsNonOrder()) writeFormatWithoutOrders(orders, fout, format);
		else writeFormatWithOrders(orders, fout, format);
		fout.close();
	}
}

/* 
 * void ReadIntoSpectralOrders(string Filename);
 * \brief augment an existing vector with information from a file
 */
void operaIOFormats::ReadIntoSpectralOrders(operaSpectralOrderVector& orders, string filename) {
	operaSpectralOrder_t format = None;
	operaistream fin(filename.c_str());
	if (fin.is_open()) {
		string dataline;
		if (getline(fin, dataline)) {
			format = FormatFromName(dataline);
		}
		//
		switch (format) {
			case None:
				throw operaException("operaIOFormats: unkown content type in "+filename+' ', operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
				break;
			case CSV:
				readOrdersFromCSV(orders, fin);
				break;
			case GainNoise:
			case Orderspacing:
			case Disp:
			case RVel:
			case Tell:
				readFormatWithoutOrders(orders, fin, format);
				break;
			default:
				readFormatWithOrders(orders, fin, format);
				break;
		}
		fin.close();
	}
	else throw operaException("operaIOFormats: could not open file "+filename+' ', operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
	if(orders.getMaxorder() != 0) orders.setCount(orders.getMaxorder() - orders.getMinorder() + 1);
}

void operaIOFormats::writeFormatWithoutOrders(const operaSpectralOrderVector& orders, ostream &fout, operaSpectralOrder_t format) {
	fout << formatnames[format] << endl;
	fout << "######################################################################" << endl;
	fout << getFormatHeader(format);
	fout << "#" << endl << "######################################################################" << endl;
	fout << getLinesFromFormatWithoutOrders(orders, format);
}

/* 
 * void readFormatWithoutOrders(string filename, operaSpectralOrder_t format)
 * \brief Reads in data that is not part of any particular spectral order.
 */
void operaIOFormats::readFormatWithoutOrders(operaSpectralOrderVector& orders, istream &fin, operaSpectralOrder_t format) {
	string dataline;
	unsigned line = 0;
	while (getline(fin, dataline)) {
		if (!dataline.empty() && dataline[0] != '#') {
			setFromLineWithoutOrders(orders, dataline, line, format);
			line++;
		}
	}
}

void operaIOFormats::writeFormatWithOrders(const operaSpectralOrderVector& orders, ostream &fout, operaSpectralOrder_t format) {
	fout << formatnames[format] << endl;
	fout << "######################################################################" << endl;
	fout << getFormatHeader(format);
	fout << "#" << endl << "######################################################################" << endl;
	fout << orders.getMaxorder() - orders.getMinorder() + 1;
	const operaSpectralOrder *firstValidOrder = NULL;
	for(unsigned order = orders.getMinorder(); order <= orders.getMaxorder(); order++) {
		firstValidOrder = orders.GetSpectralOrder(order);
		if (validFormatOrder(firstValidOrder, format)) break;
		else firstValidOrder = NULL;
	}
	if (format == Prof && firstValidOrder) {
		const operaInstrumentProfile *ip = firstValidOrder->getInstrumentProfile();
		fout << ' ' << ip->getNXPoints() << ' ' << ip->getNYPoints() << ' ' << ip->getxsize() << ' ' << ip->getXsampling() << ' ' << ip->getysize() << ' ' << ip->getYsampling();
	}
	else if (format == Polarimetry && firstValidOrder) fout << ' ' << 14 << ' ' << firstValidOrder->getPolarimetry()->getmethod();
	else if (format == ExtendedPolarimetry && firstValidOrder) fout << ' ' << getStokesParameter(firstValidOrder) << ' ' << firstValidOrder->getPolarimetry()->getmethod();
	fout << endl;
	for(unsigned order = orders.getMinorder(); order <= orders.getMaxorder(); order++) {
		const operaSpectralOrder *spectralOrder = orders.GetSpectralOrder(order);
		if (validFormatOrder(spectralOrder, format)) {
			unsigned length = sizeOfFormatOrder(spectralOrder, format);
			for (unsigned index = 0; index < length; index++) {
				fout << getLineFromFormatWithOrders(spectralOrder, format, index) << endl;
			}
		}
	}
}

/* 
 * void readFormatWithOrders(string filename, operaSpectralOrder_t format)
 * \brief Reads in a file and creates or updates the spectral orders.
 */
void operaIOFormats::readFormatWithOrders(operaSpectralOrderVector& orders, istream &fin, operaSpectralOrder_t format) {
	unsigned line = 0;
	unsigned lastorder = 0;
	unsigned index = 0;
	std::vector<unsigned> headervals;
	string dataline;
	while (getline(fin, dataline)) {
		if (!dataline.empty() && dataline[0] != '#') {
			istringstream ss(dataline);
			if (line == 0) {
				unsigned temp;
				while (ss >> temp) headervals.push_back(temp);
				line++;
			} else {
				unsigned order;
				ss >> order;
#ifdef RANGE_CHECK
				if (order > length) throw operaException("operaIOFormats: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
#endif
				if (order < orders.getMinorder() || orders.getMinorder() == 0) orders.setMinorder(order);
				if (order > orders.getMaxorder() || orders.getMaxorder() == 0) orders.setMaxorder(order);
				if (lastorder != order) index = 0;
				operaSpectralOrder *spectralOrder = orders.GetSpectralOrder(order);
				setFromLineWithOrders(spectralOrder, dataline, lastorder!=order, index, format, headervals);
				lastorder = order;
				index++;
				line++;
			}
		}
	}
}

void operaIOFormats::writeLibreEsprit(const operaSpectralOrderVector& orders, ostream &fout, operaSpectralOrder_t format) {
	unsigned rows = 0;
	for (unsigned order=orders.getMaxorder(); order>=orders.getMinorder(); order--) {
		const operaSpectralOrder *spectralOrder = orders.GetSpectralOrder(order);
		if (validFormatOrder(spectralOrder, format)) {
			unsigned length = sizeOfFormatOrder(spectralOrder, format);
			for (unsigned index = 0; index < length; index++) {
				if (validLibreEspritElement(spectralOrder, format, index)) {
					rows++;
				}
			}
		}
	}
	unsigned cols = 0;
	if (format == LibreEspritsp2Spectrum || format == LibreEspritpolSpectrum) cols = 2;
	else if (format == LibreEspritsp1Spectrum) cols = 6;
	else if (format == LibreEspritpolarimetry) cols = 5;
	else if (format == LibreEspritSNR) cols = 1;
	if(format == LibreEspritSNR) fout << "***SNR of '" << orders.getObject() << "'" << endl;
	else fout << "***Reduced spectrum of '" << orders.getObject() << "'" << endl;
	fout << rows << " " << cols << endl;
	for (unsigned order=orders.getMaxorder(); order>=orders.getMinorder(); order--) {
		const operaSpectralOrder *spectralOrder = orders.GetSpectralOrder(order);
		if (validFormatOrder(spectralOrder, format)) {
			unsigned count = 0;
			unsigned length = sizeOfFormatOrder(spectralOrder, format);
			for (unsigned index = 0; index < length; index++) {
				if (validLibreEspritElement(spectralOrder, format, index)) {
					fout << getLineFromFormatWithOrders(spectralOrder, format, index) << endl;
					count++;
				}
			}
			if (NEWLINES_BETWEEN_ORDERS && count > 0) fout << endl; // split the orders for plotting
		}
	}
}

void operaIOFormats::writeLibreEspritCenterSNR(const operaSpectralOrderVector& orders, ostream &fout) {
	fout << "***SNR of '" << orders.getObject() << "'" << endl;
	fout << orders.getMaxorder() - orders.getMinorder() + 1 << " 1" << endl;
	for (unsigned order=orders.getMaxorder(); order>=orders.getMinorder(); order--) {
		const operaSpectralOrder *spectralOrder = orders.GetSpectralOrder(order);
		if (validFormatOrder(spectralOrder, SNR)) {
			fout << getLineFromLibreEspritCenterSNR(spectralOrder) << endl;
		}
	}
}

bool operaIOFormats::centerSNROnly(const operaSpectralOrderVector& orders) {
	for(unsigned order = orders.getMinorder(); order <= orders.getMaxorder(); order++) {
		if (orders.GetSpectralOrder(order)->gethasCenterSNROnly()) return true;
	}
	return false;
}

bool operaIOFormats::validFormatOrder(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format) {
	IOFormatFlags formatflags (format);
	if (formatflags.typeIsPolarimetry()) return spectralOrder->gethasSpectralElements() && spectralOrder->getSpectralElements()->getHasWavelength() && spectralOrder->gethasPolarimetry();
	if (formatflags.typeIsCalibrated() || formatflags.isLE()) return spectralOrder->gethasSpectralElements() && spectralOrder->getSpectralElements()->getHasWavelength(); //Calibrated spectra, LE spectra, LE SNR
	if (formatflags.typeIsSpectrum()) return spectralOrder->gethasSpectralElements(); //Uncalibrated spectra
	switch (format) {
		case Aperture: return spectralOrder->gethasExtractionApertures();
		case Fcal: return spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->getSpectralEnergyDistribution()->getHasFluxCalibration() && spectralOrder->getSpectralEnergyDistribution()->getHasInstrumentThroughput();
		case Geom: return spectralOrder->gethasGeometry();
		case Wave: return spectralOrder->gethasWavelength();
		case Prof: return spectralOrder->gethasInstrumentProfile();
		case SNR: return spectralOrder->gethasSpectralElements() && spectralOrder->getSpectralElements()->getHasWavelength() && !isnan(spectralOrder->getCenterSNR()) && spectralOrder->gethasCenterSNROnly();
		default: return false;
	}
}

unsigned operaIOFormats::sizeOfFormatOrder(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format) {
	IOFormatFlags formatflags (format);
	if (formatflags.isSingle()) return 1;
	if (formatflags.typeIsPolarimetry()) return spectralOrder->getPolarimetry()->getLength();
	if (format == Fcal) return spectralOrder->getSpectralEnergyDistribution()->getFluxCalibrationElements()->getnSpectralElements();
	if (format == Prof) return spectralOrder->getInstrumentProfile()->getNYPoints() * spectralOrder->getInstrumentProfile()->getNXPoints();
	if (formatflags.typeIsSpectrum() || formatflags.typeIsCalibrated() || formatflags.isLE()) return spectralOrder->getSpectralElements()->getnSpectralElements();
	return 0;
}

stokes_parameter_t operaIOFormats::getStokesParameter(const operaSpectralOrder *spectralOrder) {
	const operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
	if (Polarimetry->getHasStokesV()) return StokesV;
	if (Polarimetry->getHasStokesQ()) return StokesQ;
	if (Polarimetry->getHasStokesU()) return StokesU;
	return StokesI;
}

bool operaIOFormats::validLibreEspritElement(const operaSpectralOrder *spectralOrder, operaSpectralOrder_t format, unsigned index) {
	switch (format) {
		case LibreEspritpolSpectrum:
		case LibreEspritsp2Spectrum: {
			const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
			return !isnan(spectralElements->getFlux(index)) && !isnan(spectralElements->getFluxVariance(index));
		}
		case LibreEspritsp1Spectrum: {
			const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
			const operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
			const operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
			return !isnan(spectralElements->getFlux(index)) && !isnan(beamElements1->getFlux(index)) && !isnan(beamElements0->getFluxVariance(index)) && !isnan(beamElements1->getFluxVariance(index)); //?
		}
		case LibreEspritpolarimetry: {
			const operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
			stokes_parameter_t stokesParameter = getStokesParameter(spectralOrder);
			return !isnan(Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(index)) && !isnan(Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(index))
				&& !isnan(Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(index)) && !isnan(Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(index));
		}
		case LibreEspritSNR:
			return !isnan(spectralOrder->getSpectralElements()->getFluxSNR(index));
		default:
			return false;
	}
}

//This doesn't look like it would work currently...
void operaIOFormats::readOrdersFromCSV(operaSpectralOrderVector& orders, istream& flist) {
	string dataline;
	unsigned line = 0;
	unsigned order = 0;
	unsigned lastorder = 0;
    char Object[48];
    char Mode[24];
	float Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance,DegreeOfPolarization,StokesParameter,NullSpectrum1,NullSpectrum2,PolarizationError;

	//operaSpectralOrder *spectralOrder = NULL;
	while (getline (flist,dataline)) {
		if (!dataline.empty()) {
			if (dataline[0] == '#') {
				sscanf(dataline.c_str(), ",%s", Object);
			} else {
				sscanf(dataline.c_str(), ",,%s", Mode);
				// pol
				// #!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance,DegreeOfPolarization,StokesParameter,NullSpectrum1,NullSpectrum2,PolarizationError
				// #!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance
				if (startsWith(Mode, "polar")) {
					sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR, &LeftBeamFlux, &RightBeamFlux, &LeftBeamFluxVariance, &RightBeamFluxvariance, &DegreeOfPolarization, &StokesParameter, &NullSpectrum1, &NullSpectrum2, &PolarizationError);
				} else if (startsWith(Mode, "pol")) {
					sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR, &LeftBeamFlux, &RightBeamFlux, &LeftBeamFluxVariance, &RightBeamFluxvariance);
				}
				// star plus sky
				// #!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance
				if (startsWith(Mode, "sp1")) {
					sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR, &LeftBeamFlux, &RightBeamFlux, &LeftBeamFluxVariance, &RightBeamFluxvariance);
				}
				// star only
				// #!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR
				if (startsWith(Mode, "sp2")) {
					sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR);
				}
				if (order != lastorder) {
					//spectralOrder = GetSpectralOrder(order);
				}
#ifdef RANGE_CHECK
				if (order > length) {
					throw operaException("operaIOFormats: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
				}
#endif
				if (order < orders.getMinorder() || orders.getMinorder() == 0) {
					orders.setMinorder(order);
				}
				if (order > orders.getMaxorder() || orders.getMaxorder() == 0) {
					orders.setMaxorder(order);
				}
				lastorder = order;
				line++;
			}
		}
	}
}

void operaIOFormats::writeCSVFromOrders(const operaSpectralOrderVector& orders, ostream& fout) {
	if (fout.good()) {
		for(unsigned order = orders.getMinorder(); order <= orders.getMaxorder(); order++) {
			const operaSpectralOrder *spectralOrder = orders.GetSpectralOrder(order);
			switch (orders.getInstrumentmode()) {
				case MODE_POLAR:
					if (spectralOrder->gethasPolarimetry()) {
						fout << "#!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance,DegreeOfPolarization,StokesParameter,NullSpectrum1,NullSpectrum2,PolarizationError" << endl;
					} else {
						fout << "#!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance" << endl;
					}
					for (unsigned order=orders.getMinorder(); order<=orders.getMaxorder(); order++) {
						if (spectralOrder->gethasSpectralElements()) {
							const operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								
								if (spectralOrder->gethasPolarimetry()) {
									double PolarizationVariance, Polarization, NullSpectrum1, NullSpectrum2;
									const operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
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
									Polarization = Polarimetry->getStokesParameter(stokesParameter)->getflux(i);
									PolarizationVariance = Polarimetry->getStokesParameter(stokesParameter)->getvariance(i);
									NullSpectrum1 = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(i);
									NullSpectrum2 = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(i);
									fout << ',' 
									<< ",,polar" << ','
									<< spectralelements->getwavelength(i) << ','
									<< spectralelements->getFlux(i) << ','
									<< spectralelements->getFlux(i) << ',' // fixme
									<< sqrt(spectralelements->getFluxVariance(i)) << ','
									<< spectralelements->getFluxSNR(i) << ','
									<< Polarization << ','
									<< stokesParameter << ','
									<< NullSpectrum1 << ','
									<< NullSpectrum2 << ','
									<< sqrt(PolarizationVariance)
									<< endl;
								} else {	// gethasPolarimetry()
									fout << ",,pol" << ','
									<< spectralelements->getwavelength(i) << ','
									<< spectralelements->getFlux(i) << ','
									<< spectralelements->getFlux(i) << ',' // fixme
									<< sqrt(spectralelements->getFluxVariance(i)) << ','
									<< spectralelements->getFluxSNR(i)
									<< endl;
								}
							}
						}
					}
					break;
				case MODE_STAR_ONLY:
					fout << "#!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR" << endl;
					for (unsigned order=orders.getMinorder(); order<=orders.getMaxorder(); order++) {
						if (spectralOrder->gethasSpectralElements()) {
							const operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								fout << ",,sp2" << ','
								<< spectralelements->getwavelength(i) << ','
								<< spectralelements->getFlux(i) << ','
								<< spectralelements->getFlux(i) << ',' // fixme
								<< sqrt(spectralelements->getFluxVariance(i)) << ','
								<< spectralelements->getFluxSNR(i)
								<< endl;
							}
						}
					}
					break;
				case MODE_STAR_PLUS_SKY:
					fout << "#!csv,"+orders.getObject()+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance" << endl;
					for (unsigned order=orders.getMinorder(); order<=orders.getMaxorder(); order++) {
						const operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
						const operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
						if (spectralOrder->gethasSpectralElements()) {
							const operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								fout << ",,sp1" << ','
								<< spectralelements->getwavelength(i) << ','
								<< (beamElements0->getFlux(i) + beamElements1->getFlux(i)) << ','
								<< (beamElements0->getFlux(i) + beamElements1->getFlux(i)) << ','// fixme
								<< sqrt(beamElements0->getFluxVariance(i)+beamElements1->getFluxVariance(i)) << ','
								<< spectralelements->getFluxSNR(i) << ','
								<< beamElements0->getFlux(i) << ','
								<< beamElements1->getFlux(i) << ','
								<< sqrt(beamElements0->getFluxVariance(i)) << ','
								<< sqrt(beamElements1->getFluxVariance(i)) << endl;
							}
						}
					}
					break;
				default:
					break;
			}
		}
	}
}

string operaIOFormats::getFormatHeader(operaSpectralOrder_t format) {
	switch (format) {
		case GainNoise: return getGainNoiseHeader();
		case Orderspacing: return getOrderSpacingHeader();
		case Disp: return getDispHeader();
		case RVel: return getRvelHeader();
		case Tell: return getTellHeader();
		case Aperture: return getApertureHeader();
		case Fcal: return getFcalHeader();
		case Geom: return getGeomHeader();
		case Wave: return getWaveHeader();
		case Prof: return getProfHeader();
		case SNR: return getSNRHeader();
		case RawSpectrum:
		case StandardSpectrum:
		case OptimalSpectrum:
		case OperaOptimalSpectrum:
			return getSpectrumHeader();
		case CalibratedRawSpectrum:
		case CalibratedStandardSpectrum:
		case CalibratedOptimalSpectrum:
		case CalibratedOperaOptimalSpectrum:
			return getCalibratedSpectrumHeader();
		case RawBeamSpectrum:
		case StandardBeamSpectrum:
		case OptimalBeamSpectrum:
		case OperaOptimalBeamSpectrum:
			return getBeamSpectrumHeader();
		case CalibratedRawBeamSpectrum:
		case CalibratedStandardBeamSpectrum:
		case CalibratedOptimalBeamSpectrum:
		case CalibratedOperaOptimalBeamSpectrum:
			return getCalibratedBeamSpectrumHeader();
		case CalibratedExtendedBeamSpectrum: return getCalibratedExtendedBeamSpectrumHeader();
		case Polarimetry: return getPolarimetryHeader();
		case ExtendedPolarimetry: return getExtendedPolarimetryHeader();
		default: break;
	}
}

string operaIOFormats::getLinesFromFormatWithoutOrders(const operaSpectralOrderVector& orders, operaSpectralOrder_t format) {
	switch (format) {
		case GainNoise: return getLinesFromGainNoise(orders);
		case Orderspacing: return getLinesFromOrderSpacing(orders);
		case Disp: return getLinesFromDisp(orders);
		case RVel: return getLinesFromRvel(orders);
		case Tell: return getLinesFromTell(orders);
		default: return string();
	}
}

string operaIOFormats::getLineFromFormatWithOrders(const operaSpectralOrder* spectralOrder, operaSpectralOrder_t format, unsigned index) {
	switch (format) {
		case Aperture: return getLineFromAperture(spectralOrder);
		case Fcal: return getLineFromFcal(spectralOrder, index);
		case Geom: return getLineFromGeom(spectralOrder);
		case Wave: return getLineFromWave(spectralOrder);
		case Prof: return getLineFromProf(spectralOrder, index);
		case SNR: return getLineFromSNR(spectralOrder);
		case RawSpectrum:
		case StandardSpectrum:
		case OptimalSpectrum:
		case OperaOptimalSpectrum:
			return getLineFromSpectrum(spectralOrder, index);
		case CalibratedRawSpectrum:
		case CalibratedStandardSpectrum:
		case CalibratedOptimalSpectrum:
		case CalibratedOperaOptimalSpectrum:
			return getLineFromCalibratedSpectrum(spectralOrder, index);
		case RawBeamSpectrum:
		case StandardBeamSpectrum:
		case OptimalBeamSpectrum:
		case OperaOptimalBeamSpectrum:
			return getLineFromBeamSpectrum(spectralOrder, index);
		case CalibratedRawBeamSpectrum:
		case CalibratedStandardBeamSpectrum:
		case CalibratedOptimalBeamSpectrum:
		case CalibratedOperaOptimalBeamSpectrum:
			return getLineFromCalibratedBeamSpectrum(spectralOrder, index);
		case CalibratedExtendedBeamSpectrum: return getLineFromCalibratedExtendedBeamSpectrum(spectralOrder, index);
		case Polarimetry: return getLineFromPolarimetry(spectralOrder, index);
		case ExtendedPolarimetry: return getLineFromExtendedPolarimetry(spectralOrder, index);
		case LibreEspritpolSpectrum:
		case LibreEspritsp2Spectrum:
			return getLineFromLibreEspritsp2Spectrum(spectralOrder, index);
		case LibreEspritsp1Spectrum: return getLineFromLibreEspritsp1Spectrum(spectralOrder, index);
		case LibreEspritpolarimetry: return getLineFromLibreEspritpolarimetry(spectralOrder, index);
		case LibreEspritSNR: return getLineFromLibreEspritSNR(spectralOrder, index);
		default: return string();
	}
}

/* 
 * void updateFromLineWithoutOrder(string dataline, unsigned linenumber, operaSpectralOrder_t format)
 * \brief Calls the appropriate function, depending on the format, to update from the line.
 */
void operaIOFormats::setFromLineWithoutOrders(operaSpectralOrderVector& orders, string dataline, unsigned linenumber, operaSpectralOrder_t format) {
	switch (format) {
		case GainNoise: return setGainNoiseFromLine(orders, dataline, line);
		case Orderspacing: return setOrderSpacingFromLine(orders, dataline);
		case Disp: return setDispFromLine(orders, dataline, line);
		case RVel: return setRvelFromLine(orders, dataline);
		case Tell: return setTellFromLine(orders, dataline);
		default: return;
	}
}

/* 
 * void updateFromLineWithOrder(string dataline, unsigned index, operaSpectralOrder_t format, const std::vector<unsigned> &headervals)
 * \brief Calls the appropriate function, depending on the format, to update from the line
 * \brief headervals contains parameters already read in from the first line, that may be used depending on the format.
 */
void operaIOFormats::setFromLineWithOrders(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, operaSpectralOrder_t format, const std::vector<unsigned> &headervals) {
	switch (format) {
		case Aperture: return setApertureFromLine(spectralOrder, dataline);
		case Fcal: return setFcalFromLine(spectralOrder, dataline, neworder);
		case Geom: return setGeomFromLine(spectralOrder, dataline);
		case Wave: return setWaveFromLine(spectralOrder, dataline);
		case Prof: return setProfFromLine(spectralOrder, dataline, neworder, headervals[3], headervals[4], headervals[5], headervals[6]); //xsize, xsampling, ysize, ysampling
		case SNR: return setSNRFromLine(spectralOrder, dataline);
		case RawSpectrum:
		case StandardSpectrum:
		case OptimalSpectrum:
		case OperaOptimalSpectrum:
			return setSpectrumFromLine(spectralOrder, dataline, neworder, index, format);
		case RawBeamSpectrum:
		case StandardBeamSpectrum:
		case OptimalBeamSpectrum:
		case OperaOptimalBeamSpectrum:
			return setBeamSpectrumFromLine(spectralOrder, dataline, neworder, format);
		case CalibratedRawSpectrum:
		case CalibratedStandardSpectrum:
		case CalibratedOptimalSpectrum:
		case CalibratedOperaOptimalSpectrum:
			return setCalibratedSpectrumFromLine(spectralOrder, dataline, neworder, index, format);
		case CalibratedRawBeamSpectrum:
		case CalibratedStandardBeamSpectrum:
		case CalibratedOptimalBeamSpectrum:
		case CalibratedOperaOptimalBeamSpectrum:
			return setCalibratedBeamSpectrumFromLine(spectralOrder, dataline, neworder, format);
		case CalibratedExtendedBeamSpectrum: return setCalibratedExtendedBeamSpectrumFromLine(spectralOrder, dataline, neworder, format);
		case Polarimetry: return setPolarFromLine(spectralOrder, dataline, neworder, index, headervals[2]); //method
		case ExtendedPolarimetry: return setExtendedPolarimetryFromLine(spectralOrder, dataline, neworder, headervals[1], headervals[2]); //StokesParameter, method
		case OrderWavelengthRange: return setWavelengthRangeFromLine(spectralOrder, dataline);
		default: return;
	}
}
				
string operaIOFormats::getGainNoiseHeader() {
	ostringstream ss;
	ss << "# Gain Noise format is:\n"
	<< "# <number of amps> <newline>\n"
	<< "# <amp> <gain> <noise> <gainerror> <bias> <datasec x1> <datasec x2> <datasec y1> <datasec y2> <newline>\n"
	<< "# ...\n"
	<< "# Note that <amp> is zero-based.\n"
	<< "# Note that gain is in units e/ADU, noise in e and bias in ADU.\n";
	return ss.str();
}

string operaIOFormats::getLinesFromGainNoise(const operaSpectralOrderVector& orders) {
	ostringstream ss;
	const GainBiasNoise* gbn = orders.getGainBiasNoise();
	ss << gbn->getAmps() << endl;
	for (unsigned i = 0; i < gbn->getAmps(); i++) {
		DATASEC_t datasec;
		gbn->getDatasec(i, datasec);
		ss << i << ' ' << gbn->getGain(i) << ' '  << gbn->getNoise(i) << ' ' << gbn->getGainError(i)  << ' ' << gbn->getBias(i) << ' '
		<< datasec.x1  << ' ' << datasec.x2  << ' ' << datasec.y1  << ' ' << datasec.y2 << endl;
	}
	return ss.str();
}

/* 
 * void setGainNoiseFromLine(string dataline, unsigned linenumber)
 * \brief Sets gain, gain error, noise, and bias for a given amp from a line read in from a gain file.
 * \brief First line contains total number of amps.
 */
void operaIOFormats::setGainNoiseFromLine(operaSpectralOrderVector& orders, string dataline, unsigned linenumber) {
	istringstream ss(dataline);
	GainBiasNoise* gbn = orders.getGainBiasNoise();
	if (linenumber == 0) {
		unsigned ampcount = 0;
		ss >> ampcount;
		gbn->setAmps(ampcount);
	} else {
		unsigned amp = 0, x1 = 0, y1 = 0, x2 = 0, y2 = 0;
		Double gain = 0.0, gainerror = 0.0, noise = 0.0, bias = 0.0;
		ss >> amp >> gain >> noise >> gainerror >> bias >> x1 >> x2 >> y1 >> y2;
		DATASEC_t datasec = {x1, x2, y1, y2};
		gbn->setGain(amp, gain.d);
		gbn->setGainError(amp, gainerror.d);
		gbn->setNoise(amp, noise.d);
		gbn->setBias(amp, bias.d);
		gbn->setDatasec(amp, datasec);
	}
}

string operaIOFormats::getOrderSpacingHeader() {
	ostringstream ss;
	ss << "# Order Spacing Polynomial, the format is:\n"
	<< "# <number of coefficients> <polynomial coefficient> <polynomial coefficienterror> ... <newline>\n";
	return ss.str();
}

string operaIOFormats::getLinesFromOrderSpacing(const operaSpectralOrderVector& orders) {
	ostringstream ss;
	const Polynomial *polynomial = orders.getOrderSpacingPolynomial();
	unsigned npar = polynomial->getOrderOfPolynomial();
	ss << npar << ' ';
	for (unsigned i = 0; i < npar; i++) {
		ss << polynomial->getCoefficient(i) << ' ' << polynomial->getCoefficientError(i) << ' ';
	}
	ss << endl;
	return ss.str();
}

/* 
 * void setOrderSpacingPolynomialFromLine(string dataline)
 * \brief Sets the order spacing polynomial from a line read in from an order spacing file.
 */
void operaIOFormats::setOrderSpacingFromLine(operaSpectralOrderVector& orders, string dataline) {
	istringstream ss(dataline);
	unsigned npar= 0;
	ss >> npar;
	if(npar < 2 || npar > 6) {
		throw operaException("operaSpectralOrder: polynomial order must be between 2 and 6", operaErrorInvalidParameter, __FILE__, __FUNCTION__, __LINE__);
	}
	Polynomial *p = orders.getOrderSpacingPolynomial();
	p->setOrderOfPolynomial(npar);
	p->setChisqr(0.0);
	for(unsigned i = 0; i < npar; i++) {
		double orderSpacingPolynomialCoefficient = 0.0, orderSpacingPolynomialError = 0.0;
		ss >> orderSpacingPolynomialCoefficient >> orderSpacingPolynomialError;
		p->setCoefficient(i, orderSpacingPolynomialCoefficient);
		p->setCoefficientError(i, orderSpacingPolynomialError);
	}
}

string operaIOFormats::getDispHeader() {
	ostringstream ss;
	ss << "# Dispersion Polynomial, the format is:\n"
	<< "# <numberOfDispersionPolynomials> <newline>\n"
	<< "# <PolynomialIndex> <MinorderOfPolynomial> <MaxorderOfPolynomial> <polynomial coefficient> <polynomial coefficienterror> ... <newline>\n";
	return ss.str();
}

string operaIOFormats::getLinesFromDisp(const operaSpectralOrderVector& orders) {
	ostringstream ss;
	ss << orders.getnumberOfDispersionPolynomials() << endl;
	for (unsigned dispIndex=0; dispIndex<orders.getnumberOfDispersionPolynomials(); dispIndex++) {
		const LaurentPolynomial *polynomial = orders.getDispersionPolynomial(dispIndex);
		ss << dispIndex << ' ' << polynomial->getMinorderOfLaurentPolynomial() << ' ' << polynomial->getMaxorderOfLaurentPolynomial() << ' ';
		unsigned npar = polynomial->getNumberOfCoefficients();
		for (unsigned i = 0; i < npar; i++) {
			ss << polynomial->getCoefficient(i) << ' ' << polynomial->getCoefficientError(i) << ' ';
		}
		ss << endl;
	}
	return ss.str();
}

/* 
 * void setDispersionPolynomialFromLine(string dataline, unsigned linenumber)
 * \brief Sets a dispersion polynomial from a line read in from a dispersion file.
 * \brief First line contains total number of dispersion polynomials.
 */
void operaIOFormats::setDispFromLine(operaSpectralOrderVector& orders, string dataline, unsigned linenumber) {
	istringstream ss(dataline);
	if (linenumber == 0) {
		unsigned NumberOfDispersionPolynomials = 0;
		ss >> NumberOfDispersionPolynomials;
		orders.setnumberOfDispersionPolynomials(NumberOfDispersionPolynomials);
	} else {
		unsigned dispIndex = 0;
		int MinorderOfLaurentPolynomial = 0, MaxorderOfLaurentPolynomial = 0;
		ss >> dispIndex >> MinorderOfLaurentPolynomial >> MaxorderOfLaurentPolynomial;
		if(MaxorderOfLaurentPolynomial-MinorderOfLaurentPolynomial < 0 || MaxorderOfLaurentPolynomial-MinorderOfLaurentPolynomial > 5) {
			throw operaException("operaSpectralOrder: polynomial order must be between 2 and 6", operaErrorInvalidParameter, __FILE__, __FUNCTION__, __LINE__); //actually, checks are 1-6...
		}
		LaurentPolynomial *p = orders.getDispersionPolynomial(dispIndex);
		p->setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);
		p->setChisqr(0.0);
		for(unsigned i = 0; i < p->getNumberOfCoefficients(); i++) {
			double dispersionPolynomialCoefficient = 0.0, dispersionPolynomialError = 0.0;
			ss >> dispersionPolynomialCoefficient >> dispersionPolynomialError;
			p->setCoefficient(i, dispersionPolynomialCoefficient);
			p->setCoefficientError(i, dispersionPolynomialError);
		}
	}
}

string operaIOFormats::getRvelHeader() {
	ostringstream ss;
	ss << "# Radial Velocity Correction (km/s) format is:\n"
	<< "# <radialvelocity> <newline>\n";
	return ss.str();
}

string operaIOFormats::getLinesFromRvel(const operaSpectralOrderVector& orders) {
	ostringstream ss;
	ss << orders.getRadialVelocityCorrection() << endl;
	return ss.str();
}

/* 
 * void setRadialVelocityCorrectionFromLine(string dataline)
 * \brief Sets the radial velocity correction from a line read in from an rvel file.
 */
void operaIOFormats::setRvelFromLine(operaSpectralOrderVector& orders, string dataline) {
	istringstream ss(dataline);
	double rvel = 0.0;
	ss >> rvel;
	orders.setRadialVelocityCorrection(rvel);
}

string operaIOFormats::getTellHeader() {
	ostringstream ss;
	ss << "# Telluric Radial Velocity Correction (km/s) format is:\n"
	<< "# <radialvelocity> <newline>\n";
	return ss.str();
}

string operaIOFormats::getLinesFromTell(const operaSpectralOrderVector& orders) {
	ostringstream ss;
	ss << orders.getTelluricRadialVelocityCorrection() << endl;
	return ss.str();
}

/* 
 * void setTelluricRVCorrectionFromLine(string dataline)
 * \brief Sets the telluric radial velocity correction from a line read in from a tell file.
 */
void operaIOFormats::setTellFromLine(operaSpectralOrderVector& orders, string dataline) {
	istringstream ss(dataline);
	double tell = 0.0;
	ss >> tell;
	orders.setTelluricRadialVelocityCorrection(tell);
}

string operaIOFormats::getApertureHeader() {
	ostringstream ss;
	ss << "# Extraction Aperture format is:\n"
	<< "#\n"
	<< "# <number of orders><newline>\n"
	<< "# <order number> <number of beams> <measured tilt> <tilt error>\n"
	<< "# <leftBackgroundIndex> <xsampling> <ysampling> <lb height> <lb width> <lb slope> <lb xcenter> <lb ycenter> <lb fluxFraction>\n"
	<< "# <rightBackgroundIndex> <xsampling> <ysampling>  <rb height> <rb width> <rb slope> <rb xcenter> <rb ycenter> <rb fluxFraction>\n"
	<< "# <beam> <xsampling> <ysampling>  <beam height> <beam width> <beam slope> <beam xcenter> <beam ycenter>  <beam fluxFraction> <newline>\n"
	<< "#\n"
	<< "# Note that leftBackgroundIndex = 0\n"
	<< "# Note that rightBackgroundIndex = 1\n"
	<< "# Note that <beam> is zero-based.\n";
	return ss.str();
}

string operaIOFormats::getLineFromAperture(const operaSpectralOrder *spectralOrder) {
	ostringstream ss;
	ss << spectralOrder->getorder() << ' ' << spectralOrder->getnumberOfBeams() << ' ' << spectralOrder->getTiltInDegreesValue() << ' ' << spectralOrder->getTiltInDegreesError() << ' '; 
	for (unsigned i = 0; i < LEFTANDRIGHT+spectralOrder->getnumberOfBeams(); i++) { // loop through LEFTANDRIGHT backgrounds, then beams
		const operaExtractionAperture *aperture;
		unsigned b = i;
		if(i < LEFTANDRIGHT) aperture = spectralOrder->getBackgroundApertures(b); //first two (LEFTANDRIGHT) iterations are background
		else {
			b = i-LEFTANDRIGHT; //rest are the beams
			aperture = spectralOrder->getExtractionApertures(b);
		}
		const Line *lineAperture = aperture->getLineAperture();
		ss << b << ' ' << aperture->getXsampling() << ' ' << aperture->getYsampling() << ' ' << lineAperture->getWidth() << ' ' << lineAperture->getLength() << ' ' << lineAperture->getSlope() << ' '
		<< lineAperture->getMidPoint()->getXcoord() << ' ' << lineAperture->getMidPoint()->getYcoord() << ' ' << aperture->getFluxFraction() << ' ';
	}
	return ss.str();
}

/* 
 * void setApertureFromLine(operaSpectralOrder *spectralOrder, string dataline)
 * \brief Takes a line from an aperture file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setApertureFromLine(operaSpectralOrder *spectralOrder, string dataline) {
	istringstream ss(dataline);
	unsigned order, beams;
	Double tiltInDegreesValue, tiltInDegreesError;
	ss >> order >> beams >> tiltInDegreesValue >> tiltInDegreesError;
	spectralOrder->sethasExtractionApertures(true);
	spectralOrder->setnumberOfBeams(beams);
	spectralOrder->setTiltInDegrees(tiltInDegreesValue.d, tiltInDegreesError.d);
	for (unsigned i = 0; i < LEFTANDRIGHT+beams; i++) { // loop through LEFTANDRIGHT backgrounds, then beams
		unsigned b, xsampling, ysampling;
		Double width, length, slope, midpointx, midpointy, fluxfraction;
		ss >> b >> xsampling >> ysampling >> width >> length >> slope >> midpointx >> midpointy >> fluxfraction;
		operaPoint point(midpointx.d, midpointy.d);
		Line LineAperture(slope.d, width.d, length.d, &point);
		operaExtractionAperture *Aperture = new operaExtractionAperture(&LineAperture, xsampling, ysampling);
		Aperture->setFluxFraction(fluxfraction.d);
		if(i < LEFTANDRIGHT) spectralOrder->setBackgroundApertures(b, Aperture); //first two (LEFTANDRIGHT) iterations are background
		else spectralOrder->setExtractionApertures(b, Aperture); //rest are the beams
	}
	spectralOrder->sethasExtractionApertures(true);
}

string operaIOFormats::getFcalHeader() {
	ostringstream ss;
	ss << "# Flux Calibration Beam Spectrum format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <nElements> <nBeams> <wavelengthForNormalization> <elementindex> <wavelength> <SpectralElements flux conversion> <flux conversion variance> <SpectralElements throughput> <throughput variance> \n"
	<< "# <beam> <BeamSED[beam] flux conversion> <flux conversion variance> <BeamSED[beam] throughput> <throughput variance> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromFcal(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
	const operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();
	const operaSpectralElements *InstrumentThroughput = spectralEnergyDistribution->getThroughputElements();
	ss << spectralOrder->getorder() << ' ' << FluxCalibration->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' '
	<< fixed << setprecision(4) << spectralEnergyDistribution->getwavelengthForNormalization() << ' ' << index << ' ' << FluxCalibration->getwavelength(index) << ' '
	<< scientific << setprecision(6) << FluxCalibration->getFlux(index) << ' ' << FluxCalibration->getFluxVariance(index) << ' '
	<< InstrumentThroughput->getFlux(index) << ' ' << InstrumentThroughput->getFluxVariance(index) << ' ';
	for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
		ss << beam << ' ' << spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFlux(index) << ' ' << spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFluxVariance(index) << ' '
		<< spectralOrder->getBeamSED(beam)->getThroughputElements()->getFlux(index) << ' ' << spectralOrder->getBeamSED(beam)->getThroughputElements()->getFluxVariance(index) << ' ';
	}
	return ss.str();
}

/* 
 * void setFcalFromLine(operaSpectralOrder *spectralOrder, string dataline)
 * \brief Takes a line from a flux calibration .fcal file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setFcalFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder) {
	istringstream ss(dataline);
	unsigned order, nElements, beams, elementindex;
	double wavelengthForNormalization;
	ss >> order >> nElements >> beams >> wavelengthForNormalization >> elementindex;
#ifdef RANGE_CHECK
	if (elementindex > nElements) {
		throw operaException("operaIOFormats: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (beams != spectralOrder->getnumberOfBeams()) beams = 0; //Ignore all beams when the number of beams are different (previously done in operaSpectralOrderVector::correctFlatField)
	if (neworder) {
		spectralOrder->setnumberOfBeams(beams);
		spectralOrder->createSpectralEnergyDistributionElements(nElements);
		operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
		spectralEnergyDistribution->getFluxCalibrationElements()->setHasWavelength(true); //we don't have checks for this when writing... should we? should this even be here?
		spectralEnergyDistribution->getThroughputElements()->setHasWavelength(true); //we don't have checks for this when writing... should we? should this even be here?
		spectralEnergyDistribution->setHasFluxCalibration(true);
		spectralEnergyDistribution->setHasInstrumentThroughput(true);
		spectralOrder->sethasSpectralEnergyDistribution(true);
	}
	operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
	spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
	Double wl, fluxcal, fcalvariance, throughput, throughputvariance;
	ss >> wl >> fluxcal >> fcalvariance >> throughput >> throughputvariance;
	operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();
	FluxCalibration->setwavelength(wl.d, elementindex);
	FluxCalibration->setFlux(fluxcal.d, elementindex);
	FluxCalibration->setFluxVariance(fcalvariance.d, elementindex);
	FluxCalibration->setphotoCenter(0.0, 0.0, elementindex);
	operaSpectralElements *InstrumentThroughput = spectralEnergyDistribution->getThroughputElements();
	InstrumentThroughput->setwavelength(wl.d, elementindex);
	InstrumentThroughput->setFlux(throughput.d, elementindex);
	InstrumentThroughput->setFluxVariance(throughputvariance.d, elementindex);
	InstrumentThroughput->setphotoCenter(0.0, 0.0, elementindex);
	for (unsigned b=0; b < beams; b++) {
		unsigned beam = 0;
		Double beamfluxcal, beamfcalvariance, beamthroughput, beamthrouputvariance;
		ss >> beam >> beamfluxcal >> beamfcalvariance >> beamthroughput >> beamthrouputvariance;
		operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);  
		operaSpectralElements *beamFluxcalibration = BeamSED->getFluxCalibrationElements();
		beamFluxcalibration->setwavelength(wl.d, elementindex);
		beamFluxcalibration->setFlux(beamfluxcal.d, elementindex);
		beamFluxcalibration->setFluxVariance(beamfcalvariance.d, elementindex);
		beamFluxcalibration->setphotoCenter(0.0, 0.0, elementindex);
		operaSpectralElements *beamThroughput = BeamSED->getThroughputElements();
		beamThroughput->setwavelength(wl.d, elementindex);
		beamThroughput->setFlux(beamthroughput.d, elementindex);
		beamThroughput->setFluxVariance(beamthrouputvariance.d, elementindex);
		beamThroughput->setphotoCenter(0.0, 0.0, elementindex);
	}
}

string operaIOFormats::getGeomHeader() {
	ostringstream ss;
	ss << "# Geometry Calibration, the format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <number of coefficients> <ndatapoints> <polynomial coefficient> <polynomial coefficienterror>... <chisqr> <YBinning> <miny><maxy> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromGeom(const operaSpectralOrder *spectralOrder) {
	ostringstream ss;
	const Polynomial *polynomial = spectralOrder->getGeometry()->getCenterPolynomial();
	ss << spectralOrder->getorder() << ' ' << polynomial->getOrderOfPolynomial() << ' ' << spectralOrder->getGeometry()->getNdatapoints() << ' ';
	for (unsigned coeff=0; coeff<polynomial->getOrderOfPolynomial(); coeff++) {
		ss << polynomial->getCoefficient(coeff) << ' ' << polynomial->getCoefficientError(coeff) << ' ';
	}
	ss << polynomial->getChisqr() << ' ' << spectralOrder->getGeometry()->getNumberofPointsToBinInYDirection() << ' ' << spectralOrder->getGeometry()->getYmin() << ' ' << spectralOrder->getGeometry()->getYmax();
	return ss.str();
}

/* 
 * void setGeometryFromLine(operaSpectralOrder *spectralOrder, string dataline)
 * \brief Takes a line of geometry polynomial coefficients from a geometry .geom file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setGeomFromLine(operaSpectralOrder *spectralOrder, string dataline) {
	istringstream ss(dataline);
	unsigned order, npar;
	int ndatapoints;
	ss >> order >> npar >> ndatapoints;
	spectralOrder->createGeometry(ndatapoints, ndatapoints);
	spectralOrder->sethasGeometry(true);
	operaGeometry *geometry = spectralOrder->getGeometry();
	Polynomial *p = geometry->getCenterPolynomial();
	p->setOrderOfPolynomial(npar);
	for (unsigned i=0; i<npar; i++) {
		Float coeff, coefferr;
		ss >> coeff >> coefferr;
		p->setCoefficient(i, coeff.f);
		p->setCoefficientError(i, coefferr.f);
	}
	Float chisqr, miny, maxy;
	unsigned ybinning;
	ss >> chisqr >> ybinning >> miny >> maxy;
	p->setChisqr(chisqr.f);
	geometry->setNumberofPointsToBinInYDirection(ybinning);
	geometry->setYmin(miny.f);
	geometry->setYmax(maxy.f);
	geometry->CalculateAndSetOrderLength();
}

string operaIOFormats::getWaveHeader() {
	ostringstream ss;
	ss << "# Wavelength Calibration, the format is:\n"
	<< "# <count of orders> <newline>\n"
	<< "# <order number> <number of coefficients> <polynomial coefficients> <polynomial coefficient error> ... <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromWave(const operaSpectralOrder *spectralOrder) {
	ostringstream ss;
	const Polynomial *polynomial = spectralOrder->getWavelength()->getWavelengthPolynomial();
	ss << spectralOrder->getorder() << ' ' << polynomial->getOrderOfPolynomial() << ' ';
	for	(unsigned coeff=0; coeff<polynomial->getOrderOfPolynomial(); coeff++) {
		ss << scientific << polynomial->getCoefficient(coeff) << ' ' << polynomial->getCoefficientError(coeff) << ' ';
	}
	return ss.str();
}

/* 
 * void setWavelengthFromLine(operaSpectralOrder *spectralOrder, string dataline)
 * \brief Takes a line of wavelength polynomial coefficients and inserts the values into the order specified by that line.
 */
void operaIOFormats::setWaveFromLine(operaSpectralOrder *spectralOrder, string dataline) {
	istringstream ss(dataline);
	unsigned order, npar;
	ss >> order >> npar;
	if(!spectralOrder->getWavelength()) spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL); //Why does a wavelgnth already exist here? Don't know, but it crashes if we don't check...
	spectralOrder->sethasWavelength(true);
	Polynomial *poly = spectralOrder->getWavelength()->getWavelengthPolynomial();
	poly->setOrderOfPolynomial(npar);
	for (unsigned i=0; i < npar; i++) {
		Float coeff, coefferr;
		ss >> coeff >> coefferr;
		poly->setCoefficient(i, coeff.f);
		poly->setCoefficientError(i, coefferr.f);
	}
}

string operaIOFormats::getProfHeader() {
	ostringstream ss;
	ss << "# Instrument Profile Calibration the format is:\n"
	<< "# <number of orders> <number of columns i> <number of rows j> <xsize> <xsampling> <ysize> <ysampling> <newline>\n"
	<< "# <order number> <i> <j> <number of coefficients> <ndatapoints> <polynomial coefficients> <chisqr> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromProf(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
	// Note: Instead of passing i and j, we pass a single index to emulate a nested 2D loop: j=[0,y) { i=[0,x) }
	// Iteration n of our actual loop will emulating the iteration j of our outer loop and iteration i of our inner loop
	// In other words, n = x*j + i, which means...
	unsigned j = index / instrumentProfile->getNXPoints(); // j = n/x (integer division means i/x = 0, since i < x)
	unsigned i = index % instrumentProfile->getNXPoints(); // i = n mod x (since x*j mod x = 0 and i mod x = i)
	const PolynomialCoeffs_t *pp = instrumentProfile->getipPolyModelCoefficients(i,j);
	ss << spectralOrder->getorder() << ' ' << i << ' ' << j << ' ' << pp->orderofPolynomial << ' ' << '0'/*instrumentProfile->getnDataPoints()*/ << ' ';
	for	(unsigned coeff=0; coeff<pp->orderofPolynomial; coeff++) {
		ss << pp->p[coeff] << ' ';
	}
	ss << instrumentProfile->getchisqrMatrixValue(i, j);
	return ss.str();
}

/* 
 * void setInstrumentProfileFromLine(operaSpectralOrder *spectralOrder, string dataline, bool existingorder, unsigned xsize, unsigned xsampling, unsigned ysize, unsigned ysampling)
 * \brief Takes a line from instrument profile and inserts the values into the order specified by that line.
 */
void operaIOFormats::setProfFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned xsize, unsigned xsampling, unsigned ysize, unsigned ysampling) {
	istringstream ss(dataline);
	unsigned order, npar, col, row, ndatapoints;
	ss >> order >> col >> row >> npar >> ndatapoints;
	if(npar < 1 || npar > 6) throw operaException("operaSpectralOrder: profile polynomial order must be between 1 and 6, got "+itos(npar), operaErrorInvalidParameter, __FILE__, __FUNCTION__, __LINE__);
	
	if (neworder) {
		spectralOrder->setInstrumentProfileVector(xsize, xsampling, ysize, ysampling, 1);
		spectralOrder->sethasInstrumentProfile(true);
	}
	PolynomialCoeffs_t pp;
	pp.orderofPolynomial = npar;
	for (unsigned i=0; i < npar; i++) {
		Float ProfileCoeff;
		ss >> ProfileCoeff;
		pp.p[i] = ProfileCoeff.f;
	}
	Float chisqr;
	ss >> chisqr;
	operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
	instrumentProfile->setipPolyModelCoefficients(pp, col, row);
	instrumentProfile->setchisqrMatrixValue(chisqr.f, col, row);
}

string operaIOFormats::getSNRHeader() {
	ostringstream ss;
	ss << "# SNR format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <center wavelength> <center SNR per spectral bin> <center SNR per CCD pixel> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromSNR(const operaSpectralOrder *spectralOrder) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << spectralOrder->getorder() << ' ' << fixed << setprecision(4) << spectralElements->getwavelength(spectralElements->getnSpectralElements() / 2) << ' '
	<< scientific << spectralOrder->getCenterSNR() << ' ' << spectralOrder->getCenterSNR()/sqrt(spectralOrder->getsnrSpectralBinSize());
	return ss.str();
}

/* 
 * void setSNRFromLine(operaSpectralOrder *spectralOrder, string dataline)
 * \brief Takes a line from an SNR table and inserts the values into the order specified by that line.
 */
void operaIOFormats::setSNRFromLine(operaSpectralOrder *spectralOrder, string dataline) {
	istringstream ss(dataline);
	unsigned order;
	Float wl = 0.0, centersnr = 0.0, snrperpix = 0.0;
	ss >> order >> wl >> centersnr >> snrperpix;
	spectralOrder->sethasSNR(true);
	spectralOrder->sethasCenterSNROnly(true);
	spectralOrder->setCenterSNR(centersnr.f);
	spectralOrder->setsnrSpectralBinSize((centersnr.f*centersnr.f)/(snrperpix.f*snrperpix.f));
}

string operaIOFormats::getSpectrumHeader() {
	ostringstream ss;
	ss << "# Spectrum format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <distance> <flux> <variance> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromSpectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << spectralOrder->getorder() << ' ' << spectralElements->getdistd(index) << ' ' << spectralElements->getFlux(index) << ' ' << spectralElements->getFluxVariance(index);
	return ss.str();
}

/* 
 * void setSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, unsigned index, operaSpectralOrder_t format)
 * \brief Takes a line from a spectrum .s file and inserts the values into the order specified by that line.
 * \brief The raw spectrum has distances rather than wavelength.
 */
void operaIOFormats::setSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, operaSpectralOrder_t format) {
	istringstream ss(dataline);
	unsigned order;
	ss >> order;
	if (neworder) {
		spectralOrder->createSpectralElements(MAXSPECTRALELEMENTSPERORDER, format);
		spectralOrder->sethasSpectralElements(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasDistance(true);
	}
	Double dist, flux, fluxvariance;
	ss >> dist >> flux >> fluxvariance;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setnSpectralElements(index+1);
	spectralElements->setdistd(dist.d, index);
	spectralElements->setFlux(flux.d, index);
	spectralElements->setFluxVariance(fluxvariance.d, index);
}

string operaIOFormats::getCalibratedSpectrumHeader() {
	ostringstream ss;
	ss << "# Calibrated Spectrum format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <wavelength> <flux> <variance> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromCalibratedSpectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << spectralOrder->getorder() << ' ' << fixed << setprecision(4) << spectralElements->getwavelength(index) << ' '
	<< scientific << spectralElements->getFlux(index) << ' ' << spectralElements->getFluxVariance(index);
	return ss.str();
}

/* 
 * void setCalibratedSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, unsigned index, operaSpectralOrder_t format)
 * \brief Takes a line from a calibrated spectrum .s file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setCalibratedSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, operaSpectralOrder_t format) {
	istringstream ss(dataline);
	unsigned order;
	ss >> order;
	if (neworder) {
		spectralOrder->createSpectralElements(MAXSPECTRALELEMENTSPERORDER, format);
		spectralOrder->sethasSpectralElements(true);
		spectralOrder->setSpectrumType(format);
		spectralOrder->sethasWavelength(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasWavelength(true);
	}
	Double wavelength, flux, fluxvariance;
	ss >> wavelength >> flux >> fluxvariance;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setnSpectralElements(index+1);
	spectralElements->setwavelength(wavelength.d, index);
	spectralElements->setFlux(flux.d, index);
	spectralElements->setFluxVariance(fluxvariance.d, index);
}

string operaIOFormats::getBeamSpectrumHeader() {
	ostringstream ss;
	ss << "# Beam Spectrum format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <nElements> <nBeams> <elementindex> <SpectralElements photoCenterX> <SpectralElements photoCenterY> <SpectralElements dist> <SpectralElements flux> <SpectralElements flux variance> <XCorrelation>\n"
	<< "# <beam> <BeamElements[beam] photoCenterX> <BeamElements[beam] photoCenterY> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromBeamSpectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ' << index << ' '
	<< spectralElements->getphotoCenterX(index) << ' '  << spectralElements->getphotoCenterY(index) << ' '
	<< spectralElements->getdistd(index) << ' '
	<< spectralElements->getFlux(index) << ' ' << spectralElements->getFluxVariance(index) << ' '
	<< spectralElements->getXCorrelation(index) << ' ';
	for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
		const operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
		ss << beam << ' '
		<< beamElements->getphotoCenterX(index) << ' ' << beamElements->getphotoCenterY(index) << ' '
		<< beamElements->getFlux(index) << ' ' << beamElements->getFluxVariance(index) << ' ';
	}
	return ss.str();
}

/* 
 * void setBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, operaSpectralOrder_t format)
 * \brief Takes a line from a beam spectrum .e file and inserts the values into the order specified by that line.
 * \brief The raw spectrum has distances rather than wavelength.
 */
void operaIOFormats::setBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, operaSpectralOrder_t format) {
	istringstream ss(dataline);
	unsigned order, nElements, beams, elementindex;
	ss >> order >> nElements >> beams >> elementindex;
#ifdef RANGE_CHECK
	if (elementindex > nElements) {
		throw operaException("operaIOFormats: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (neworder) {
		spectralOrder->createSpectralElements(nElements, format);
		spectralOrder->createBeamsAndBackgrounds(nElements, beams, format);
		spectralOrder->setnumberOfBeams(beams);
		spectralOrder->setSpectrumType(format);
		spectralOrder->sethasSpectralElements(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasDistance(true);
		spectralElements->setHasXCorrelation(true);
	}
	Double photoCenterX, photoCenterY, dist, flux, variance, xcorrelation;
	ss >> photoCenterX >> photoCenterY >> dist >> flux >> variance >> xcorrelation;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setdistd(dist.d, elementindex);
	spectralElements->setFlux(flux.d, elementindex);
	spectralElements->setFluxVariance(variance.d, elementindex);
	spectralElements->setXCorrelation(xcorrelation.d, elementindex);
	spectralElements->setphotoCenter(photoCenterX.d, photoCenterY.d, elementindex);
	for (unsigned b=0; b < beams; b++) {
		unsigned beam = 0;
		Double beamphotoCenterX, beamphotoCenterY, beamflux, beamvariance;
		ss >> beam >> beamphotoCenterX >> beamphotoCenterY >> beamflux >> beamvariance;
		operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
		beamElements->setFlux(beamflux.d, elementindex);
		beamElements->setFluxVariance(beamvariance.d, elementindex);
		beamElements->setphotoCenter(beamphotoCenterX.d, beamphotoCenterY.d, elementindex);
	}
}

string operaIOFormats::getCalibratedBeamSpectrumHeader() {
	ostringstream ss;
	ss << "# Calibrated Beam Spectrum format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <SpectralElements flux> <SpectralElements flux variance><cross correlation> \n"
	<< "# <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromCalibratedBeamSpectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ' << index << ' '
	<< fixed << setprecision(4) << spectralElements->getwavelength(index) << ' '
	<< scientific << spectralElements->getFlux(index) << ' ' << spectralElements->getFluxVariance(index) << ' ' << spectralElements->getXCorrelation(index) << ' ';
	for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
		const operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
		ss << beam << ' ' << beamElements->getFlux(index) << ' ' << beamElements->getFluxVariance(index) << ' ';             
	}
	return ss.str();
}

/* 
 * void setCalibratedBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, operaSpectralOrder_t format)
 * \brief Takes a line from a calibrated beam spectrum .e file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setCalibratedBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, operaSpectralOrder_t format) {
	istringstream ss(dataline);
	unsigned order, nElements, beams, elementindex;
	ss >> order >> nElements >> beams >> elementindex;
#ifdef RANGE_CHECK
	if (elementindex > nElements) {
		throw operaException("operaIOFormats: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (neworder) {
		spectralOrder->createSpectralElements(nElements, format);
		spectralOrder->createBeamsAndBackgrounds(nElements, beams, format);
		spectralOrder->setnumberOfBeams(beams);
		spectralOrder->setSpectrumType(format);
		spectralOrder->sethasSpectralElements(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasWavelength(true);
		spectralElements->setHasXCorrelation(true);
	}
	Double wl, flux, variance, xcorrelation;
	ss >> wl >> flux >> variance >> xcorrelation;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setwavelength(wl.d, elementindex);
	spectralElements->setFlux(flux.d, elementindex);
	spectralElements->setFluxVariance(variance.d, elementindex);
	spectralElements->setXCorrelation(xcorrelation.d, elementindex);
	for (unsigned b=0; b < beams; b++) {
		unsigned beam = 0;
		Double beamflux, beamvariance;
		ss >> beam >> beamflux >> beamvariance;
		operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
		beamElements->setFlux(beamflux.d, elementindex);
		beamElements->setFluxVariance(beamvariance.d, elementindex);
	}
}

string operaIOFormats::getCalibratedExtendedBeamSpectrumHeader() {
	ostringstream ss;
	ss << "# Calibrated Extended Beam Spectrum (as output from operaExtraction) format is:\n"
	<< "# <number of orders> <newline>\n"
	<< "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <wavelength telluric corrected> <heliocentric wavelength correction> <crosscorrelation>\n"
	<< "# <rawFlux> <rawFlux variance> <normalizedFlux> <normalizedFlux variance> <fcalFlux> <fcalFlux variance>\n"
	<< "# <beamindex> <beam rawFlux> <beam rawFlux variance> <beam normalizedFlux> <beam normalizedFlux variance> <beam fcalFlux> <beam fcalFlux variance> ... <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromCalibratedExtendedBeamSpectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	double relativeVariance = spectralElements->getrawFluxVariance(index) / (spectralElements->getrawFlux(index) * spectralElements->getrawFlux(index));
	ss << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ' << index << ' '
	<< fixed << setprecision(8) << spectralElements->getwavelength(index) << ' '
	<< spectralElements->gettell(index) << ' ' << spectralElements->getrvel(index) << ' '
	<< scientific << spectralElements->getXCorrelation(index) << ' ' << spectralElements->getrawFlux(index) << ' ' << spectralElements->getrawFluxVariance(index) << ' '
	<< spectralElements->getnormalizedFlux(index) << ' ' << relativeVariance * spectralElements->getnormalizedFlux(index) * spectralElements->getnormalizedFlux(index) << ' '
	<< spectralElements->getfcalFlux(index) << ' ' << relativeVariance * spectralElements->getfcalFlux(index) * spectralElements->getfcalFlux(index) << ' ';
	
	if (spectralOrder->getnumberOfBeams() > 1) {
		for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
			const operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
			double relativeBeamVariance = beamElements->getrawFluxVariance(index) / (beamElements->getrawFlux(index) * beamElements->getrawFlux(index));
			ss << beam << ' ' << beamElements->getrawFlux(index) << ' ' << beamElements->getrawFluxVariance(index) << ' '
			<< beamElements->getnormalizedFlux(index) << ' ' << relativeBeamVariance * beamElements->getnormalizedFlux(index) * beamElements->getnormalizedFlux(index) << ' '
			<< beamElements->getfcalFlux(index) << ' ' << relativeBeamVariance * beamElements->getfcalFlux(index) * beamElements->getfcalFlux(index) << ' ';
		}
	}
	return ss.str();
}

/* 
 * void setCalibratedExtendedBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, operaSpectralOrder_t format)
 * \brief Takes a line from a calibrated extended beam spectrum .spc file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setCalibratedExtendedBeamSpectrumFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, operaSpectralOrder_t format) {
	istringstream ss(dataline);
	unsigned order, nElements, beams, elementindex;
	ss >> order >> nElements >> beams >> elementindex;
#ifdef RANGE_CHECK
	if (elementindex > nElements) {
		throw operaException("operaIOFormats: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (neworder) {
		spectralOrder->createSpectralElements(nElements, format, true);
		spectralOrder->createBeamsAndBackgrounds(nElements, beams, format, true);
		spectralOrder->setnumberOfBeams(beams);
		spectralOrder->setSpectrumType(format);
		spectralOrder->sethasSpectralElements(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasWavelength(true);
		spectralElements->setHasXCorrelation(true);
		spectralElements->setHasExtendedBeamFlux(true);
	}
	Double wl, tell, rvel, xcorrelation, rawFlux, rawFluxError, normalizedFlux, normalizedFluxError, fcalFlux, fcalFluxError;
	ss >> wl >> tell >> rvel >> xcorrelation >> rawFlux >> rawFluxError >> normalizedFlux >> normalizedFluxError >> fcalFlux >> fcalFluxError;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setwavelength(wl.d, elementindex);
	spectralElements->settell(tell.d, elementindex);
	spectralElements->setrvel(rvel.d, elementindex);
	spectralElements->setXCorrelation(xcorrelation.d, elementindex);
	spectralElements->setFlux(rawFlux.d, elementindex);
	spectralElements->setFluxVariance(rawFluxError.d, elementindex);
	spectralElements->setrawFlux(rawFlux.d,elementindex);
	spectralElements->setrawFluxVariance(rawFluxError.d, elementindex);
	spectralElements->setnormalizedFlux(normalizedFlux.d, elementindex);
	spectralElements->setnormalizedFluxVariance(normalizedFluxError.d, elementindex);
	spectralElements->setfcalFlux(fcalFlux.d, elementindex);
	spectralElements->setfcalFluxVariance(fcalFluxError.d, elementindex);
	for (unsigned b=0; b < beams; b++) {
		unsigned beam = b;
		if (beams > 1) ss >> beam >> rawFlux >> rawFluxError >> normalizedFlux >> normalizedFluxError >> fcalFlux >> fcalFluxError;
		operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
		beamElements->setFlux(rawFlux.d, elementindex);
		beamElements->setFluxVariance(rawFluxError.d, elementindex);
		beamElements->setrawFlux(rawFlux.d,elementindex);
		beamElements->setrawFluxVariance(rawFluxError.d, elementindex);
		beamElements->setnormalizedFlux(normalizedFlux.d, elementindex);
		beamElements->setnormalizedFluxVariance(normalizedFluxError.d, elementindex);
		beamElements->setfcalFlux(fcalFlux.d, elementindex);
		beamElements->setfcalFluxVariance(fcalFluxError.d, elementindex);
	}
}

string operaIOFormats::getPolarimetryHeader() {
	ostringstream ss;
	ss << "# Polarimetry format is:\n"
	<< "# <number of orders> <cols> <method> <newline>\n"
	<< "# <order number> <StokesParameter_t> <length> <distance> <wavelength> <crosscorrelation> <Stokes(Q,U,V) flux> <Stokes(Q,U,V) variance> <StokesI flux> <StokesI variance> <degree of polarization flux> <degree of polarization variance> <first null polarization> <first null polarization variance> <second null polarization> <second null polarization variance> <newline>\n"
	<< "# ...\n";
	return ss.str();
}

string operaIOFormats::getLineFromPolarimetry(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	stokes_parameter_t StokesParameter = getStokesParameter(spectralOrder); //would be faster to only call this once per order, difference probably negligable vs I/O though?
	const operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
	const operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
	const operaStokesVector *StokesVector = Polarimetry->getStokesVector();
	const operaStokesVector *DegreeOfPolarization = Polarimetry->getDegreeOfPolarization();
	ss << spectralOrder->getorder() << ' ' << StokesParameter << ' ' << Polarimetry->getLength() << ' '
	<< scientific << setprecision(5) << SpectralElements->getdistd(index) << ' '
	<< fixed << setprecision(4) << SpectralElements->getwavelength(index) << ' '
	<< scientific << SpectralElements->getXCorrelation(index) << ' '
	<< StokesVector->getStokesParameterFlux(StokesParameter, index) << ' ' << StokesVector->getStokesParameterVariance(StokesParameter, index) << ' '
	<< StokesVector->getStokesParameterFlux(StokesI, index) << ' ' << StokesVector->getStokesParameterVariance(StokesI, index) << ' '
	<< DegreeOfPolarization->getStokesParameterFlux(StokesParameter, index) << ' ' << DegreeOfPolarization->getStokesParameterVariance(StokesParameter, index) << ' ';
	if (Polarimetry->getHasFirstNullPolarization()) {
		const operaStokesVector *FirstNullPolarization = Polarimetry->getFirstNullPolarization();
		ss << FirstNullPolarization->getStokesParameterFlux(StokesParameter, index) << ' ' << FirstNullPolarization->getStokesParameterVariance(StokesParameter, index) << ' ';
	}
	else ss << "0.0 0.0 ";
	if (Polarimetry->getHasSecondNullPolarization()) {
		const operaStokesVector *SecondNullPolarization = Polarimetry->getSecondNullPolarization();
		ss << SecondNullPolarization->getStokesParameterFlux(StokesParameter, index) << ' ' << SecondNullPolarization->getStokesParameterVariance(StokesParameter, index) << ' ';
	}
	else ss << "0.0 0.0 ";
	return ss.str();
}

/* 
 * void setPolarFromLine(operaSpectralOrder *spectralOrder, string dataline, unsigned index, unsigned method)
 * \brief Takes a line from a polar .p file and inserts the values into the order specified by that line.
 * \brief The raw spectrum has distances rather than wavelength.
 */
void operaIOFormats::setPolarFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned index, unsigned method) {
	istringstream ss(dataline);
	unsigned order = 0, stokespar = 0, nElements = 0;
	ss >> order >> stokespar >> nElements;
	stokes_parameter_t StokesParameter = stokes_parameter_t(stokespar);
	if (neworder) {
		spectralOrder->createPolarimetry(nElements);
		spectralOrder->sethasPolarimetry(true);
		spectralOrder->createSpectralElements(nElements, None);
		spectralOrder->sethasSpectralElements(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasWavelength(true);
		spectralElements->setHasDistance(true);
		spectralElements->setHasXCorrelation(true);
		operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
		Polarimetry->setHasWavelength(true);
		Polarimetry->setmethod((method_t)method);
	}
	Double distance, wavelength, crosscorrelation, QUVFlux, QUVVar, IFlux, IVar, DegPolarFlux, DegPolarVar, FirstNullPolar, FirstNullPolarVar, SecondNullPolar, SecondNullPolarVar;
	ss >> distance >> wavelength >> crosscorrelation >> QUVFlux >> QUVVar >> IFlux >> IVar >> DegPolarFlux >> DegPolarVar >> FirstNullPolar >> FirstNullPolarVar >> SecondNullPolar >> SecondNullPolarVar;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setdistd(distance.d, index);
	spectralElements->setwavelength(wavelength.d, index);
	spectralElements->setXCorrelation(crosscorrelation.d, index);
	spectralElements->setFlux(IFlux.d, index);
	spectralElements->setFluxVariance(IVar.d, index);
	operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
	Polarimetry->setwavelength(wavelength.d, index);
	Polarimetry->setStokesParameter(StokesParameter, QUVFlux.d, QUVVar.d, index);
	Polarimetry->setStokesParameter(StokesI, IFlux.d, IVar.d, index);
	Polarimetry->setDegreeOfPolarization(StokesParameter, DegPolarFlux.d, DegPolarVar.d, index);
	Polarimetry->setFirstNullPolarization(StokesParameter, FirstNullPolar.d, FirstNullPolarVar.d, index);
	Polarimetry->setSecondNullPolarization(StokesParameter, SecondNullPolar.d, SecondNullPolarVar.d, index);
	Polarimetry->setHasStokesI(true);
	if (StokesParameter == StokesQ) {
		Polarimetry->setHasStokesQ(true);
		Polarimetry->setHasDegreeOfStokesQ(true);
	}
	if (StokesParameter == StokesU) {
		Polarimetry->setHasStokesU(true);
		Polarimetry->setHasDegreeOfStokesU(true);
	}
	if (StokesParameter == StokesV) {
		Polarimetry->setHasStokesV(true);
		Polarimetry->setHasDegreeOfStokesV(true);
	}
}

string operaIOFormats::getExtendedPolarimetryHeader() {
	ostringstream ss;
	ss << "# Extended Polarimetry format is:\n"
	<< "# <number of orders> <StokesParameter_t> <method> <newline>\n"
	<< "# <order number> <nElements> <elementindex> <wavelength> <wavelength telluric corrected> <heliocentric wavelength correction> <crosscorrelation>\n"
	<< "# <StokesI flux> <StokesI variance> <normalized StokesI flux> <calibrated StokesI flux> <normalized StokesI flux variance> <calibrated StokesI flux variance>\n"
	<< "# <degree of polarization> <degree of polarization variance> <first null polarization> <second null polarization> <newline>\n";
	return ss.str();
}

string operaIOFormats::getLineFromExtendedPolarimetry(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	stokes_parameter_t StokesParameter = getStokesParameter(spectralOrder); //would be faster to only call this once per order, difference probably negligable vs I/O though?
	const operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
	const operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
	const operaStokesVector *DegreeOfPolarization = Polarimetry->getDegreeOfPolarization();
	double relativeVariance = SpectralElements->getrawFluxVariance(index) / (SpectralElements->getrawFlux(index) * SpectralElements->getrawFlux(index));
	ss << spectralOrder->getorder() << ' ' << Polarimetry->getLength() << ' ' << index << ' '
	<< fixed << setprecision(8) << SpectralElements->getwavelength(index) << ' '
	<< SpectralElements->gettell(index) << ' ' << SpectralElements->getrvel(index) << ' '
	<< scientific << SpectralElements->getXCorrelation(index) << ' '
	<< SpectralElements->getrawFlux(index) << ' ' << SpectralElements->getrawFluxVariance(index) << ' '
	<< SpectralElements->getnormalizedFlux(index) << ' ' << SpectralElements->getfcalFlux(index) << ' '
	<< relativeVariance * SpectralElements->getnormalizedFlux(index) * SpectralElements->getnormalizedFlux(index) << ' '
	<< relativeVariance * SpectralElements->getfcalFlux(index) * SpectralElements->getfcalFlux(index) << ' '
	<< DegreeOfPolarization->getStokesParameterFlux(StokesParameter, index) << ' ' << DegreeOfPolarization->getStokesParameterVariance(StokesParameter, index) << ' ';
	if (Polarimetry->getHasFirstNullPolarization()) {
		const operaStokesVector *FirstNullPolarization = Polarimetry->getFirstNullPolarization();
		ss  << FirstNullPolarization->getStokesParameterFlux(StokesParameter, index) << ' ';
	}
	else ss << "0.0 ";
	if (Polarimetry->getHasSecondNullPolarization()) {
		const operaStokesVector *SecondNullPolarization = Polarimetry->getSecondNullPolarization();
		ss << SecondNullPolarization->getStokesParameterFlux(StokesParameter, index) << ' ';
	}
	else ss << "0.0 ";
	return ss.str();
}

/* 
 * void setExtendedPolarimetryFromLine(operaSpectralOrder *spectralOrder, string dataline, unsigned StokesParameter, unsigned method)
 * \brief Takes a line from a extended polarimetry .pol file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setExtendedPolarimetryFromLine(operaSpectralOrder *spectralOrder, string dataline, bool neworder, unsigned stokespar, unsigned method) {
	istringstream ss(dataline);
	unsigned order, nElements, index;
	ss >> order >> nElements >> index;
	stokes_parameter_t StokesParameter = stokes_parameter_t(stokespar);
	if (neworder) {
		spectralOrder->createPolarimetry(nElements);
		spectralOrder->sethasPolarimetry(true);
		spectralOrder->createSpectralElements(nElements, None, true);
		spectralOrder->sethasSpectralElements(true);
		operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
		spectralElements->setHasWavelength(true);
		spectralElements->setHasDistance(true);
		spectralElements->setHasXCorrelation(true);
		operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
		Polarimetry->setHasWavelength(true);
		Polarimetry->setmethod((method_t)method);
	}
	Double wavelength, tell, rvel, crosscorrelation, IFlux, IVar, normalizedFlux, fcalFlux, normalizedFluxFluxError, fcalFluxError, DegPolarFlux, DegPolarVar, FirstNullPolar, SecondNullPolar;
	ss >> wavelength >> tell >> rvel >> crosscorrelation >> IFlux >> IVar >> normalizedFlux >> fcalFlux >> normalizedFluxFluxError >> fcalFluxError >> DegPolarFlux >> DegPolarVar >> FirstNullPolar >> SecondNullPolar;
	operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	spectralElements->setwavelength(wavelength.d, index);
	spectralElements->settell(tell.d, index);
	spectralElements->setrvel(rvel.d, index);
	spectralElements->setXCorrelation(crosscorrelation.d, index);
	spectralElements->setFlux(IFlux.d, index);
	spectralElements->setFluxVariance(IVar.d, index);
	spectralElements->setrawFlux(IFlux.d,index);
	spectralElements->setrawFluxVariance(IVar.d, index);
	spectralElements->setnormalizedFlux(normalizedFlux.d, index);
	spectralElements->setnormalizedFluxVariance(normalizedFluxFluxError.d, index);
	spectralElements->setfcalFlux(fcalFlux.d, index);
	spectralElements->setfcalFluxVariance(fcalFluxError.d, index);
	operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
	Polarimetry->setwavelength(wavelength.d, index);
	Polarimetry->setStokesParameter(StokesParameter, 0.0, 0.0, index);
	Polarimetry->setStokesParameter(StokesI, IFlux.d, IVar.d, index);
	Polarimetry->setDegreeOfPolarization(StokesParameter, DegPolarFlux.d, DegPolarVar.d, index);
	Polarimetry->setFirstNullPolarization(StokesParameter, FirstNullPolar.d, 0.0, index);
	Polarimetry->setSecondNullPolarization(StokesParameter, SecondNullPolar.d, 0.0, index);
	Polarimetry->setHasStokesI(true);
	if (StokesParameter == StokesQ) {
		Polarimetry->setHasStokesQ(true);
		Polarimetry->setHasDegreeOfStokesQ(true);
	}
	if (StokesParameter == StokesU) {
		Polarimetry->setHasStokesU(true);
		Polarimetry->setHasDegreeOfStokesU(true);
	}
	if (StokesParameter == StokesV) {
		Polarimetry->setHasStokesV(true);
		Polarimetry->setHasDegreeOfStokesV(true);
	}
}

/* 
 * void setWavelengthRangeFromLine(operaSpectralOrder *spectralOrder, string dataline)
 * \brief Takes a line from a wavelength range file and inserts the values into the order specified by that line.
 */
void operaIOFormats::setWavelengthRangeFromLine(operaSpectralOrder *spectralOrder, string dataline) {
	istringstream ss(dataline);
	unsigned order;
	double wl0, wlf;
	ss >> order >> wl0 >> wlf;
	spectralOrder->setminwavelength(wl0);
	spectralOrder->setmaxwavelength(wlf);
	spectralOrder->sethasWavelengthRange(true);
}

string operaIOFormats::getLineFromLibreEspritsp2Spectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << fixed << setprecision(4) << spectralElements->getwavelength(index) << ' '
	<< scientific << spectralElements->getFlux(index) << ' ' << sqrt(spectralElements->getFluxVariance(index));
	return ss.str();
}

string operaIOFormats::getLineFromLibreEspritsp1Spectrum(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	const operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
	const operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
	ss << fixed << setprecision(4) << spectralElements->getwavelength(index) << ' '
	<< scientific << spectralElements->getFlux(index) << ' ' << beamElements0->getFlux(index) << ' ' << beamElements1->getFlux(index) << ' '
	<< sqrt(spectralElements->getFluxVariance(index)) << ' ' << sqrt(beamElements0->getFluxVariance(index)) << ' ' << sqrt(beamElements1->getFluxVariance(index));
	return ss.str();
}

string operaIOFormats::getLineFromLibreEspritpolarimetry(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	const operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
	stokes_parameter_t stokesParameter = getStokesParameter(spectralOrder);
	ss << fixed << setprecision(4) << spectralElements->getwavelength(index) << ' '
	<< scientific << spectralElements->getFlux(index) << ' '
	<< Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(index) << ' '
	<< Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(index) << ' '
	<< Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(index) << ' '
	<< sqrt(Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(index));
	return ss.str();
}

string operaIOFormats::getLineFromLibreEspritSNR(const operaSpectralOrder *spectralOrder, unsigned index) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << fixed << setprecision(4) << spectralElements->getwavelength(index) << ' ' << scientific << spectralElements->getFluxSNR(index);
	return ss.str();
}

string operaIOFormats::getLineFromLibreEspritCenterSNR(const operaSpectralOrder *spectralOrder) {
	ostringstream ss;
	const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
	ss << fixed << setprecision(4) << spectralElements->getwavelength(spectralElements->getnSpectralElements()/2) << ' ' << scientific << spectralOrder->getCenterSNR();
	return ss.str();
}
