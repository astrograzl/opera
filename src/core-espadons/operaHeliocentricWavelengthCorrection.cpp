/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaHeliocentricWavelengthCorrection
 Version: 1.0
 Description: Apply Heliocentric velocity wavelength correction 
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
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"               // for anglecoord_t and timecoord_t
#include "libraries/operaHelio.h"					// for sexigesimal conversion
#include "../sofa/20120301_a/c/src/sofa.h"           // for sofa routines
#include "../sofa/20120301_a/c/src/sofam.h"          // for sofa routines
#include "core-espadons/operaHeliocentricWavelengthCorrection.h"

template <typename T> inline int sign(const T& value) {
    return value < 0 ? -1 : value > 0;
}

#ifndef NOTPROVIDED
#define NOTPROVIDED -999
#endif

//#ifndef TWO_PI
//#define TWO_PI 6.141592653589793238462643383279  //Wrong Value!! Change for TWOPI from operaLibCommon.h
//#endif

/*! \file operaHeliocentricWavelengthCorrection.cpp */

using namespace std;

/*! 
 * operaHeliocentricWavelengthCorrection
 * \author Eder Martioli & Lison Malo
 * \brief Calculate and apply Heliocentric velocity wavelength correction.
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


void GenerateHeliocentricWavelengthCorrectionPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display);

int main(int argc, char *argv[])
{
	int opt;
    
	string inputWaveFile;
	string outputRVelFile;
	   
    /*
     * Parameters for Heliocentric wavelength correction
     */
    double JDTime = 0.0;
    double MJDTime = 0.0;
    double etime=0.0;

    skycoord_t object_coords = {0, 0, 0.0, 0, 0, 0.0};
    double ra, dec;
    geocoord_t observatory_coords = {19, 49, 41.86, -155, 28, 18.00};
    hacoord_t ha_start = {0, 0, 0.0};

    double observatory_elevation = 4200;
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;
    bool maxorderprovided = false;
    
	bool debug=false, verbose=false, trace=false;
    
	struct option longopts[] = {
		{"inputWaveFile",				1, NULL, 'i'},    // input wavelength calibration file (.wcal)
        {"outputRVelFile",				1, NULL, 'o'},    // output radial velocity correction file (.rvel)
 		
		{"JDTime",						1, NULL, 'J'},    // time in julian date
		{"MJDTime",						1, NULL, 'D'},    // time in modified julian date
		{"object_coords",				1, NULL, 'K'},    // object sky coordinates
		{"observatory_coords",			1, NULL, 'G'},    // observatory geographic coordinates
		{"observatory_elevation",		1, NULL, 'E'},    // observatory elevation in meters
        {"ha_start",		            1, NULL, 'F'},    // Hour angle at start
        {"etime",		                1, NULL, 'L'},    // Exposure time (shutter open)

		{"ordernumber",					1, NULL, 'O'},
		{"minorder",					1, NULL, 'M'},
		{"maxorder",					1, NULL, 'X'},
        
		{"plot",						0, NULL, 'p'},
		{"verbose",						0, NULL, 'v'},
		{"debug",						0, NULL, 'd'},
		{"trace",						0, NULL, 't'},
		{"help",						0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:J:D:K:G:E:F:L:O:M:X:p::v::d::t::h",
							 longopts, NULL))  != -1)
	{
		switch(opt)
		{
			case 'i':
				inputWaveFile = optarg;
				break;
			case 'o':
				outputRVelFile = optarg;
				break;
                
            case 'J':		// time in julian date
				JDTime = atof(optarg);
				break;
            case 'D':		// time in modified julian date
				MJDTime = atof(optarg);
				break;
            case 'K':		// sky coordinates "RA Dec"
                if (strlen(optarg)) {
					struct time_coord sexigesimal = {0,0,0};
					sscanf(optarg, "%lf %lf", &ra, &dec);
                    ra = ra * (24.0/360.0);
                    
					dec_to_sexigesimal(ra, &sexigesimal);
					object_coords.ra_h = (int)sexigesimal.hh;
					object_coords.ra_m = (int)sexigesimal.mm;
					object_coords.ra_s = sexigesimal.ss;
					dec_to_sexigesimal(dec, &sexigesimal);
					object_coords.dec_d = (int)sexigesimal.hh;
					object_coords.dec_m = (int)sexigesimal.mm;
					object_coords.dec_s = sexigesimal.ss;
                    //sscanf(optarg, "%d:%d:%lf %d:%d:%lf", &object_coords.ra_h, &object_coords.ra_m, &object_coords.ra_s, &object_coords.dec_d, &object_coords.dec_m, &object_coords.dec_s);
				}
				break;
            case 'G':		// geographical coordinates "latitude longitude"
                if (strlen(optarg)) {
                    sscanf(optarg, "%d:%d:%lf %d:%d:%lf",&observatory_coords.latitude_d, &observatory_coords.latitude_m, &observatory_coords.latitude_s, &observatory_coords.longitude_d, &observatory_coords.longitude_m, &observatory_coords.longitude_s);
                }
                break;
			case 'E':
				observatory_elevation = atof(optarg);
				break;
            case 'F':		// HA value
                if (strlen(optarg)) {
                    sscanf(optarg, "%d:%d:%lf",&ha_start.ha_d, &ha_start.ha_m, &ha_start.ha_s);
                }    
                break;
            case 'L':
                etime = atof(optarg);
                break;
			case 'O':
				ordernumber = atoi(optarg);
				break;
			case 'M':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break;
			case 'v':
				verbose = true;
				break;
			case 'p':
				break;
			case 'd':
				debug = true;
				break;
			case 't':
				trace = true;
				break;
			case 'h':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}
	
	/*Start the module here*/
	
	try {
        cout.precision(6);
        cout << fixed;
        
		// we need an input wavelength calibration file ...
		if (inputWaveFile.empty()) {
			throw operaException("operaHeliocentricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an outputRVelFile file ...        
		if (outputRVelFile.empty()) {
			throw operaException("operaHeliocentricWavelengthCorrection: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
  
		if (verbose) {
			cout << "operaHeliocentricWavelengthCorrection: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaHeliocentricWavelengthCorrection: outputRVelFile = " << outputRVelFile << endl;
			cout << "operaHeliocentricWavelengthCorrection: sky coordinates RA = " << object_coords.ra_h << ":" << object_coords.ra_m << ":"<< object_coords.ra_s  << " Dec=" << object_coords.dec_d << ":" << object_coords.dec_m << ":"<< object_coords.dec_s << "\n";
			cout << "operaHeliocentricWavelengthCorrection: geographic coordinates Latitude = " << observatory_coords.latitude_d << ":" << observatory_coords.latitude_m << ":"<< observatory_coords.latitude_s  << " Longitude=" << observatory_coords.longitude_d << ":" << observatory_coords.longitude_m << ":"<< observatory_coords.longitude_s << "\n";
            cout << "operaHeliocentricWavelengthCorrection: observatory_elevation = " << observatory_elevation << " m" << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaHeliocentricWavelengthCorrection: ordernumber = " << ordernumber << endl;
            }
		}
        
		operaSpectralOrderVector spectralOrders(inputWaveFile);
        

        if(!minorderprovided) {
            minorder = spectralOrders.getMinorder();
        }
        if(!maxorderprovided) {
            maxorder = spectralOrders.getMaxorder();
        }
        
        if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}
		
        if (verbose)
			cout << "operaHeliocentricWavelengthCorrection: minorder = "<< minorder << " maxorder = " << maxorder << endl;

        double rvCorrection, timeCorrection;

        double dbl_latitude = sign((double)(observatory_coords.latitude_d)) * (fabs((double)(observatory_coords.latitude_d)) + (double)(observatory_coords.latitude_m)/60.0 + (double)(observatory_coords.latitude_s)/3600.0);
        double dbl_longitude = sign((double)(observatory_coords.longitude_d)) * (fabs((double)(observatory_coords.longitude_d)) + (double)(observatory_coords.longitude_m)/60.0 + (double)(observatory_coords.longitude_s)/3600.0);

        double dbl_ha =sign((double)(ha_start.ha_d)) * (fabs((double)(ha_start.ha_d)) + (double)(ha_start.ha_m)/60.0 + (double)(ha_start.ha_s)/3600.0);
        
        double latitudeInRadians =  dbl_latitude*TWOPI/360; // Latitude (-PI/2, PI/2)
        double longitudeInRadians = dbl_longitude*TWOPI/360; // Longitude (-PI, PI)
        
        //if(longitudeInRadians < 0 && longitudeInRadians > -TWOPI/2.0) { // WEST
        //    longitudeInRadians = TWOPI + longitudeInRadians;
        //}
        
        if(debug) {
            cout << endl;
            cout << "operaHeliocentricWavelengthCorrection: Latitude (deg, (+) North, (-) South) = " << dbl_latitude << endl;
            cout << "operaHeliocentricWavelengthCorrection: Longitude (deg, +East) = " << longitudeInRadians * (360/TWOPI) << endl;

        }
        
        double earthEquatorialRadius;   //     a    double      equatorial radius (meters, Note 2)
        double earthFlattening;         //     f    double      flattening (Note 2)
        
        iauEform (WGS84,&earthEquatorialRadius,&earthFlattening);
        
        if(debug) {
            cout << endl;
            cout << "operaHeliocentricWavelengthCorrection: (ref: IAU-SOFA WGS84) earthFlattening = " << earthFlattening  << endl; // f = 1.0 / 298.257223563 
            cout << "operaHeliocentricWavelengthCorrection: (ref: IAU-SOFA WGS84) earthEquatorialRadius = " << earthEquatorialRadius/1000 << " km" << endl; // a = 6378137.0
        }
        
        double elong = longitudeInRadians; // elong   double     longitude (radians, east +ve)
        double phi = latitudeInRadians;     // phi   double     latitude (geodetic, radians, Note 4)
        double xyz[3];  // output geocentric coordinates
        
        iauGd2gce(earthEquatorialRadius,earthFlattening, elong, phi, observatory_elevation, xyz);
        
        if(debug) {
            cout << endl;            
            cout << "operaHeliocentricWavelengthCorrection: x_geo = " << xyz[0]/1000 << " km" << endl;
            cout << "operaHeliocentricWavelengthCorrection: y_geo = " << xyz[1]/1000 << " km" << endl;
            cout << "operaHeliocentricWavelengthCorrection: z_geo = " << xyz[2]/1000 << " km" << endl;
        }
        
        if(debug) {
            cout << endl;            
            cout << "operaHeliocentricWavelengthCorrection: Right Ascension (hr, 0-24) = " << ra <<  " hr" << endl;
            cout << "operaHeliocentricWavelengthCorrection: Right Ascension (deg, +East) = " << ra*(360/24) <<  " deg" << endl;
        }
        
        double JDTime1 = 2400000.5;
        
        if(MJDTime == 0 && JDTime != 0) {
            MJDTime = JDTime - JDTime1;
        } else if (MJDTime != 0 && JDTime == 0) {
            JDTime = MJDTime + JDTime1;
        } else {
            throw operaException("operaHeliocentricWavelengthCorrection: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        if (verbose || debug) {
            cout << "operaHeliocentricWavelengthCorrection: time in JD = " << JDTime << endl;
            cout << "operaHeliocentricWavelengthCorrection: time in MJD = " << MJDTime << endl;
        }
        
        double mjd = MJDTime;

/* Commented part of the code not required anymore. (LMa Mar1,2015)
 */
//        int iyear, imonth, iday;
//        double fractionday;
//        iauJd2cal(JDTime1,mjd, &iyear, &imonth, &iday, &fractionday);
        
//        double deltat;
//        iauDat(iyear, imonth, iday,fractionday, &deltat);
        
//        double tt1, tt2;
//        // Time scale transformation:  Universal Time, UT1, to Terrestrial Time, TT.
//        iauUt1tt(JDTime1,mjd,deltat,&tt1,&tt2);
                
//        if(debug) {
//            cout << endl;
//            cout << "operaHeliocentricWavelengthCorrection:  JDTime= " << JDTime << " JDTime1= " << JDTime1 << " mjd= " << mjd << endl;
//            cout << "operaHeliocentricWavelengthCorrection:  UT = " << iyear << " / " << imonth << " / " << (double)iday + fractionday << endl;
//            cout << "operaHeliocentricWavelengthCorrection:  deltat = " << deltat << endl;
//            cout << "operaHeliocentricWavelengthCorrection:  TT1 = " << tt1 << " TT2 = " << tt2 << endl;
//            cout << endl;
//        }
        
        
//        double gmst = 0;
        
        /*
         * It seems like SOFA routines for GMST calculation gives a 
         * fairly different value than traditional formulae. For this reason, 
         * below we keep both ways for calculating the GMST. This will 
         * need to be checked.
         */
//        bool useSOFAiaugmst = true;

//        if(useSOFAiaugmst) {
//            gmst = iauGmst00(JDTime1,mjd,tt1,tt2)*(24/TWOPI);
//
//            if(debug) {
//                struct time_coord gmst_sexigesimal = {0,0,0};
//                dec_to_sexigesimal(gmst, &gmst_sexigesimal);
//                cout << "operaHeliocentricWavelengthCorrection: Greenwich Mean Sidereal Time (GMST, iauGmst00) = " << (int)gmst_sexigesimal.hh <<  ":" << (int)gmst_sexigesimal.mm << ":" << gmst_sexigesimal.ss<< endl;
//                cout << endl;
//            }
//            // Greenwich apparent sidereal time (radians)
//            // double gst = iauGst00a(JDTime1,mjd,tt1,tt2); // alternative way to calculate gst
//            // double gst = iauGst06a(JDTime1,mjd,tt1,tt2);
//
//        } else {
//            //  Below is a simple algorithm for computing apparent sidereal time obtained at USNO website.
//            //  http://aa.usno.navy.mil/faq/docs/GAST.php
//
//            double julian2000 = 2451545.0;
//            double julianCenturies = (JDTime - julian2000) / 36525;  /* centuries since J2000 */
//            gmst = 18.697374558 + 24.06570982441908 * julianCenturies * 36525;
//            gmst = ((gmst/24) - (double)floor(gmst/24))*24;
//            if(debug) {
//                struct time_coord gmst_sexigesimal = {0,0,0};
//                dec_to_sexigesimal(gmst, &gmst_sexigesimal);
//                cout << "operaHeliocentricWavelengthCorrection: Greenwich Mean Sidereal Time (GMST, formula) = " << (int)gmst_sexigesimal.hh <<  ":" << (int)gmst_sexigesimal.mm << ":" << gmst_sexigesimal.ss<< endl;
//            }
//        }

//        double ee00b = iauEe00b(JDTime1,mjd);
//        double gst = iauAnp(gmst*(TWOPI/24) + ee00b); // This is equivalent to routine iauGst00a
        
//        if(debug) {
//            struct time_coord gst_sexigesimal = {0,0,0};
//            dec_to_sexigesimal(gst*(24/TWOPI), &gst_sexigesimal);
            
//            cout << "operaHeliocentricWavelengthCorrection: Greenwich apparent Sidereal Time (GST, iauGst06a) = " << (int)gst_sexigesimal.hh <<  ":" << (int)gst_sexigesimal.mm << ":" << gst_sexigesimal.ss<< endl;
//            cout << "operaHeliocentricWavelengthCorrection: Greenwich apparent Sidereal Time (GST, iauGst06a) = " << gst*(24/TWOPI) <<  " hr" << endl;
//            cout << "operaHeliocentricWavelengthCorrection: Greenwich apparent Sidereal Time (GST, iauGst06a) = " << gst*(360/TWOPI) <<  " deg" << endl;
//            cout << endl;
//        }
        
//        // Earth rotation angle (IAU 2000 model)
//        double era = iauEra00(JDTime1,mjd);
        
//        if(debug) {
//            cout << "operaHeliocentricWavelengthCorrection: Earth rotation angle (ERA, iauEra00) = " << era*(24/TWOPI) <<  " hr" << endl;
//            cout << "operaHeliocentricWavelengthCorrection: Earth rotation angle (ERA, iauEra00) = " << era*(360/TWOPI) <<  " deg" << endl;
//            cout << endl;
//        }
        
//        double lst = gst + longitudeInRadians;
        
//        struct time_coord lst_sexigesimal = {0,0,0};
//        dec_to_sexigesimal(lst*(24/TWOPI), &lst_sexigesimal);

//        if(debug) {
//            cout << "operaHeliocentricWavelengthCorrection: Local Sidereal Time (LST) = " << (int)lst_sexigesimal.hh <<  ":" << (int)lst_sexigesimal.mm << ":" << lst_sexigesimal.ss<< endl;
//            cout << "operaHeliocentricWavelengthCorrection: Local Sidereal Time (LST) := GST + longitude = " << lst*(24/TWOPI) <<  " hr" << endl;
//            cout << endl;
//        }
        
        //double ha = (lst - ra*(TWOPI/24)); // in radians
        double ha = dbl_ha + etime/3600.;
        struct time_coord ha_sexigesimal = {0,0,0};
        dec_to_sexigesimal(ha*(24/TWOPI), &ha_sexigesimal);
        
        if(debug) {
            cout << "operaHeliocentricWavelengthCorrection: Hour Angle (HA) = " << (int)ha_sexigesimal.hh <<  ":" << (int)ha_sexigesimal.mm << ":" << ha_sexigesimal.ss<< endl;            
            cout << "operaHeliocentricWavelengthCorrection: Hour Angle (HA):= LST - RA = " << ha*(24/TWOPI) <<  " hr" << endl;
            cout << endl;            
        }
      
        /*
         * Below we calculate the Heliocentric radial velocity correction using operaHelio library
         */
        heliocentric_correction(JDTime,ra,dec,ha,dbl_latitude,observatory_elevation, &timeCorrection, &rvCorrection);
        
        if(verbose) {
            cout << "operaHeliocentricWavelengthCorrection: TOTAL velocity correction (operaHelio library) = " << rvCorrection <<  " km/s" << endl;
            cout << endl;
        }

        spectralOrders.setBarycentricRadialVelocityCorrection(rvCorrection);

		// output a rvCorrection file
		spectralOrders.WriteSpectralOrders(outputRVelFile, RVel);
        
	}
	catch (operaException e) {
		cerr << "operaHeliocentricWavelengthCorrection: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaHeliocentricWavelengthCorrection: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --inputWaveFile=<WAVE_FILE>"
	" --outputWaveFile=<WAVE_FILE>"
	" --object_coords=<\"RA_DBL_VALUE DEC_DBL_VALUE\">"
	" --observatory_coords=<\"UNS_VALUE:UNS_VALUE:DBL_VALUE UNS_VALUE:UNS_VALUE:DBL_VALUE\">"
	" --observatory_elevation=<DBL_VALUE>"
    " --JDTime=<DBL_VALUE>"
    " --MJDTime=<DBL_VALUE>"
    " --ha_start=<\"HA_DBL_VALUE\">"
	" --etime=<INT_VALUE>"
	" Example: "+string(modulename)+" ----inputWaveFile=/Users/lisonmalo/opera//calibrations/13BQ08-Nov25/OLAPAa_sp2_Normal.wcal.gz --observatory_coords=\"19:49:36 -155:28:18\" --object_coords=\"130.806167 3.398667\" --observatory_elevation=4207 --MJDTime=56622.6666429  --ha_start=1:17:50.06  --etime=30.0  --outputRVelFile=/Users/lisonmalo/opera//calibrations/13BQ08-Nov25/1671579i.rvel.gz -v\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file \n"
	"  -o, --outputWaveFile=<WAVE_FILE>, Output wavelength calibration file to store final solution\n"
	"  -J, --JDTime=<DBL_VALUE>, Time in julian date \n"
	"  -K, --object_coords=<\"RA_DBL_VALUE DEC_DBL_VALUE\">, Object sky coordinates \n"
	"  -G, --observatory_coords=<\"UNS_VALUE:UNS_VALUE:DBL_VALUE UNS_VALUE:UNS_VALUE:DBL_VALUE\">, Observatory geographic coordinates <latitude longitude>\n"
	"  -E, --observatory_elevation=<DBL_VALUE>, Observatory elevation above sea level in meters\n"
    "  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	"  -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n\n";
}

