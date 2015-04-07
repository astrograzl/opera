/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaBarycentricWavelengthCorrection
 Version: 1.0
 Description: Apply Barycentric velocity wavelength correction 
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
#include "libraries/operaArgumentHandler.h"

template <typename T> inline int sign(const T& value) {
    return value < 0 ? -1 : value > 0;
}

//#ifndef TWO_PI
//#define TWO_PI 6.141592653589793238462643383279  //Wrong Value!! Change for TWOPI from operaLibCommon.h
//#endif

/*! \file operaBarycentricWavelengthCorrection.cpp */

using namespace std;

/*! 
 * operaBarycentricWavelengthCorrection
 * \author Eder Martioli & Lison Malo
 * \brief Calculate and apply Barycentric velocity wavelength correction.
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

void GenerateBarycentricWavelengthCorrectionPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display);

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
	string inputWaveFile;
	string outputRVelFile;
    double JDTime = 0.0;
    double MJDTime = 0.0;
    double etime=0.0;
    string observatory_coords_s;
    string object_coords_s;
    string ha_start_s;
    double observatory_elevation = 0.0;
	
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "input wavelength calibration file (.wcal)");
	args.AddRequiredArgument("outputRVelFile", outputRVelFile, "output radial velocity correction file (.rvel)");
    args.AddOptionalArgument("JDTime", JDTime, 0.0, "time in julian date");
    args.AddOptionalArgument("MJDTime", MJDTime, 0.0, "time in modified julian date");
    args.AddOptionalArgument("object_coords", object_coords_s, "0.0 0.0", "object sky coordinates \"RA Dec\"");
    args.AddOptionalArgument("observatory_coords", observatory_coords_s, "19:49:41.86 -155:28:18.00", "observatory geographic coordinates \"latitude longitude\"");
    args.AddOptionalArgument("observatory_elevation", observatory_elevation, 4200, "observatory elevation in meters");
    args.AddOptionalArgument("ha_start", ha_start_s, "00:00:00", "hour angle at start");
    args.AddOptionalArgument("etime", etime, 0.0, "exposure time (shutter open)");
	//" Example: "+string(modulename)+" ----inputWaveFile=/Users/lisonmalo/opera//calibrations/13BQ08-Nov25/OLAPAa_sp2_Normal.wcal.gz --observatory_coords=\"19:49:36 -155:28:18\" --object_coords=\"130.806167 3.398667\" --observatory_elevation=4207 --MJDTime=56622.6666429  --ha_start=1:17:50.06  --etime=30.0  --outputRVelFile=/Users/lisonmalo/opera//calibrations/13BQ08-Nov25/1671579i.rvel.gz -v\n\n"

	try {
		args.Parse(argc, argv);
		
		double ra = 0.0, dec = 0.0;
		skycoord_t object_coords = {0, 0, 0.0, 0, 0, 0.0};
		geocoord_t observatory_coords = {0, 0, 0.0, 0, 0, 0.0};
		hacoord_t ha_start = {0, 0, 0.0};
		
		if (!object_coords_s.empty()) {
			struct time_coord sexigesimal = {0,0,0};
			sscanf(object_coords_s.c_str(), "%lf %lf", &ra, &dec);
			ra = ra * (24.0/360.0);
			
			dec_to_sexigesimal(ra, &sexigesimal);
			object_coords.ra_h = (int)sexigesimal.hh;
			object_coords.ra_m = (int)sexigesimal.mm;
			object_coords.ra_s = sexigesimal.ss;
			dec_to_sexigesimal(dec, &sexigesimal);
			object_coords.dec_d = (int)sexigesimal.hh;
			object_coords.dec_m = (int)sexigesimal.mm;
			object_coords.dec_s = sexigesimal.ss;
		}
		if (!observatory_coords_s.empty()) {
			sscanf(observatory_coords_s.c_str(), "%d:%d:%lf %d:%d:%lf", &observatory_coords.latitude_d, &observatory_coords.latitude_m, &observatory_coords.latitude_s, &observatory_coords.longitude_d, &observatory_coords.longitude_m, &observatory_coords.longitude_s);
		}
		if (!ha_start_s.empty()) {
			sscanf(ha_start_s.c_str(), "%d:%d:%lf", &ha_start.ha_d, &ha_start.ha_m, &ha_start.ha_s);
		}
		
        cout.precision(6);
        cout << fixed;
        
		// we need an input wavelength calibration file ...
		if (inputWaveFile.empty()) {
			throw operaException("operaBarycentricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an outputRVelFile file ...        
		if (outputRVelFile.empty()) {
			throw operaException("operaBarycentricWavelengthCorrection: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
  
		if (args.verbose) {
			cout << "operaBarycentricWavelengthCorrection: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaBarycentricWavelengthCorrection: outputRVelFile = " << outputRVelFile << endl;
			cout << "operaBarycentricWavelengthCorrection: sky coordinates RA = " << object_coords.ra_h << ":" << object_coords.ra_m << ":"<< object_coords.ra_s  << " Dec=" << object_coords.dec_d << ":" << object_coords.dec_m << ":"<< object_coords.dec_s << "\n";
			cout << "operaBarycentricWavelengthCorrection: geographic coordinates Latitude = " << observatory_coords.latitude_d << ":" << observatory_coords.latitude_m << ":"<< observatory_coords.latitude_s  << " Longitude=" << observatory_coords.longitude_d << ":" << observatory_coords.longitude_m << ":"<< observatory_coords.longitude_s << "\n";
            cout << "operaBarycentricWavelengthCorrection: observatory_elevation = " << observatory_elevation << " m" << endl;
		}
        
		operaSpectralOrderVector spectralOrders(inputWaveFile);

        double rvCorrection, timeCorrection;

        double dbl_latitude = sign((double)(observatory_coords.latitude_d)) * (fabs((double)(observatory_coords.latitude_d)) + (double)(observatory_coords.latitude_m)/60.0 + (double)(observatory_coords.latitude_s)/3600.0);
        double dbl_longitude = sign((double)(observatory_coords.longitude_d)) * (fabs((double)(observatory_coords.longitude_d)) + (double)(observatory_coords.longitude_m)/60.0 + (double)(observatory_coords.longitude_s)/3600.0);

        double dbl_ha =sign((double)(ha_start.ha_d)) * (fabs((double)(ha_start.ha_d)) + (double)(ha_start.ha_m)/60.0 + (double)(ha_start.ha_s)/3600.0);
        
        double latitudeInRadians =  dbl_latitude*TWOPI/360; // Latitude (-PI/2, PI/2)
        double longitudeInRadians = dbl_longitude*TWOPI/360; // Longitude (-PI, PI)
        
        //if(longitudeInRadians < 0 && longitudeInRadians > -TWOPI/2.0) { // WEST
        //    longitudeInRadians = TWOPI + longitudeInRadians;
        //}
        
        if(args.debug) {
            cout << endl;
            cout << "operaBarycentricWavelengthCorrection: Latitude (deg, (+) North, (-) South) = " << dbl_latitude << endl;
            cout << "operaBarycentricWavelengthCorrection: Longitude (deg, +East) = " << longitudeInRadians * (360/TWOPI) << endl;

        }
        
        double earthEquatorialRadius;   //     a    double      equatorial radius (meters, Note 2)
        double earthFlattening;         //     f    double      flattening (Note 2)
        
        iauEform (WGS84,&earthEquatorialRadius,&earthFlattening);
        
        if(args.debug) {
            cout << endl;
            cout << "operaBarycentricWavelengthCorrection: (ref: IAU-SOFA WGS84) earthFlattening = " << earthFlattening  << endl; // f = 1.0 / 298.257223563 
            cout << "operaBarycentricWavelengthCorrection: (ref: IAU-SOFA WGS84) earthEquatorialRadius = " << earthEquatorialRadius/1000 << " km" << endl; // a = 6378137.0
        }
        
        double elong = longitudeInRadians; // elong   double     longitude (radians, east +ve)
        double phi = latitudeInRadians;     // phi   double     latitude (geodetic, radians, Note 4)
        double xyz[3];  // output geocentric coordinates
        
        iauGd2gce(earthEquatorialRadius,earthFlattening, elong, phi, observatory_elevation, xyz);
        
        if(args.debug) {
            cout << endl;            
            cout << "operaBarycentricWavelengthCorrection: x_geo = " << xyz[0]/1000 << " km" << endl;
            cout << "operaBarycentricWavelengthCorrection: y_geo = " << xyz[1]/1000 << " km" << endl;
            cout << "operaBarycentricWavelengthCorrection: z_geo = " << xyz[2]/1000 << " km" << endl;
        }
        
        if(args.debug) {
            cout << endl;            
            cout << "operaBarycentricWavelengthCorrection: Right Ascension (hr, 0-24) = " << ra <<  " hr" << endl;
            cout << "operaBarycentricWavelengthCorrection: Right Ascension (deg, +East) = " << ra*(360/24) <<  " deg" << endl;
        }
        
        double JDTime1 = 2400000.5;
        
        if(MJDTime == 0 && JDTime != 0) {
            MJDTime = JDTime - JDTime1;
        } else if (MJDTime != 0 && JDTime == 0) {
            JDTime = MJDTime + JDTime1;
        } else {
            throw operaException("operaBarycentricWavelengthCorrection: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        if (args.verbose) {
            cout << "operaBarycentricWavelengthCorrection: time in JD = " << JDTime << endl;
            cout << "operaBarycentricWavelengthCorrection: time in MJD = " << MJDTime << endl;
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
                
//        if(args.debug) {
//            cout << endl;
//            cout << "operaBarycentricWavelengthCorrection:  JDTime= " << JDTime << " JDTime1= " << JDTime1 << " mjd= " << mjd << endl;
//            cout << "operaBarycentricWavelengthCorrection:  UT = " << iyear << " / " << imonth << " / " << (double)iday + fractionday << endl;
//            cout << "operaBarycentricWavelengthCorrection:  deltat = " << deltat << endl;
//            cout << "operaBarycentricWavelengthCorrection:  TT1 = " << tt1 << " TT2 = " << tt2 << endl;
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
//            if(args.debug) {
//                struct time_coord gmst_sexigesimal = {0,0,0};
//                dec_to_sexigesimal(gmst, &gmst_sexigesimal);
//                cout << "operaBarycentricWavelengthCorrection: Greenwich Mean Sidereal Time (GMST, iauGmst00) = " << (int)gmst_sexigesimal.hh <<  ":" << (int)gmst_sexigesimal.mm << ":" << gmst_sexigesimal.ss<< endl;
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
//            if(args.debug) {
//                struct time_coord gmst_sexigesimal = {0,0,0};
//                dec_to_sexigesimal(gmst, &gmst_sexigesimal);
//                cout << "operaBarycentricWavelengthCorrection: Greenwich Mean Sidereal Time (GMST, formula) = " << (int)gmst_sexigesimal.hh <<  ":" << (int)gmst_sexigesimal.mm << ":" << gmst_sexigesimal.ss<< endl;
//            }
//        }

//        double ee00b = iauEe00b(JDTime1,mjd);
//        double gst = iauAnp(gmst*(TWOPI/24) + ee00b); // This is equivalent to routine iauGst00a
        
//        if(args.debug) {
//            struct time_coord gst_sexigesimal = {0,0,0};
//            dec_to_sexigesimal(gst*(24/TWOPI), &gst_sexigesimal);
            
//            cout << "operaBarycentricWavelengthCorrection: Greenwich apparent Sidereal Time (GST, iauGst06a) = " << (int)gst_sexigesimal.hh <<  ":" << (int)gst_sexigesimal.mm << ":" << gst_sexigesimal.ss<< endl;
//            cout << "operaBarycentricWavelengthCorrection: Greenwich apparent Sidereal Time (GST, iauGst06a) = " << gst*(24/TWOPI) <<  " hr" << endl;
//            cout << "operaBarycentricWavelengthCorrection: Greenwich apparent Sidereal Time (GST, iauGst06a) = " << gst*(360/TWOPI) <<  " deg" << endl;
//            cout << endl;
//        }
        
//        // Earth rotation angle (IAU 2000 model)
//        double era = iauEra00(JDTime1,mjd);
        
//        if(args.debug) {
//            cout << "operaBarycentricWavelengthCorrection: Earth rotation angle (ERA, iauEra00) = " << era*(24/TWOPI) <<  " hr" << endl;
//            cout << "operaBarycentricWavelengthCorrection: Earth rotation angle (ERA, iauEra00) = " << era*(360/TWOPI) <<  " deg" << endl;
//            cout << endl;
//        }
        
//        double lst = gst + longitudeInRadians;
        
//        struct time_coord lst_sexigesimal = {0,0,0};
//        dec_to_sexigesimal(lst*(24/TWOPI), &lst_sexigesimal);

//        if(args.debug) {
//            cout << "operaBarycentricWavelengthCorrection: Local Sidereal Time (LST) = " << (int)lst_sexigesimal.hh <<  ":" << (int)lst_sexigesimal.mm << ":" << lst_sexigesimal.ss<< endl;
//            cout << "operaBarycentricWavelengthCorrection: Local Sidereal Time (LST) := GST + longitude = " << lst*(24/TWOPI) <<  " hr" << endl;
//            cout << endl;
//        }
        
        //double ha = (lst - ra*(TWOPI/24)); // in radians
        double ha = dbl_ha + etime/3600.;
        struct time_coord ha_sexigesimal = {0,0,0};
        dec_to_sexigesimal(ha*(24/TWOPI), &ha_sexigesimal);
        
        if(args.debug) {
            cout << "operaBarycentricWavelengthCorrection: Hour Angle (HA) = " << (int)ha_sexigesimal.hh <<  ":" << (int)ha_sexigesimal.mm << ":" << ha_sexigesimal.ss<< endl;            
            cout << "operaBarycentricWavelengthCorrection: Hour Angle (HA):= LST - RA = " << ha*(24/TWOPI) <<  " hr" << endl;
            cout << endl;            
        }
      
        /*
         * Below we calculate the Heliocentric radial velocity correction using operaHelio library
         */
        heliocentric_correction(JDTime,ra,dec,ha,dbl_latitude,observatory_elevation, &timeCorrection, &rvCorrection);
        
     //   if(args.debug) {
            cout << "operaBarycentricWavelengthCorrection: TOTAL velocity correction (operaHelio library) = " << rvCorrection <<  " km/s" << endl;            
            cout << endl;            
      //  }

        /*
         * Below we implement SOFA routines to calculate the Barycentric
         *  radial velocity correction, instead of using the operaHelio library.
         */
        
        double raInRadians = ra/HRS_IN_RADIAN;
        double decInRadians = dec/DEG_IN_RADIAN;
        double haInRadians = ha/HRS_IN_RADIAN;
        /*
         * 1. The diurnal rotation of the Earth
         */
        
        // Angular rotation rate (2.PI/T)
        // According to IERS Numerical Standards (IAG 1999) w=7.2921150(1)e-5 rad/s
        double angularRotationRate = TWOPI/SIDEREALDAY_IN_SECONDS;

        if(args.debug)
            cout << "operaBarycentricWavelengthCorrection: angularRotationRate = " << angularRotationRate << " rad/s" << endl;
        
        // Geocentric radius at altitude H, in units of km
        double  geocentricRadius = sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2])/1000;
        
        if(args.debug)
            cout << "operaBarycentricWavelengthCorrection: geocentric radius = " << geocentricRadius << " km" << endl; 
        
        
        // Projected component of observer velocity to star
        // at declination = decInRadians and at hour angle = haInRadians, in units of km/s
        double diurnal_rvcorr = - angularRotationRate * geocentricRadius * cos(decInRadians) * sin(haInRadians);

        if(args.verbose) {
            cout << endl;
            cout << "\noperaBarycentricWavelengthCorrection: DIURNAL velocity correction (SOFA) = " << diurnal_rvcorr << " km/s" << endl;
        }
        
        /*
         * Below we will use the Earth position and velocity, Barycentric and barycentric,
         * with respect to the Barycentric Celestial Reference System given by SOFA library
         * to calculate the following components of the radial velocity:
         *
         * -> The motion of the Earth around the Earth-Moon barycenter (EMB)
         * -> The Barycentric motion of the EMB
         * -> The motion of the Sun around the center of mass of the solar system (CMSS)
         */
        /*
         *  Below it follows the definition of the position/velocity vectors given in 
         *  the SOFA library:
         *
         **        pvh[0][0]  x       }
         **        pvh[0][1]  y       } Barycentric position, AU
         **        pvh[0][2]  z       }
         **
         **        pvh[1][0]  xdot    }
         **        pvh[1][1]  ydot    } Barycentric velocity, AU/d
         **        pvh[1][2]  zdot    }
         **
         **        pvb[0][0]  x       }
         **        pvb[0][1]  y       } barycentric position, AU
         **        pvb[0][2]  z       }
         **
         **        pvb[1][0]  xdot    }
         **        pvb[1][1]  ydot    } barycentric velocity, AU/d
         **        pvb[1][2]  zdot    }
         **     The vectors are with respect to the Barycentric Celestial
         **     Reference System.  The time unit is one day in TDB.
         *
         */
         
         double pvh[2][3]; // Barycentric position(AU) / velocity (AU/d)
         double pvb[2][3]; // baricentric position(AU) / velocity (AU/d)
        
        iauEpv00(JDTime1,mjd,pvh,pvb);
        
        /*     
         * Below we calculate the projection of the Earth velocity vector on the direction of observation
         */
        double VxProjection = (pvb[1][0] * AU_IN_KILOMETERS/DAY_IN_SECONDS)*(cos(decInRadians)*cos(raInRadians));
        double VyProjection = (pvb[1][1] * AU_IN_KILOMETERS/DAY_IN_SECONDS)*(cos(decInRadians)*sin(raInRadians));
        double VzProjection = (pvb[1][2] * AU_IN_KILOMETERS/DAY_IN_SECONDS)*(sin(decInRadians));

        double barycentric_rvcorr = VxProjection + VyProjection + VzProjection;
        
		if(args.debug) {
			cout << "operaBarycentricWavelengthCorrection: alpha = " << raInRadians << endl;
			cout << "operaBarycentricWavelengthCorrection: delta = " << decInRadians << endl;
			cout << "operaBarycentricWavelengthCorrection: Vr.x = " << VxProjection << " km/s  " << endl;
			cout << "operaBarycentricWavelengthCorrection: Vr.y = " << VyProjection << " km/s  " << endl;
			cout << "operaBarycentricWavelengthCorrection: Vr.z = " << VzProjection << " km/s  " << endl;
		}
        
        if(args.verbose) {
            cout << endl;
            cout << "operaBarycentricWavelengthCorrection: BARYCENTRIC velocity correction (SOFA) = " << barycentric_rvcorr << " km/s  " << endl;
        }
        double sofaRVCorrection = barycentric_rvcorr + diurnal_rvcorr;
        
		if(args.verbose) {
            cout << endl;
            cout << "operaBarycentricWavelengthCorrection: TOTAL velocity correction (SOFA library) = " << sofaRVCorrection <<  " km/s" << endl;
		}
		
		// Set SOFA radial velocity correction         
		spectralOrders.setBarycentricRadialVelocityCorrection(sofaRVCorrection);
		// Uncomment below to set operaHelio radial velocity correction instead of SOFA
        // spectralOrders.setBarycentricRadialVelocityCorrection(rvCorrection);

		// output a rvCorrection file
		spectralOrders.WriteSpectralOrders(outputRVelFile, RVel);
        
	}
	catch (operaException e) {
		cerr << "operaBarycentricWavelengthCorrection: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaBarycentricWavelengthCorrection: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
