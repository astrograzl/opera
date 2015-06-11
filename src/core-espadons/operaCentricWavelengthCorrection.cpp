/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Library name: operaCentricWavelengthCorrection
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: May/2015
 
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

#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaHelio.h"					// for sexigesimal conversion
#include "libraries/operaArgumentHandler.h"
#include "../sofa/20120301_a/c/src/sofam.h"          // for sofa routines

template <typename T> inline int sign(const T& value) {
    return value < 0 ? -1 : value > 0;
}

using namespace std;

/*! 
 * operaCentricWavelengthCorrection
 * \author Eder Martioli & Lison Malo
 * \brief Calculate and apply Heliocentric or Barycentric velocity wavelength correction.
 * \file operaCentricWavelengthCorrection.cpp
 * \ingroup libraries
 */

int CentricWavelengthCorrection(int argc, char *argv[], const string modulename, const bool barycentric)
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
			throw operaException(modulename+": ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an outputRVelFile file ...        
		if (outputRVelFile.empty()) {
			throw operaException(modulename+": ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
  
		if (args.verbose) {
			cout << modulename << ": inputWaveFile = " << inputWaveFile << endl;
			cout << modulename << ": outputRVelFile = " << outputRVelFile << endl;
			cout << modulename << ": sky coordinates RA = " << object_coords.ra_h << ":" << object_coords.ra_m << ":"<< object_coords.ra_s  << " Dec=" << object_coords.dec_d << ":" << object_coords.dec_m << ":"<< object_coords.dec_s << "\n";
			cout << modulename << ": geographic coordinates Latitude = " << observatory_coords.latitude_d << ":" << observatory_coords.latitude_m << ":"<< observatory_coords.latitude_s  << " Longitude=" << observatory_coords.longitude_d << ":" << observatory_coords.longitude_m << ":"<< observatory_coords.longitude_s << "\n";
            cout << modulename << ": observatory_elevation = " << observatory_elevation << " m" << endl;
            cout << modulename << ": etime = " << etime << endl;
            cout << modulename << ": ha_start = " << ha_start.ha_d << ":" << ha_start.ha_m << ":" << ha_start.ha_s << endl;
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
            cout << modulename << ": Latitude (deg, (+) North, (-) South) = " << dbl_latitude << endl;
            cout << modulename << ": Longitude (deg, +East) = " << longitudeInRadians * (360/TWOPI) << endl;
        }
        
        double earthEquatorialRadius;   //     a    double      equatorial radius (meters, Note 2)
        double earthFlattening;         //     f    double      flattening (Note 2)
        
        iauEform (WGS84,&earthEquatorialRadius,&earthFlattening);
        
        if(args.debug) {
            cout << endl;
            cout << modulename << ": (ref: IAU-SOFA WGS84) earthFlattening = " << earthFlattening  << endl; // f = 1.0 / 298.257223563 
            cout << modulename << ": (ref: IAU-SOFA WGS84) earthEquatorialRadius = " << earthEquatorialRadius/1000 << " km" << endl; // a = 6378137.0
        }
        
        double elong = longitudeInRadians; // elong   double     longitude (radians, east +ve)
        double phi = latitudeInRadians;     // phi   double     latitude (geodetic, radians, Note 4)
        double xyz[3];  // output geocentric coordinates
        
        iauGd2gce(earthEquatorialRadius,earthFlattening, elong, phi, observatory_elevation, xyz);
        
        if(args.debug) {
            cout << endl;            
            cout << modulename << ": x_geo = " << xyz[0]/1000 << " km" << endl;
            cout << modulename << ": y_geo = " << xyz[1]/1000 << " km" << endl;
            cout << modulename << ": z_geo = " << xyz[2]/1000 << " km" << endl;
            cout << endl;            
            cout << modulename << ": Right Ascension (hr, 0-24) = " << ra <<  " hr" << endl;
            cout << modulename << ": Right Ascension (deg, +East) = " << ra*(360/24) <<  " deg" << endl;
        }
        
        double JDTime1 = 2400000.5;
        
        if(MJDTime == 0 && JDTime != 0) MJDTime = JDTime - JDTime1;
        else if (MJDTime != 0 && JDTime == 0) JDTime = MJDTime + JDTime1;
        else throw operaException(modulename+": ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
        
        if (args.verbose) {
            cout << modulename << ": time in JD = " << JDTime << endl;
            cout << modulename << ": time in MJD = " << MJDTime << endl;
        }
                
        //double ha = (lst - ra*(TWOPI/24)); // in radians
        double ha = dbl_ha + etime/3600.;
        struct time_coord ha_sexigesimal = {0,0,0};
        dec_to_sexigesimal(ha*(24/TWOPI), &ha_sexigesimal);
        
        if(args.debug) {
            cout << modulename << ": Hour Angle (HA) = " << (int)ha_sexigesimal.hh <<  ":" << (int)ha_sexigesimal.mm << ":" << ha_sexigesimal.ss<< endl;            
            cout << modulename << ": Hour Angle (HA):= LST - RA = " << ha*(24/TWOPI) <<  " hr" << endl;
            cout << endl;            
        }
      
        /*
         * Below we calculate the Heliocentric radial velocity correction using operaHelio library
         */
        heliocentric_correction(JDTime,ra,dec,ha,dbl_latitude,observatory_elevation, &timeCorrection, &rvCorrection);
        
        if(args.verbose) {
            cout << modulename << ": TOTAL velocity correction (operaHelio library) = " << rvCorrection <<  " km/s" << endl;
            cout << endl;
        }
        
        if(barycentric) {
			double mjd = MJDTime;
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

			if(args.debug) cout << modulename << ": angularRotationRate = " << angularRotationRate << " rad/s" << endl;
			
			// Geocentric radius at altitude H, in units of km
			double geocentricRadius = sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2])/1000;
			
			if(args.debug) cout << modulename << ": geocentric radius = " << geocentricRadius << " km" << endl; 
			
			
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
				cout << modulename << ": alpha = " << raInRadians << endl;
				cout << modulename << ": delta = " << decInRadians << endl;
				cout << modulename << ": Vr.x = " << VxProjection << " km/s  " << endl;
				cout << modulename << ": Vr.y = " << VyProjection << " km/s  " << endl;
				cout << modulename << ": Vr.z = " << VzProjection << " km/s  " << endl;
			}
			
			if(args.verbose) {
				cout << endl;
				cout << modulename << ": BARYCENTRIC velocity correction (SOFA) = " << barycentric_rvcorr << " km/s  " << endl;
			}
			double sofaRVCorrection = barycentric_rvcorr + diurnal_rvcorr;
			
			if(args.verbose) {
				cout << endl;
				cout << modulename << ": TOTAL velocity correction (SOFA library) = " << sofaRVCorrection <<  " km/s" << endl;
			}
			
			// Set SOFA radial velocity correction         
			spectralOrders.setRadialVelocityCorrection(sofaRVCorrection);
			// Uncomment below to set operaHelio radial velocity correction instead of SOFA
			// spectralOrders.setRadialVelocityCorrection(rvCorrection);
		}
        else spectralOrders.setRadialVelocityCorrection(rvCorrection);

		// output a rvCorrection file
		spectralOrders.WriteSpectralOrders(outputRVelFile, RVel);
	}
	catch (operaException e) {
		cerr << modulename << ": " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << modulename << ": " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
