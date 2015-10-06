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

#include "libraries/operaIOFormats.h"
#include "libraries/operaHelio.h"					// for sexigesimal conversion
#include "libraries/operaArgumentHandler.h"
#include "../sofa/20120301_a/c/src/sofam.h"          // for sofa routines
#include "core-espadons/operaCentricWavelengthCorrection.h"

using namespace std;

class Sexigesimal {
	public:
		Sexigesimal() : sign(1), h(0), m(0), s(0) { }
		bool SetFromString(string str);
		double ToDecimal() const;
		friend istream& operator>>(istream& in, Sexigesimal& var);
		friend ostream& operator<<(ostream& out, const Sexigesimal& var);
	private:
		int sign;
		int h;
		int m;
		double s;
};

bool Sexigesimal::SetFromString(string str) {
	if(str.empty()) return false;
	istringstream ss(str);
	if(str[0] == '-' || str[0] == '+') ss.get();
	int h2, m2;
	double s2;
	if(!(ss >> h2) || ss.get() != ':' || !(ss >> m2) || ss.get() != ':' || !(ss >> s2)) return false;
	sign = (str[0] == '-' ? -1 : 1);
	h = h2;
	m = m2;
	s = s2;
	return true;
}

double Sexigesimal::ToDecimal() const {
	return sign*h + m/60.0 + s/3600.0;
}

istream& operator>>(istream& in, Sexigesimal& var) {
	string temp;
	in >> temp;
	if(!var.SetFromString(temp)) in.setstate(ios_base::failbit);
	return in;
}

ostream& operator<<(ostream& out, const Sexigesimal& var) {
	if(var.sign < 0) out << '-';
	out << var.h << ":" << var.m << ":" << var.s;
	return out;
}

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
    string observatory_coords_s = "19:49:41.86 -155:28:18.00";
    string object_coords_s = "0.0 0.0";
    string ha_start_s = "00:00:00";
    double observatory_elevation = 4200;
	
    args.AddRequiredArgument("inputWaveFile", inputWaveFile, "input wavelength calibration file (.wcal)");
    args.AddRequiredArgument("outputRVelFile", outputRVelFile, "output radial velocity correction file (.rvel)");
    args.AddOptionalArgument("JDTime", JDTime, 0.0, "time in julian date (use either this or MJDTime)");
    args.AddOptionalArgument("MJDTime", MJDTime, 0.0, "time in modified julian date (use either this or JDTime)");
    args.AddRequiredArgument("object_coords", object_coords_s, "object sky coordinates \"RA Dec\"");
    args.AddRequiredArgument("observatory_coords", observatory_coords_s, "observatory geographic coordinates \"latitude longitude\"");
    args.AddRequiredArgument("observatory_elevation", observatory_elevation, "observatory elevation in meters");
    args.AddRequiredArgument("ha_start", ha_start_s, "hour angle at start");
    args.AddRequiredArgument("etime", etime, "exposure time (shutter open)");

	try {
		args.Parse(argc, argv);
		
		// we need an input wavelength calibration file ...
		if (inputWaveFile.empty()) {
			throw operaException(modulename+": ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an outputRVelFile file ...        
		if (outputRVelFile.empty()) {
			throw operaException(modulename+": ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		double ra = 0.0, dec = 0.0;
		Sexigesimal observatory_latitude, observatory_longitude, ha_start;
		
		if (!object_coords_s.empty()) {
			istringstream ss(object_coords_s);
			ss >> ra >> dec;
			ra = ra * (24.0/360.0); //convert right acension from degrees to hour angle
		}
		if (!observatory_coords_s.empty()) {
			istringstream ss(observatory_coords_s);
			ss >> observatory_latitude >> observatory_longitude;
		}
		if (!ha_start_s.empty()) {
			ha_start.SetFromString(ha_start_s);
		}
		
        cout.precision(6);
        cout << fixed;
  
		if (args.verbose) {
			cout << modulename << ": inputWaveFile = " << inputWaveFile << endl;
			cout << modulename << ": outputRVelFile = " << outputRVelFile << endl;
			cout << modulename << ": sky coordinates RA = " << ra << " hrs, Dec = " << dec << " deg" << endl;
			cout << modulename << ": geographic coordinates Latitude = " << observatory_latitude  << " Longitude = " << observatory_longitude << "\n";
            cout << modulename << ": observatory_elevation = " << observatory_elevation << " m" << endl;
            cout << modulename << ": etime = " << etime << endl;
            cout << modulename << ": ha_start = " << ha_start << endl;
		}
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);
		
		// Convert latitude, longitude, and ha from sexigesimal to decimal
        double dbl_latitude = observatory_latitude.ToDecimal();
        double dbl_longitude = observatory_longitude.ToDecimal();
        double dbl_ha = ha_start.ToDecimal();
        
        // Convert from MJDTime to JDTime or vice-versa
        if(MJDTime == 0 && JDTime != 0) MJDTime = JDTime - JDTime1;
        else if (MJDTime != 0 && JDTime == 0) JDTime = MJDTime + JDTime1;
        else throw operaException(modulename+": ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        
        // Move time and HA to the middle of the exposure
        const double midexp = etime/7200.0; //Convert from seconds to hours, then divide by 2 to get the middle of the exposure
		JDTime += midexp;
        MJDTime += midexp;
        dbl_ha += midexp;
        if (args.verbose) {
            cout << modulename << ": mid-exposure etime = " << etime/2.0 << endl;
            cout << modulename << ": time in JD = " << JDTime << endl;
            cout << modulename << ": time in MJD = " << MJDTime << endl;
            cout << modulename << ": Hour Angle (HA) = " << dbl_ha << " hr" << endl;
        }
        
        if(barycentric) {
			// Calculate the radial velocity correction using the SOFA library
			double sofaRVCorrection = CalculateSOFABarycentricCorrection(dbl_latitude, dbl_longitude, observatory_elevation, MJDTime, ra, dec, dbl_ha, args.debug, modulename);
			if(args.verbose) {
				cout << modulename << ": TOTAL velocity correction (SOFA library) = " << sofaRVCorrection <<  " km/s" << endl;
				cout << endl;
			}
			spectralOrders.setRadialVelocityCorrection(sofaRVCorrection); // Set SOFA radial velocity correction         
		}
        else {
			// Calculate the Heliocentric radial velocity correction using operaHelio library
			double rvCorrection, timeCorrection;
			heliocentric_correction(JDTime,ra,dec,dbl_ha,dbl_latitude,observatory_elevation, &timeCorrection, &rvCorrection);
			if(args.verbose) {
				cout << modulename << ": TOTAL velocity correction (operaHelio library) = " << rvCorrection <<  " km/s" << endl;
				cout << endl;
			}
			spectralOrders.setRadialVelocityCorrection(rvCorrection); // set operaHelio radial velocity correction instead of SOFA
		}

		// output rvCorrection file
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputRVelFile, RVel);
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

double CalculateSOFABarycentricCorrection(double latitudeInDeg, double longitudeInDeg, double elevation, double MJDTime, double raInHrs, double decInDeg, double haInHrs, bool debug, string modulename)
{
	double latitudeInRadians =  latitudeInDeg*TWOPI/360; // Latitude (-PI/2, PI/2)
	double longitudeInRadians = longitudeInDeg*TWOPI/360; // Longitude (-PI, PI)
	
	//if(longitudeInRadians < 0 && longitudeInRadians > -TWOPI/2.0) { // WEST
	//    longitudeInRadians = TWOPI + longitudeInRadians;
	//}
	
	if(debug) {
		cout << endl;
		cout << modulename << ": Latitude (deg, (+) North, (-) South) = " << latitudeInDeg << endl;
		cout << modulename << ": Longitude (deg, +East) = " << longitudeInDeg << endl;
	}
	
	double earthEquatorialRadius;   //     a    double      equatorial radius (meters, Note 2)
	double earthFlattening;         //     f    double      flattening (Note 2)
	
	iauEform (WGS84,&earthEquatorialRadius,&earthFlattening);
	
	if(debug) {
		cout << endl;
		cout << modulename << ": (ref: IAU-SOFA WGS84) earthFlattening = " << earthFlattening  << endl; // f = 1.0 / 298.257223563 
		cout << modulename << ": (ref: IAU-SOFA WGS84) earthEquatorialRadius = " << earthEquatorialRadius/1000 << " km" << endl; // a = 6378137.0
	}
	
	double elong = longitudeInRadians; // elong   double     longitude (radians, east +ve)
	double phi = latitudeInRadians;     // phi   double     latitude (geodetic, radians, Note 4)
	double xyz[3];  // output geocentric coordinates
	
	iauGd2gce(earthEquatorialRadius,earthFlattening, elong, phi, elevation, xyz);
	
	if(debug) {
		cout << endl;            
		cout << modulename << ": x_geo = " << xyz[0]/1000 << " km" << endl;
		cout << modulename << ": y_geo = " << xyz[1]/1000 << " km" << endl;
		cout << modulename << ": z_geo = " << xyz[2]/1000 << " km" << endl;
		cout << endl;            
		cout << modulename << ": Right Ascension (hr, 0-24) = " << raInHrs <<  " hr" << endl;
		cout << modulename << ": Right Ascension (deg, +East) = " << raInHrs*(360/24) <<  " deg" << endl;
	}
	
	double mjd = MJDTime;
	/*
	 * Below we implement SOFA routines to calculate the Barycentric
	 *  radial velocity correction, instead of using the operaHelio library.
	 */
	
	double raInRadians = raInHrs/HRS_IN_RADIAN;
	double decInRadians = decInDeg/DEG_IN_RADIAN;
	double haInRadians = haInHrs/HRS_IN_RADIAN;
	/*
	 * 1. The diurnal rotation of the Earth
	 */
	
	// Angular rotation rate (2.PI/T)
	// According to IERS Numerical Standards (IAG 1999) w=7.2921150(1)e-5 rad/s
	double angularRotationRate = TWOPI/SIDEREALDAY_IN_SECONDS;

	if(debug) cout << modulename << ": angularRotationRate = " << angularRotationRate << " rad/s" << endl;
	
	// Geocentric radius at altitude H, in units of km
	double geocentricRadius = sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2])/1000;
	
	if(debug) cout << modulename << ": geocentric radius = " << geocentricRadius << " km" << endl; 
	
	
	// Projected component of observer velocity to star
	// at declination = decInRadians and at hour angle = haInRadians, in units of km/s
	double diurnal_rvcorr = - angularRotationRate * geocentricRadius * cos(decInRadians) * sin(haInRadians);

	if(debug) {
		cout << endl;
		cout << modulename << ": DIURNAL velocity correction (SOFA) = " << diurnal_rvcorr << " km/s" << endl;
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
	
	if(debug) {
		cout << modulename << ": alpha = " << raInRadians << endl;
		cout << modulename << ": delta = " << decInRadians << endl;
		cout << modulename << ": Vr.x = " << VxProjection << " km/s  " << endl;
		cout << modulename << ": Vr.y = " << VyProjection << " km/s  " << endl;
		cout << modulename << ": Vr.z = " << VzProjection << " km/s  " << endl;
	}
	
	if(debug) {
		cout << modulename << ": BARYCENTRIC velocity correction (SOFA) = " << barycentric_rvcorr << " km/s  " << endl;
	}
	return barycentric_rvcorr + diurnal_rvcorr;
}
