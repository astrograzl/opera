#ifndef OPERAQUERYIMAGEINFO_H
#define OPERAQUERYIMAGEINFO_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaQueryImageInfo
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

/*! \brief Query image information. */
/*! \file operaQueryImageInfo.h */
/*! \ingroup tools */

#define MAXCONFIGVALUES 100
#define MAXNDIRS 100
#define MAXNFILES 1000
#define MAXDIRNAMESIZE 1000
#define MAX_MODE_CHANGES 24
#define MAX_ETYPES_CONFIG 7

/* prototypes */
unsigned operaQueryImageInfoConfigurationAccess(char *par, char *values[], unsigned maxvalues, operaErrorCode *errorcode);
unsigned SplitValues(char *parvalues, char *values[], unsigned maxvalues);
void cleanFITSHeaderValue(char *value_from_FITS,char *output_clean_value);

#endif


