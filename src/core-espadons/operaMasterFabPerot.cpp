/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterFabPerot
 Version: 1.0
 Description: This module creates a master Fabry-Perot image.
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

#include "core-espadons/operaMasterCalibration.h"

/*! \file operaMasterFabPerot.cpp */

/*!
 * operaMasterFabPerot
 * \author Doug Teeple
 * \brief Creates a master Align FITS image from a list of input Align FITS file names.
 * \arg argc
 * \arg argv
 * \note --images=[...]* --output=... [--pick=\<posint\>0\>] 
 * \note Pick one image or median stack images.
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */
int main(int argc, char *argv[])
{
	return MasterCalibrationCreation(argc, argv, "operaMasterFabPerot", "align");
}	
