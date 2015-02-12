#ifndef OPERAPOLAR_H
#define OPERAPOLAR_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaPolar
 Version: 1.0
 Description: Create a polarization spectrum. 
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

/*! \brief Create a polarization spectrum. */
/*! \file operaPolar.h */
/*! \ingroup core */

/* prototypes */
static void printUsageSyntax(char *prgname);


/*!
 * \brief Plot the degree of polarization.
 * \details This function plots and creates the gnuplot script to plot the degree of polarization for every spectral element.
 * \details It can also display that plot after printing it.
 * \param gnuScriptFileName Output gnuplot script file name
 * \param outputPlotEPSFileName EPS plot file name
 * \param dataFileName Name of the data file holding the plot data
 * \param display Boolean value to display the plot on the screen
 * \param minorder minimum order to include in the plot range
 * \param maxorder maximum order to include in the plot range
 * \param StokesParameter A stokes_parameter_t value
 * \return void
 */
void GeneratePolarimetryPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, bool display, unsigned minorder, unsigned maxorder, stokes_parameter_t StokesParameter);

/*!
 * \brief Image plot of the degree of polarization or the flux for a given Stokes parameter.
 * \details This function plots and creates a gnuplot script to produce a 3D plot of the degree of polarization or the Stokes flux.
 * \details It can also display that plot after printing it.
 * \param gnuScriptFileName Output gnuplot script file name
 * \param outputPlotEPSFileName EPS plot file name
 * \param dataFileName Name of the data file holding the plot data
 * \param plotContinuum Boolean value to plot flux instead of degree of polarization
 * \param display Boolean value to display the plot on the screen
 * \param StokesParameter A stokes_parameter_t value
 * \return void
 */
void GeneratePolarization3DPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, bool plotContinuum, bool display, stokes_parameter_t StokesParameter);

#endif
