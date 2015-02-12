#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
"""
    Created on Oct 21 2014

    Description: A wrapper to run Opera ESPaDOnS reduction pipeline.
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    ./opera.py --night=13BQ04-Sep20 -pvt
"""

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import sys,os
import espadonspipeline
import espadons

parser = OptionParser()
parser.add_option("-N", "--night", dest="night", help="night directory",type='string',default="")
parser.add_option("-T", "--product", dest="product", help='target product: "CALIBRATIONS", "OBJECTS" (default) or "LIBRE-ESPRIT"',type='string',default="OBJECTS")
parser.add_option("-a", action="store_true", dest="cleanall", help="JUST clean all products",default=False)
parser.add_option("-c", action="store_true", dest="clean", help="clean products",default=False)
parser.add_option("-s", action="store_true", dest="simulate", help="simulate",default=False)
parser.add_option("-p", action="store_true", dest="plot", help="plots",default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose",default=False)
parser.add_option("-t", action="store_true", dest="trace", help="trace",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with opera.py -h ";sys.exit(1);

pipelinehomedir = "/Users/edermartioli/opera-1.0/"
datarootdir ="/data/espadons/"
productrootdir = "/Users/edermartioli/Reductions/Espadons/"

if options.verbose:
    print 'PIPELINE HOME DIR: ', pipelinehomedir
    print 'DATA ROOT DIR: ', datarootdir
    print 'PRODUCT ROOT DIR: ', productrootdir
    print 'NIGHT: ', options.night

"""
Set up directories:
"""
Dirs = espadons.Directories(pipelinehomedir,datarootdir,productrootdir,options.night)
"""
Set up config files:
"""
config = espadons.ConfigFiles(Dirs)

"""
Set up Espadons keywords:
"""
keywords = espadons.Keywords()

"""
Set up modes available for reduction:
"""
allowanyreadout = False

forcecalibration = False # This is set to "True" for calibrations even when there is no object files
if(options.product == "CALIBRATIONS") :
    forcecalibration = True

modes = espadons.ReductionModes(Dirs, keywords, allowanyreadout, forcecalibration)

modes.displayOverallStats()

for mode in modes.getInstReadModes() :
    intrumentmode = mode[0]
    readoutspeed = mode[1]
    modes.displayModeStats(intrumentmode,readoutspeed)
    input = [options.night,intrumentmode,readoutspeed,options.clean,options.simulate,options.plot,options.verbose,options.trace,allowanyreadout,options.product,options.cleanall]
    espadonspipeline.executePipeline(input, Dirs, config, keywords)

modes.cleanModes()



