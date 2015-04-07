# -*- coding: iso-8859-1 -*-
"""
Created on Oct 21 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Nov 03, 2014
"""
import sys
import espadons

############# MAIN FUNCTION TO EXECUTE PIPELINE ###############
def executePipeline(input, Dirs, config, keywords):

    """ 
    First assign input variables
    """
    night = input[0]
    intrumentmode = input[1]
    readoutspeed = input[2]
    clean = input[3]
    simulate = input[4]
    plot = input[5]
    verbose = input[6]
    trace = input[7]
    allowanyreadout = input[8]
    specificProduct = input[9]
    cleanall = input[10]
        
    """
    Set up instrument mode and readout mode
    """
    Instmode = espadons.InstMode(intrumentmode)
    Readmode = espadons.ReadoutMode(readoutspeed)
    
    """
    Set up default calibrations
    """
    DefaultCal = espadons.DefaultCalibration(Dirs, Instmode, Readmode)
    
    """
    CALIBRATION: Create calibration product file names, dependencies and command lines: 
    """
    productFilenames = espadons.setProductFilenames(Dirs,night,Instmode,Readmode,DefaultCal,keywords,allowanyreadout)
    plotFilenames = espadons.setPlotFilenames(Dirs,night,Instmode,Readmode,plot)
    productDependencies = espadons.setDependencies(productFilenames,Instmode)
    productCommands = espadons.setPipelineCommands(productFilenames,Dirs,night,Instmode,Readmode,keywords,config,plotFilenames,allowanyreadout,verbose)
   
    """
    Configure calibration product targets
    """
    products = espadons.Products(productFilenames,plotFilenames,productDependencies,productCommands,trace)
    
    """
    REDUCTION: Create object product file names, dependencies and command lines:
    """
    
    objproductFilenames = espadons.setObjectTargets(productFilenames, Dirs, night, Instmode, Readmode, DefaultCal, keywords, config, verbose)
    
    """
    Configure reduction product targets
    """
    products.addTargets(objproductFilenames)
        
    """
    Activate simulation mode
    """    
    if (simulate) :
        print "\n--- START SIMULATION ---\n"
        products.setSimulation()
    
    
    """
    ****** RUN PIPELINE ********
    """
    if (cleanall) :
        if(verbose and bool(simulate) == False) :
            print "Removing the following products: "
            products.displayTargets()
        
        products.cleanAll()
    else :
        #products.executeTarget("INSTRUMENTPROFILEPRODUCT") # This is an example how to execute a given target
        # Note that any target will trigger a novel process until the final product can be produced
            
        if (specificProduct in "LIBRE-ESPRIT") :
            products.executeTargetsWithSubstrings(["LESPC","LEPOL"])
        elif (specificProduct in "OBJECTS") :
            products.executeTargetsWithSubstrings(["OPSPC","OPPOL"])
        elif (specificProduct in "CALIBRATIONS") :
            products.executeTargetsWithSubstrings(["GEOMETRY","WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"])
        elif (specificProduct == "") :
            products.executeAllTargets()
        else :
            products.executeTargetsWithSubstrings([specificProduct])

        if clean :
            if (specificProduct in "LIBRE-ESPRIT") :
                products.removeTargets(["LESPC","LEPOL"])
            elif (specificProduct in "OBJECTS") :
                products.removeTargets(["OPSPC","OPPOL"])
            elif (specificProduct in "CALIBRATIONS") :
                products.removeTargets(["GEOMETRY","WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"])
            elif (specificProduct != "") :
                products.removeTargets([specificProduct])
            
            if(verbose and bool(simulate) == False) :
                print "Removing the following products: "
                products.displayTargets()

            products.cleanAll()
    #---------------------------
    
    """
    Reset simulation mode
    """
    if (simulate) :
        products.resetSimulation()
        print "\n --- END SIMULATION --- \n"

############# END ####################


