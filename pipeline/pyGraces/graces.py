# -*- coding: utf-8 -*-
"""
Created on May 28 2015
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Jul 24 2015
"""
import sys, os
import subprocess
#import pipelinebrain
#import gracesUtils

########## DIRECTORIES ############
class Directories :
    
    'Common base class for directories'
    def __init__(self, pipelinehomedir,datarootdir,productrootdir,night):
        self.OP_HOME=pipelinehomedir
        self.DATAROOTDIR=datarootdir
        self.PRODUCTROOTDIR=productrootdir
        self.PIPELINEDIR=self.OP_HOME+"/pipeline/"
        self.DEFAULTCALIBRATIONDIR=self.OP_HOME+"/DefaultCalibration/"
        self.CONFIGDIR=self.OP_HOME+"/config/"
        self.STANDARDSDIR=self.CONFIGDIR+"/standardstars/"
        self.EXE=self.OP_HOME+"/bin/"
        self.DATADIR=self.DATAROOTDIR+"/"+night+"/"
        self.PRODUCTDIR=self.PRODUCTROOTDIR+"/"+night+"/"
        
        self.createNightDir(self.PRODUCTDIR)
    
    def createNightDir(self, inputdir):
        if os.path.exists(inputdir) :
            print "Product directory: ",inputdir," (already exists)"
        else :
            print "Product directory: ",inputdir," (created new directory)"
            createdircommand = "mkdir " + inputdir
            os.system(createdircommand)

##################################

######### CONFIG FILES ###########
class ConfigFiles :
    'Common base class for config files'
    standarddata = {}
    
    def __init__(self, Dirs):
        self.BADPIXELMASK = Dirs.CONFIGDIR + "badpix_olapa-a.fits.fz"
        self.THARATLASLINES = Dirs.CONFIGDIR + "thar_MM201006.dat.gz"
        self.THARATLASSPECTRUM = Dirs.CONFIGDIR + "LovisPepe_ThArAtlas.dat.gz"
        self.WAVEFIRSTGUESS = Dirs.CONFIGDIR + "wcal_ref.dat.gz"
        self.SOLARTYPEWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
        #self.ATYPEWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
        self.ATYPEWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForUncalContinuumDetection.txt"
        self.ATLASWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForRefContinuumDetection.txt"
        self.TELLURICWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForTelluricAbsorption.txt"
        self.TELLURICLINES = Dirs.CONFIGDIR + "opera_HITRAN08-extracted.par.gz"
        self.TELLURICSPECTRUM = Dirs.CONFIGDIR + "KPNO_atmtrans.dat.gz"
        self.LEORDERWAVELENGTH = Dirs.CONFIGDIR + "LE_order_wavelength.dat"
        
        #self.STANDARDLISTFILE = Dirs.STANDARDSDIR + "operaStandardStars.dat"
        self.STANDARDLISTFILE = Dirs.STANDARDSDIR + "operaStandardStarsForGRACES.dat"
        self.readStandards(Dirs,self.STANDARDLISTFILE)
        
        self.SYNTHETICSPECTRUM = Dirs.CONFIGDIR + "/syntethicSpectra/Teff5500.spec"
        self.RVXCORRWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForRVxcorr.txt"
        
        #self.OLAPAFLATRESPONSE = Dirs.CONFIGDIR + "flat_resp_olapa.s"
        self.OLAPAFLATRESPONSE = Dirs.CONFIGDIR + "flat_resp_olapa.fits.gz"
        self.EEV1FLATRESPONSE = Dirs.CONFIGDIR + "flat_resp_eev1.s"
    
    # Function below reads list of standard stars for which there is available calibration data
    def readStandards(self, Dirs, stdListFilename) :
        stdlistfile = open(stdListFilename)
        for line in stdlistfile :
            if(line[0] != '#') :
                name = line.split(' ')[0]
                datafilename = Dirs.STANDARDSDIR + name + "_operaFluxCal.dat"
                ra = line.split(' ')[1]
                dec = line.split(' ')[2]
                pmra = line.split(' ')[3]
                pmdec = line.split(' ')[4]
                vmag = line.split(' ')[5]
                teff = line.split(' ')[6]
                rv = line.split(' ')[7]
                stype = line.split(' ')[8]
                if os.path.exists(datafilename) :
                    self.standarddata[name] = [datafilename,ra,dec,pmra,pmdec,vmag,teff,rv,stype]
        
        stdlistfile.close()
    
    # Function to check if a given standard star exist
    def hasStandard(self, StandardName) :
        if StandardName in self.standarddata :
            return True
        else :
            return False

    # Function below returns the corresponding data for a standard star in the list
    def getStandardDataFile(self, key) :
        return self.standarddata.get(key, 0)[0]
    
    def getStandardRA(self, key) :
        return self.standarddata.get(key, 0)[1]
    
    def getStandardDec(self, key) :
        return self.standarddata.get(key, 0)[2]
    
    def getStandardPMRA(self, key) :
        return self.standarddata.get(key, 0)[3]
    
    def getStandardPMDec(self, key) :
        return self.standarddata.get(key, 0)[4]
    
    def getStandardVmag(self, key) :
        return self.standarddata.get(key, 0)[5]
    
    def getStandardTeff(self, key) :
        return self.standarddata.get(key, 0)[6]
    
    def getStandardRV(self, key) :
        return self.standarddata.get(key, 0)[7]
    
    def getStandardSType(self, key) :
        return self.standarddata.get(key, 0)[8]
#-------------------------------------------
##################################

######### KEYWORDS ###########
class Keywords :
    'Common base class for keywords'
    def __init__(self):
        self.BIASKEYWORD="BIAS"
        self.FLATKEYWORD="FLAT"
        self.COMPKEYWORD="COMPARISON"
        self.OBJECTKEYWORD="OBJECT"
##############################

######### DEFAULT CALIBRATION ###########
class DefaultCalibration :    
    'Common base class for default calibration'
    def __init__(self, Dirs, Instmode, Readmode):
        self.DEFAULTMASTERBIAS = Dirs.DEFAULTCALIBRATIONDIR + Readmode.READOUTSPEEDSHORTNAME + "_masterbias.fits.gz"
        self.DEFAULTGAIN = Dirs.DEFAULTCALIBRATIONDIR + Readmode.READOUTSPEEDSHORTNAME + ".gain.gz"
        self.DEFAULTMASTERCOMP = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_mastercomp.fits.gz"
        self.DEFAULTMASTERFLAT = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_masterflat.fits.gz"
        self.DEFAULTMASTERFLUXCAL = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME + ".fcal.gz"
##################################

######### INSTRUMENT MODE ###########
class InstMode :
    'Common base class for all instrument modes'
    mode = 0
    def __init__(self, intrumentmode):
        self.mode = int(intrumentmode)

        self.INSTRUME="GRACES"
        self.WAVELENGTHFORNORMALIZATION=548
        self.APERTUREHEIGHT = 0.6923
        
        if (self.mode == 1) :
            self.INSTRUMENTMODESHORTNAME="StarOnly"
            self.INSTRUMENTMODEKEY="FOURSLICE"
            self.STARPLUSKYMODEFLAG=0
            self.INVERTSKYFIBERFLAG=""
            self.SPCMODULE="operaStarOnly"
            self.recenterIPUsingSliceSymmetry="0"
            self.NUMBEROFSLICES="4"
            self.NUMBEROFBEAMS="1"
            self.SPECTRALRESOLUTION="60000"
            self.LESPECTRUMTYPE=21
            self.ORDSPCAPERTURE=32
            self.SPACINGREFERENCEORDERNUMBER=47
            self.SPACINGREFERENCEORDERSEPARATION=55
            self.GEOMAPERTURE=30
            self.GEOMMAXNORDERS=40
            self.GEOMMINORDERTOUSE=21
            self.REFERENCELINEWIDTH=1.8
            self.IPXSIZE="34"
            self.IPYSIZE="7"
            self.IPXSAMPLING="3"
            self.IPYSAMPLING="7"
            self.TILTANGLE="-2.0"
            self.CONSTANTTILTFLAG="1"
            self.APERAPERTURE=32
            self.APERGAP=0
            self.WAVEUNCALLINEWIDTH=1.6
            self.RADIALVELOCITYRANGE=5.0
            self.RADIALVELOCITYSTEP=0.15
            self.RADIALVELOCITYSEARCHRANGE=200
            self.RADIALVELOCITYSEARCHSTEP=0.5
            self.GRACESFLATRESPONSE = "graces_StarOnly_flatresp.fits.gz"

        elif (self.mode == 2) :
            self.INSTRUMENTMODESHORTNAME="StarPlusSky"
            self.INSTRUMENTMODEKEY="TWOSLICE"
            self.STARPLUSKYMODEFLAG=1
            self.INVERTSKYFIBERFLAG="--starplusskyInvertSkyFiber=1"
            self.SPCMODULE="operaStarPlusSky"
            self.recenterIPUsingSliceSymmetry="1"
            self.NUMBEROFSLICES="4"
            self.NUMBEROFBEAMS="2"
            self.SPECTRALRESOLUTION="40000"
            self.LESPECTRUMTYPE=20
            self.ORDSPCAPERTURE=32
            self.SPACINGREFERENCEORDERNUMBER=45
            self.SPACINGREFERENCEORDERSEPARATION=52
            self.GEOMAPERTURE=32
            self.GEOMMAXNORDERS=40
            self.GEOMMINORDERTOUSE=21
            self.REFERENCELINEWIDTH=2.5
            self.IPXSIZE="32"
            self.IPYSIZE="5"
            self.IPXSAMPLING="5"
            self.IPYSAMPLING="5"
            self.TILTANGLE="-1.0"
            self.CONSTANTTILTFLAG="1"
            self.APERAPERTURE=32
            self.APERGAP=0
            self.WAVEUNCALLINEWIDTH=1.8
            self.RADIALVELOCITYRANGE=5.0
            self.RADIALVELOCITYSTEP=0.15
            self.RADIALVELOCITYSEARCHRANGE=200
            self.RADIALVELOCITYSEARCHSTEP=0.5
            self.GRACESFLATRESPONSE = "graces_StarPlusSky_flatresp.fits.gz"

##################################

######### READOUT MODE ###########
class ReadoutMode :
    'Common base class for all readout modes'
    mode = 0
    def __init__(self, readoutmode):
        self.mode = int(readoutmode)
        if (self.mode == 1) :
            self.READOUTSPEEDSHORTNAME="fast"
            self.READOUTSPEED="Fast: 4.70e noise, 1.60e/ADU, 32s"
            self.DEFAULTGAIN="1.6"
            self.DEFAULTNOISE="4.14"
        elif (self.mode == 2) :
            self.READOUTSPEEDSHORTNAME="normal"
            self.READOUTSPEED="Normal: 4.20e noise, 1.30e/ADU, 38s"
            self.DEFAULTGAIN="1.3"
            self.DEFAULTNOISE="3.8"
        elif (self.mode == 3) :
            self.READOUTSPEEDSHORTNAME="slow"
            self.READOUTSPEED="Slow: 2.90e noise, 1.20e/ADU, 60s"
            self.DEFAULTGAIN="1.1"
            self.DEFAULTNOISE="2.98"
##################################

######### OBJECT PRODUCT FILE NAMES, DEPENDENCIES, AND  COMMAND LINES ###########
def setObjectTargets(products, Dirs, night, Instmode, Readmode, DefaultCal, Keywords, config, verbose) :
    
    emptystring = ""
    
    verstr = ""
    if (verbose) :
        verstr = " --verbose"
    
    objproducts = {}
    objdependencies = {}
    objcommands = {}
    
    commandline = ObjectListCommand(Dirs, Instmode, Readmode, Keywords.OBJECTKEYWORD)
    objectfiles = subprocess.check_output(commandline,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
    objfilelist = objectfiles.split()

    polquad = []
    polextkey = []
    polcount = 0
    pkey = ""
    ptarget = ""
    
    listofstdfcalkey = []
    listofstdfcaltargets = []
    
    ### Create TARGET for master flux calibration ######
    #masterfcalkey = "MASTERFLUXCALIBRATION"
    #INSTCONFIGPREFIX = Dirs.PRODUCTDIR + night + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
    #masterfcaltarget = INSTCONFIGPREFIX + ".fcal.gz"
    #objproducts[masterfcalkey] = masterfcaltarget
    ####################################################
    
    for file in objfilelist:
        ## Caution: the action below strongly depends on the file name format!!
        if "fits.gz" in file :
            basename = file[-22:-8]
        else :
            basename = file[-19:-5]
        
        ### Get header data info #########
        objectcommand = Dirs.EXE +"operagetheader --keyword=OBJECT " + file
        objectname = subprocess.check_output(objectcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        exptimecommand = Dirs.EXE +"operagetheader --keyword=EXPTIME " + file
        exptime = subprocess.check_output(exptimecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        racommand = Dirs.EXE +"operagetheader --keyword=RA_DEG " + file
        absra_center = subprocess.check_output(racommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        deccommand = Dirs.EXE +"operagetheader --keyword=DEC_DEG " + file
        absdec_center = subprocess.check_output(deccommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        mjcommand = Dirs.EXE +"operagetheader --keyword=MJDATE " + file
        mjdate = subprocess.check_output(mjcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        airmasscommand = Dirs.EXE +"operagetheader --keyword=AIRMASS " + file
        airmass = subprocess.check_output(airmasscommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        expnumcommand = Dirs.EXE +"operagetheader --keyword=EXPNUM " + file
        expnum = subprocess.check_output(expnumcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        outsideTemperaturecommand = Dirs.EXE +"operagetheader --keyword=TEMPERAT " + file
        outsideTemperature = subprocess.check_output(outsideTemperaturecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        windspeedcommand = Dirs.EXE +"operagetheader --keyword=WINDSPED " + file
        windspeed = subprocess.check_output(windspeedcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        relativeHumididycommand = Dirs.EXE +"operagetheader --keyword=RELHUMID " + file
        relativeHumididy = subprocess.check_output(relativeHumididycommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        atmoPressurecommand = Dirs.EXE +"operagetheader --keyword=PRESSURE " + file
        atmoPressure = subprocess.check_output(atmoPressurecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        startHAcommand = Dirs.EXE +"operagetheader --keyword=HA " + file
        startHA = subprocess.check_output(startHAcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        ########################################


        ### Create TARGETS for Extraction #########
        extkey = "EXTRACT" + basename
        exttarget = Dirs.PRODUCTDIR + basename + ".e.gz"
        objproducts[extkey] = exttarget
        objdependencies[extkey] = ["MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
        objcommands[extkey] = objectExtractionCommand(Dirs, exttarget, file, products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"], Instmode) + verstr
        ########################################

        ### Create TARGETS for telluric wavelength correction #########
        tellkey = "TELL" + basename
        telltarget = Dirs.PRODUCTDIR + basename + ".tell.gz"
        objproducts[tellkey] = telltarget
        objdependencies[tellkey] = [extkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
        objcommands[tellkey] = TelluricWaveCommand(Dirs, telltarget, exttarget, products["WAVELENGTHPRODUCT"], products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode,False) + verstr
        ###############################################################

        ### Create TARGETS for barycentric wavelength correction ######
        rvelkey = "RVEL" + basename
        rveltarget = Dirs.PRODUCTDIR + basename + ".rvel.gz"
        objproducts[rvelkey] = rveltarget
        objdependencies[rvelkey] = ["WAVELENGTHPRODUCT"]
        objcommands[rvelkey] = HeliocentricWaveCommand(Dirs, rveltarget, products["WAVELENGTHPRODUCT"], absra_center, absdec_center, mjdate, exptime, startHA) + verstr
        ###############################################################

        ### Create TARGETS for Radial Velocity #########
        # -- create a function to obtain synthetic spectrum:
        #    queryTargetSyntheticSpectrum(absra_center,absdec_center)

        # Note that the order of input quantities below is important.
        headerinfostrg = str(expnum) + ' ' + str(mjdate) + ' ' + str(exptime) + \
        ' ' + str(airmass) + ' ' + str(outsideTemperature) + ' ' + str(windspeed) + \
        ' ' + str(relativeHumididy) + ' ' + str(atmoPressure)

        radvelkey = "RADIALVELOCITY" + basename
        radveltarget = Dirs.PRODUCTDIR + basename + ".rv"
        objproducts[radvelkey] = radveltarget
        objdependencies[radvelkey] = [extkey,tellkey,"FLATFLUXCALIBRATIONSPECTRUM",rvelkey]
        objcommands[radvelkey] = RadialVelocityCommand(Dirs, radveltarget, exttarget, telltarget, rveltarget, products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, headerinfostrg, False) + verstr
        ###############################################################

        ### Create TARGETS for calibrated spectrum *.spc ##############
        spckey = "OPSPC" + basename
        spctarget = Dirs.PRODUCTDIR + basename + ".spc.gz"
        objproducts[spckey] = spctarget
        #objdependencies[spckey] = [extkey,tellkey,rvelkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
        objdependencies[spckey] = [extkey,tellkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
        #objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],objproducts["MASTERFLUXCALIBRATION"], rveltarget, telltarget, products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        #objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],emptystring, rveltarget, telltarget,products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],emptystring, "", telltarget,products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        
        ### Create TARGETS for flux calibration Standards ######
        standardname = ""
        
        if config.hasStandard(objectname.replace(" ", "")) :
            standardname = objectname.replace(" ", "")
        elif config.hasStandard(objectname) :
            standardname = objectname

        if standardname :
            stdfcalkey = "STDFCAL" + basename
            stdfcaltarget = Dirs.PRODUCTDIR + basename + ".fcal.gz"
            
            listofstdfcalkey.append(stdfcalkey)
            listofstdfcaltargets.append(stdfcaltarget)
            
            objproducts[stdfcalkey] = stdfcaltarget
            objdependencies[stdfcalkey] = [extkey,"APERTUREPRODUCT","WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
            objcommands[stdfcalkey] = CreateFcalCommand(Dirs, stdfcaltarget, exttarget, config.getStandardDataFile(standardname), products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, products["APERTUREPRODUCT"], products["WAVELENGTHPRODUCT"], exptime) + verstr
            
            # -- The target below produces a LE-style flat response for each observation of
            #    standard stars. However it could be limited to the Moon spectrum only.
            #    E. Martioli -- Feb 22 2015
            flatrespkey = "FLATRESPONSE" + basename
            flatresptarget = Dirs.PRODUCTDIR + basename + "_flat_resp.fits.gz"
            objproducts[flatrespkey] = flatresptarget
            objdependencies[flatrespkey] = [extkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
            objcommands[flatrespkey] = CreateFlatResponseCommand(Dirs, flatresptarget, exttarget, file, config.getStandardDataFile(standardname), products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, products["APERTUREPRODUCT"], products["WAVELENGTHPRODUCT"]) + verstr

        ###############################################################

        ### Create TARGETS for LE formats #############################
        lespcnwkey = "LESPCNW" + basename
        lespcnwtarget = Dirs.PRODUCTDIR + basename + ".inw.s.gz"
        objproducts[lespcnwkey] = lespcnwtarget
        objdependencies[lespcnwkey] = [spckey]
        objcommands[lespcnwkey] = GenLEFormatsCommand(Dirs, config, lespcnwtarget, spctarget, Instmode.LESPECTRUMTYPE, objectname, 2, 4) + verstr
        #
        lespcukey = "LESPCU" + basename
        lespcutarget = Dirs.PRODUCTDIR + basename + ".iu.s.gz"
        objproducts[lespcukey] = lespcutarget
        objdependencies[lespcukey] = [spckey]
        objcommands[lespcukey] = GenLEFormatsCommand(Dirs, config, lespcutarget, spctarget, Instmode.LESPECTRUMTYPE, objectname, 3, 3) + verstr
        #
        lespcuwkey = "LESPCUW" + basename
        lespcuwtarget = Dirs.PRODUCTDIR + basename + ".iuw.s.gz"
        objproducts[lespcuwkey] = lespcuwtarget
        objdependencies[lespcuwkey] = [spckey]
        objcommands[lespcuwkey] = GenLEFormatsCommand(Dirs, config, lespcuwtarget, spctarget, Instmode.LESPECTRUMTYPE, objectname, 3, 4) + verstr
        #
        lespcnkey = "LESPCN" + basename
        lespcntarget = Dirs.PRODUCTDIR + basename + ".in.s.gz"
        objproducts[lespcnkey] = lespcntarget
        objdependencies[lespcnkey] = [spckey]
        objcommands[lespcnkey] = GenLEFormatsCommand(Dirs, config, lespcntarget, spctarget, Instmode.LESPECTRUMTYPE, objectname, 2, 3) + verstr
        #

        fitsLEproductkey = "LEFITSPRODUCT" + basename
        fitsLEproducttarget = Dirs.PRODUCTDIR + basename + "i.fits.gz"
        objproducts[fitsLEproductkey] = fitsLEproducttarget
        objdependencies[fitsLEproductkey] = [lespcnwkey,lespcukey,lespcuwkey,lespcnkey,rvelkey]
        objcommands[fitsLEproductkey] = CreateLEFITSProductCommand(Dirs, fitsLEproducttarget, lespcutarget, lespcntarget, lespcuwtarget, lespcnwtarget, file, rveltarget) + verstr

        ###############################################################

# E. Martioli Mar 16 2015 - I have commented the lines below to avoid using the flux calibration
# generated by opera. Instead, we're using a default flat response.
### Create DEPENDENCIES and COMMANDLINE for master flux calibration ######
#    if listofstdfcalkey :
#    listofstdfcalkey.append("WAVELENGTHPRODUCT")
#    listofstdfcalkey.append("FLATEXTRACTEDSPECTRUM")
#    objdependencies[masterfcalkey] = listofstdfcalkey
#    objcommands[masterfcalkey] = MasterFcalCommand(Dirs, masterfcaltarget, listofstdfcaltargets, products["FLATEXTRACTEDSPECTRUM"], products["WAVELENGTHPRODUCT"]) + verstr
#   else :
#       defaultMasterFCal = [DefaultCal.DEFAULTMASTERFLUXCAL]
#       objdependencies[masterfcalkey] = ["WAVELENGTHPRODUCT","FLATEXTRACTEDSPECTRUM"]
#       objcommands[masterfcalkey] = MasterFcalCommand(Dirs, masterfcaltarget, defaultMasterFCal, products["FLATEXTRACTEDSPECTRUM"], products["WAVELENGTHPRODUCT"]) + verstr
##########################################################################

    return objproducts, objdependencies, objcommands
#-------------------------------------------

######### PRODUCT FILE NAMES ###########
def setProductFilenames(Dirs,night,Instmode,Readmode,DefaultCal,Keywords,allowanyreadout) :
    
    INSTCONFIGPREFIX = Dirs.PRODUCTDIR + night + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
    
    products = {}
    
    biascheck = testBiasListCommand(Dirs, Instmode, Readmode, Keywords)
    if subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["BIASLIST"] = INSTCONFIGPREFIX + "_bias.list"
        products["MASTERBIAS"] = INSTCONFIGPREFIX + "_masterbias.fits.gz"
    else :
        products["MASTERBIAS"] = DefaultCal.DEFAULTMASTERBIAS
    
    flatcheck = testCalibrationListCommand(Dirs, Instmode, Readmode, Keywords.FLATKEYWORD, allowanyreadout)
    if subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["FLATLIST"] = INSTCONFIGPREFIX + "_flat.list"
        products["MASTERFLAT"] = INSTCONFIGPREFIX + "_masterflat.fits.gz"
    else :
        products["MASTERFLAT"] = DefaultCal.DEFAULTMASTERFLAT
    
    compcheck = testCalibrationListCommand(Dirs, Instmode, Readmode, Keywords.COMPKEYWORD, allowanyreadout)
    if subprocess.check_output(compcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["COMPLIST"] = INSTCONFIGPREFIX + "_comp.list"
        products["MASTERCOMP"] = INSTCONFIGPREFIX + "_mastercomp.fits.gz"
    else :
        products["MASTERCOMP"] = DefaultCal.DEFAULTMASTERCOMP

    if subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') and \
        subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["GAINPRODUCT"] = INSTCONFIGPREFIX + ".gain.gz"
    else :
        products["GAINPRODUCT"] = DefaultCal.DEFAULTGAIN

    products["ORDERSPACINGPRODUCT"] = INSTCONFIGPREFIX + ".ordp.gz"
    products["GEOMETRYPRODUCT"] = INSTCONFIGPREFIX + ".geom.gz"
    products["INSTRUMENTPROFILEPRODUCT"] = INSTCONFIGPREFIX + ".prof.gz"
    products["APERTUREPRODUCT"] = INSTCONFIGPREFIX + ".aper.gz"
    products["COMPEXTRACTEDSPECTRUM"] = INSTCONFIGPREFIX + "_comp.e.gz"
    products["FLATEXTRACTEDSPECTRUM"] = INSTCONFIGPREFIX + "_flat.e.gz"
    products["FIRSTWAVELENGTHPRODUCT"] = INSTCONFIGPREFIX + ".wcar.gz"
    products["WAVELENGTHPRODUCT"] = INSTCONFIGPREFIX + ".wcal.gz"
    products["FLATFLUXCALIBRATIONSPECTRUM"] = INSTCONFIGPREFIX + "_flat.fcal.gz"
    #products["OBJECTLIST"] = INSTCONFIGPREFIX + "_object.list"
    
    return products
##################################

######### PRODUCT FILE NAMES ###########
def setPlotFilenames(Dirs,night,Instmode,Readmode, plotbool) :
    
    INSTCONFIGPREFIX = Dirs.PRODUCTDIR + night + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
    
    plots = {}
    
    if plotbool :
        plots["ORDERSPACINGPRODUCT"] = {'ORDSPCPLOTFILE': INSTCONFIGPREFIX + "_spcplot.eps", \
            'ORDSPCDATAFILE': INSTCONFIGPREFIX + "_spcplot.dat", \
                'ORDSPCSCRIPTFILE': INSTCONFIGPREFIX + "_spcplot.gnu"}
        
        plots["GEOMETRYPRODUCT"] =  {'GEOMPLOTFILE': INSTCONFIGPREFIX + "_geomplot.eps", \
            'GEOMDATAFILE': INSTCONFIGPREFIX + "_geomplot.dat", \
                'GEOMSCRIPTFILE': INSTCONFIGPREFIX + "_geomplot.gnu"}
    
        plots["INSTRUMENTPROFILEPRODUCT"] =  {'PROFPLOTFILE': INSTCONFIGPREFIX + "_profplot.eps", \
            'PROFDATAFILE': INSTCONFIGPREFIX + "_profplot.dat", \
                'PROFSCRIPTFILE': INSTCONFIGPREFIX + "_profplot.gnu"}

        plots["APERTUREPRODUCT"] =  {'APERPLOTFILE': INSTCONFIGPREFIX + "_aperplot.eps", \
            'APERDATAFILE': INSTCONFIGPREFIX + "_aperplot.dat", \
            'APERSCRIPTFILE': INSTCONFIGPREFIX + "_aperplot.gnu", \
            'APERTILTPLOTFILE': INSTCONFIGPREFIX + "_tiltplot.eps", \
            'APERTILTDATA1FILE': INSTCONFIGPREFIX + "_tiltplot1.dat", \
            'APERTILTDATA2FILE': INSTCONFIGPREFIX + "_tiltplot2.dat", \
            'APERTILTSCRIPTFILE': INSTCONFIGPREFIX + "_tiltplot.gnu"}
    
        plots["FIRSTWAVELENGTHPRODUCT"] =  {'WAVEORDSPLOTFILE': INSTCONFIGPREFIX + "_waveordsplot.eps", \
            'WAVESPECPLOTFILE': INSTCONFIGPREFIX + "_wavespecplot.eps", \
            'WAVESPECSCRIPTFILE': INSTCONFIGPREFIX + "_wavespecplot.gnu", \
            'WAVEORDSCRIPTFILE': INSTCONFIGPREFIX + "_waveordplot.gnu", \
            'WAVEORDSDATAFILE': INSTCONFIGPREFIX + "_waveordsplot.dat", \
            'WAVEATLASDATAFILE': INSTCONFIGPREFIX + "_waveatlasplot.dat", \
            'WAVECOMPDATAFILE': INSTCONFIGPREFIX + "_wavecompplot.dat", \
            'WAVELINESDATAFILE': INSTCONFIGPREFIX + "_wavelinesplot.dat"}

    else :
        plots["ORDERSPACINGPRODUCT"] = {'ORDSPCPLOTFILE': "",'ORDSPCDATAFILE': "",'ORDSPCSCRIPTFILE': ""}
        plots["GEOMETRYPRODUCT"] = {'GEOMPLOTFILE': "",'GEOMDATAFILE': "",'GEOMSCRIPTFILE': ""}
        plots["INSTRUMENTPROFILEPRODUCT"] = {'PROFPLOTFILE': "",'PROFDATAFILE': "",'PROFSCRIPTFILE': ""}
        plots["APERTUREPRODUCT"] = {'APERPLOTFILE': "",'APERDATAFILE': "",'APERSCRIPTFILE': "",\
            'APERTILTPLOTFILE': "",'APERTILTDATA1FILE': "",'APERTILTDATA2FILE': "",'APERTILTSCRIPTFILE': ""}
        plots["FIRSTWAVELENGTHPRODUCT"] = {'WAVEORDSPLOTFILE': "",'WAVESPECPLOTFILE': "",'WAVESPECSCRIPTFILE': "", \
            'WAVEORDSCRIPTFILE': "",'WAVEORDSDATAFILE': "",'WAVEATLASDATAFILE': "", \
            'WAVECOMPDATAFILE': "", 'WAVELINESDATAFILE': ""}

    return plots
##################################


######### DEPENDENCIES ###########
def setDependencies(products,Instmode) :
    dependencies = {}
    
    # Template for a target is: targets[key] = filename
    if "BIASLIST" in products :
        dependencies["BIASLIST"] = []
        dependencies["MASTERBIAS"] = ["BIASLIST"]
    else :
        dependencies["MASTERBIAS"] = []
    
    if "FLATLIST" in products :
        dependencies["FLATLIST"] = []
        dependencies["MASTERFLAT"] = ["FLATLIST"]
    else :
        dependencies["MASTERFLAT"] = []
    
    if "COMPLIST" in products :
        dependencies["MASTERCOMP"] = ["COMPLIST","MASTERBIAS"]
        dependencies["COMPLIST"] = []
    else :
        dependencies["MASTERCOMP"] = ["MASTERBIAS"]
    
    if "FLATLIST" in products and "BIASLIST" in products :
        dependencies["GAINPRODUCT"] = ["BIASLIST","FLATLIST"]
    else :
        dependencies["GAINPRODUCT"] = []
    
    dependencies["ORDERSPACINGPRODUCT"] = ["GAINPRODUCT","MASTERBIAS","MASTERFLAT"]
    dependencies["GEOMETRYPRODUCT"] = ["GAINPRODUCT","MASTERBIAS","MASTERFLAT","ORDERSPACINGPRODUCT"]
    
    if (Instmode.mode == 1 or Instmode.mode == 3) :
        dependencies["INSTRUMENTPROFILEPRODUCT"] = ["GEOMETRYPRODUCT","GAINPRODUCT","MASTERBIAS","MASTERFLAT"]
    elif (Instmode.mode == 2) :
        dependencies["INSTRUMENTPROFILEPRODUCT"] = ["GEOMETRYPRODUCT","GAINPRODUCT","MASTERBIAS","MASTERFLAT","MASTERCOMP"]

    dependencies["APERTUREPRODUCT"] = ["GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","ORDERSPACINGPRODUCT"]
    dependencies["COMPEXTRACTEDSPECTRUM"] = ["MASTERCOMP","MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
    dependencies["FLATEXTRACTEDSPECTRUM"] = ["MASTERFLAT","MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
    dependencies["FIRSTWAVELENGTHPRODUCT"] = ["GEOMETRYPRODUCT","COMPEXTRACTEDSPECTRUM"]
    dependencies["WAVELENGTHPRODUCT"] = ["FIRSTWAVELENGTHPRODUCT","COMPEXTRACTEDSPECTRUM"]
    dependencies["FLATFLUXCALIBRATIONSPECTRUM"] = ["FLATEXTRACTEDSPECTRUM","WAVELENGTHPRODUCT"]
    
#dependencies["OBJECTLIST"]=[]

    return dependencies
##################################

######### COMMAND LINES ###########
def setPipelineCommands(products,Dirs,night,Instmode,Readmode,keywords,config,plots,allowanyreadout,verbose) :
    commands = {}
    
    verstr = ""
    if (verbose) :
        verstr = " --verbose"
    
    command = ""
    
    if "BIASLIST" in products :
        commands["BIASLIST"] = BiasListCommand(Dirs, Instmode, Readmode, keywords, products["BIASLIST"])
        commands["MASTERBIAS"] = MasterCalibrationCommand(Dirs, "operaMasterBias", products["MASTERBIAS"], products["BIASLIST"]) + verstr
    
    if "FLATLIST" in products :
        commands["FLATLIST"] = CalibrationListCommand(Dirs, Instmode, Readmode, keywords.FLATKEYWORD, products["FLATLIST"],allowanyreadout)
        commands["MASTERFLAT"] = MasterCalibrationCommand(Dirs, "operaMasterFlat", products["MASTERFLAT"], products["FLATLIST"]) + verstr
    
    if "COMPLIST" in products :
        commands["COMPLIST"] = CalibrationListCommand(Dirs, Instmode, Readmode, keywords.COMPKEYWORD, products["COMPLIST"],allowanyreadout)
        commands["MASTERCOMP"] = MasterComparisonCommand(Dirs, products["MASTERCOMP"], products["COMPLIST"], config.BADPIXELMASK, products["MASTERBIAS"]) + verstr

    if "FLATLIST" in products and "BIASLIST" in products :
        commands["GAINPRODUCT"] = GainCommand(Dirs,products["GAINPRODUCT"], products["BIASLIST"], products["FLATLIST"], config.BADPIXELMASK, Readmode) + verstr
    
    commands["ORDERSPACINGPRODUCT"] = OrderSpacingCommand(Dirs,products["ORDERSPACINGPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,Instmode,plots["ORDERSPACINGPRODUCT"]) + verstr
    commands["GEOMETRYPRODUCT"] = GeometryCommand(Dirs,products["GEOMETRYPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["ORDERSPACINGPRODUCT"],Instmode,plots["GEOMETRYPRODUCT"]) + verstr
    commands["INSTRUMENTPROFILEPRODUCT"] = InstrumentProfileCommand(Dirs,products["INSTRUMENTPROFILEPRODUCT"], products["GEOMETRYPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],products["MASTERCOMP"],"",2,config.BADPIXELMASK,Instmode,plots["INSTRUMENTPROFILEPRODUCT"]) + verstr
    commands["APERTUREPRODUCT"] = ApertureCommand(Dirs,products["APERTUREPRODUCT"],products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["ORDERSPACINGPRODUCT"],Instmode,plots["APERTUREPRODUCT"]) + verstr
    commands["COMPEXTRACTEDSPECTRUM"] = compRawExtractionCommand(Dirs,products["COMPEXTRACTEDSPECTRUM"],products["MASTERCOMP"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"]) + verstr
    commands["FLATEXTRACTEDSPECTRUM"] = calibrationExtractionCommand(Dirs,products["FLATEXTRACTEDSPECTRUM"],products["MASTERFLAT"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"]) + verstr
    commands["FIRSTWAVELENGTHPRODUCT"] = WavelengthCommand(Dirs, products["FIRSTWAVELENGTHPRODUCT"], products["GEOMETRYPRODUCT"], products["COMPEXTRACTEDSPECTRUM"], Instmode, config ,plots["FIRSTWAVELENGTHPRODUCT"]) + verstr
    commands["WAVELENGTHPRODUCT"] = StitchOrdersCommand(Dirs, products["WAVELENGTHPRODUCT"], products["COMPEXTRACTEDSPECTRUM"], products["FIRSTWAVELENGTHPRODUCT"]) + verstr
    commands["FLATFLUXCALIBRATIONSPECTRUM"] = FlatFluxCalibrationCommand(Dirs, products["FLATFLUXCALIBRATIONSPECTRUM"],products["FLATEXTRACTEDSPECTRUM"], Instmode, products["WAVELENGTHPRODUCT"]) + verstr

#commands["OBJECTLIST"] = ObjectListCommandToFile(Dirs, Instmode, Readmode, keywords.OBJECTKEYWORD, products["OBJECTLIST"])

    return commands
##################################

class Products :
    'Common base class for pipeline products'
    targets = {}
    dependencies = {}
    commands = {}
    
    plots = {}
    existenceStatus = {}
    
    trace = False
    simulation = False
    
    def __init__(self,Targets, Plots, Dependencies, Commands, Trace):
        self.targets = Targets
        self.plots = Plots
        self.dependencies = Dependencies
        self.commands = Commands
        self.setExistenceStatus()
        self.trace = Trace
    
    # Function below sets processing in simulation mode
    def setSimulation(self) :
        self.simulation = True
    #-------------------------------------------
    
    # Function below resets simulation back to non-simulation mode
    def resetSimulation(self) :
        self.simulation = False
        self.setExistenceStatus()
    #-------------------------------------------
    
    # Function below returns the product file name <string> associated to a target
    def getTarget(self, key) :
        return self.targets[key]
    #-------------------------------------------
    
    # Function to display targets
    def displayTargets(self) :
        for targetItem in self.targets.items():
            print '"',targetItem[0],'" :', targetItem[1]
    #-------------------------------------------
    
    # Function to add targets from input dictionary.
    # Input is a vector of three dictionaries: target{}, dependencies{}, commands{}
    def addTargets(self, inputDictVec) :
        # input: inputDictVec = [target{}, dependencies{}, commands{}]
        
        self.targets.update(inputDictVec[0])
        
        newkeys = []
        for item in inputDictVec[0].items():
            newkeys.append(item[0])
        
        dictbool = dict.fromkeys(newkeys, False)
        self.existenceStatus.update(dictbool)
        
        self.dependencies.update(inputDictVec[1])
        self.commands.update(inputDictVec[2])
        
        self.setExistenceStatus()
    #-------------------------------------------

    # Function below returns the list of dependencies to produce product file
    def getDependencies(self, key) :
        return self.dependencies[key]
    #-------------------------------------------
    
    # Function below returns the command line to produce product file
    def getCommandLine(self, key) :
        return self.commands[key]
    #-------------------------------------------
    
    # Function below executes a single target
    def executeTarget(self, key) :
        if (self.simulation == False) :
            self.setExistenceStatus()
        # If target exists then set target existenceStatus to True and exit
        if self.existenceStatus[key] == True :
            return True
        # Otherwise execute command and recheck target existence and exit
        else :
            if self.resolveDependencies(key) == True :
                if (self.simulation == True) :
                    if key in self.commands :
                        print self.getCommandLine(key)
                        self.existenceStatus[key] = True
                else :
                    try :
                        if(self.trace == True) :
                            print self.getCommandLine(key)
                        
                        if key in self.commands :
                            #subprocess.check_output(self.getCommandLine(key),stderr=subprocess.STDOUT,shell=True)
                            os.system(self.getCommandLine(key))
                        
                        if os.path.exists(self.targets[key]) :
                            self.existenceStatus[key] = True
                    except :
                        print "Error: can\'t execute command: ",self.getCommandLine(key)

        return self.getExistenceStatus(key)
    #-------------------------------------------
    # Function below executes all targets matching any substring given in the input list of substrings
    def executeTargetsWithSubstrings(self,substrings) :
        self.setExistenceStatus()
        targetItems = self.targets.items()
        for targetItem in targetItems :
            for str in substrings :
                if (str in targetItem[0]) or (str in targetItem[1]) :
                    self.executeTarget(targetItem[0])
    
        if (self.simulation == True) :
            self.setExistenceStatus()

    # Function below returns the status of a target <bool>
    def getExistenceStatus(self, key) :
        return self.existenceStatus[key]
    #-------------------------------------------
    
    # Function below sets the status of each target based on existence check
    def setExistenceStatus(self) :
        # Convert to list of tuples
        targetItems = self.targets.items()
        for targetItem in targetItems:
            if os.path.exists(targetItem[1]) :
                self.existenceStatus[targetItem[0]] = True
            else :
                self.existenceStatus[targetItem[0]] = False
    #-------------------------------------------

    # Function below execute every step where ExistenceStatus is FALSE
    def executeAllTargets(self) :
        self.setExistenceStatus()
        targetItems = self.targets.items()
        for targetItem in targetItems:
            if (self.existenceStatus[targetItem[0]] == False) :
                self.executeTarget(targetItem[0])
        if (self.simulation == True) :
            self.setExistenceStatus()
    #-------------------------------------------

    # Function to check dependencies.
    # Return False if any dependency is missing.
    def checkDependencies(self, key) :
        dependenciesstatus = True
        for dkey in self.dependencies[key] :
            if self.existenceStatus[dkey] == False :
                dependenciesstatus = False
        return dependenciesstatus
    #-------------------------------------------

    # Function to resolve dependencies
    def resolveDependencies(self, key) :
        if (self.checkDependencies(key) == False) :
            for dkey in self.dependencies[key] :
                if self.existenceStatus[dkey] == False :
                    self.executeTarget(dkey)

        return self.checkDependencies(key)
    #-------------------------------------------

    # Function to delete all targets
    def cleanTargets(self) :
        targetItems = self.targets.items()
        for targetItem in targetItems:
            if os.path.exists(targetItem[1]) :
                try :
                    if(self.simulation == True) :
                        print "rm " + targetItem[1]
                    else :
                        os.remove(targetItem[1])
                except:
                    print "Error: can\'t remove file: ", targetItem[1]
                else:
                    self.existenceStatus[targetItem[0]] = False
    
        return
    #-------------------------------------------

    # Function to delete all plots
    def cleanPlots(self) :
        for plotItem in self.plots.items():
            for plotItemItem in plotItem[1].items() :
                if os.path.exists(plotItemItem[1]) :
                    try :
                        if(self.simulation == True) :
                            print "rm " + plotItemItem[1]
                        else :
                            os.remove(plotItemItem[1])
                    except:
                        print "Error: can\'t remove file: ", plotItemItem[1]
        return
    #-------------------------------------------

    # Function to delete all targets and plots
    def cleanAll(self) :
        self.cleanPlots()
        self.cleanTargets()
        return
    #-------------------------------------------
    
    # Function to remove targets matching any substring given in the input list of substrings
    def removeTargets(self,substrings) :
        targetItems = self.targets.items()
        for targetItem in targetItems :
            for str in substrings :
                if (str in targetItem[0]) or (str in targetItem[1]) :
                    self.targets.pop(targetItem[0], None)
        return
    #-------------------------------------------

##################################

#### Function to generate a command line for BIAS list: ####
def BiasListCommand(Dirs, Instmode, Readmode, keywords,listfilename) :
    commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
    ' -q "INSTRUME EREADSPD OBSTYPE"' + \
    ' INSTRUME="'+ Instmode.INSTRUME +'" EREADSPD="'+Readmode.READOUTSPEED+'" OBSTYPE='+keywords.BIASKEYWORD + \
    ' > ' + listfilename
    return commandline
###########################################

#### Function to test BIAS list: ####
def testBiasListCommand(Dirs, Instmode, Readmode, keywords) :
    commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
    ' -q "INSTRUME EREADSPD OBSTYPE"' + \
    ' INSTRUME="'+ Instmode.INSTRUME + '" EREADSPD="'+Readmode.READOUTSPEED+'" OBSTYPE='+keywords.BIASKEYWORD
    return commandline
###########################################

#### Function to generate a command line for FLAT, COMP lists: ####
def CalibrationListCommand(Dirs, Instmode, Readmode, keyword, listfilename, allowanyreadout) :
    commandline = ""
    if allowanyreadout :
        commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
        ' -q "INSTRUME GSLICER OBSTYPE"' + \
        ' INSTRUME="'+ Instmode.INSTRUME +'" GSLICER="'+Instmode.INSTRUMENTMODEKEY+'" OBSTYPE='+keyword + \
        ' > ' + listfilename
    else :
        commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
        ' -q "INSTRUME GSLICER EREADSPD OBSTYPE"' + \
        ' INSTRUME="'+ Instmode.INSTRUME +'" GSLICER="'+Instmode.INSTRUMENTMODEKEY+'" EREADSPD="'+Readmode.READOUTSPEED+'" OBSTYPE='+keyword + \
        ' > ' + listfilename
    return commandline
###########################################

#### Function to test FLAT, COMP lists: ####
def testCalibrationListCommand(Dirs, Instmode, Readmode, keyword, allowanyreadout) :
    commandline = ""
    if allowanyreadout :
        commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
        ' -q "INSTRUME GSLICER OBSTYPE"' + \
        ' INSTRUME="'+ Instmode.INSTRUME +'" GSLICER="'+Instmode.INSTRUMENTMODEKEY+'" OBSTYPE='+keyword
    else :
        commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
        ' -q "INSTRUME GSLICER EREADSPD OBSTYPE"' + \
        ' INSTRUME="'+ Instmode.INSTRUME +'" GSLICER="'+Instmode.INSTRUMENTMODEKEY+'" EREADSPD="'+Readmode.READOUTSPEED+'" OBSTYPE='+keyword
    return commandline
###########################################

#### Function to generate a command line for mastercalibrations: ####
def MasterCalibrationCommand(Dirs, bin, product, list) :
    commandline = Dirs.EXE + bin + " --output=" + product + " --imagelistfile=" + list
    return commandline
###########################################

#### Function to generate a command line for mastercomparison: ####
def MasterComparisonCommand(Dirs, product, list, badpix, masterbias) :
    commandline = Dirs.EXE + 'operaMasterComparison --output=' + product + \
    ' --imagelistfile=' + list + ' --badpixelmask=' + badpix + ' --masterbias=' + masterbias + \
    ' --combineMethod=1 --saturationLimit=65535 --outputExposureTime=40 ' +\
    ' --biasConstant=0 --truncateOuputFluxToSaturation=1 --expTimeFITSKeyword=EXPTIME'
    return commandline
###########################################

#### Function to generate a command line for operaGain: ####
def GainCommand(Dirs, product, biaslist, flatlist, badpix, Readmode) :
    commandline = Dirs.EXE + 'operaGain --output=' + product + \
    ' --biaslistfile=' + biaslist + ' --flatlistfile=' + flatlist + \
    ' --badpixelmask=' + badpix + ' --defaultgain=' +  Readmode.DEFAULTGAIN + ' --defaultnoise=' + Readmode.DEFAULTNOISE + \
    ' --DATASEC="1 2048 1 4608" --numberofamplifiers=1 --subwindow="100 800 500 4000" --gainMinPixPerBin=1000 --gainMaxNBins=100 --gainLowestCount=1000 --gainHighestCount=30000'
    return commandline
###########################################

#### Function to generate a command line for operaOrderSpacingCalibration: ####
def OrderSpacingCommand(Dirs, product, gainproduct, masterbias, masterflat, badpix, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["ORDSPCPLOTFILE"] + ' --datafilename=' + plots["ORDSPCDATAFILE"] + ' --scriptfilename=' + plots["ORDSPCSCRIPTFILE"]
    
    commandline = Dirs.EXE + 'operaOrderSpacingCalibration --orderspacingoutput=' + product + \
    ' --inputGainFile=' + gainproduct + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + \
    ' --badpixelmask=' + badpix + ' --aperture=' + str(Instmode.ORDSPCAPERTURE) + ' --referenceOrderNumber=' + str(Instmode.SPACINGREFERENCEORDERNUMBER) + \
    ' --referenceOrderSeparation=' + str(Instmode.SPACINGREFERENCEORDERSEPARATION) + \
    ' --numberOfsamples=30 --sampleCenterPosition=2300 --subformat="8 2040 3 4600"' + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaGeometryCalibration: ####
def GeometryCommand(Dirs, product, gainproduct, masterbias, masterflat, badpix, orderspacing, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["GEOMPLOTFILE"] + ' --datafilename=' + plots["GEOMDATAFILE"] + ' --scriptfilename=' + plots["GEOMSCRIPTFILE"]
    
    commandline = Dirs.EXE + 'operaGeometryCalibration --outputGeomFile=' + product + \
    ' --inputGainFile=' + gainproduct + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + \
    ' --badpixelmask=' + badpix + ' --inputOrderSpacing=' + orderspacing + ' --aperture=' + str(Instmode.GEOMAPERTURE) + \
    ' --maxorders=' + str(Instmode.GEOMMAXNORDERS) + ' --minordertouse=' + str(Instmode.GEOMMINORDERTOUSE) + \
    ' --recenterIPUsingSliceSymmetry=' + str(Instmode.recenterIPUsingSliceSymmetry) + ' --totalNumberOfSlices=' + str(Instmode.NUMBEROFSLICES) + \
    ' --subformat="8 2040 3 4600" --detectionMethod=2 --FFTfilter=0 --nsamples=5  --orderOfTracingPolynomial=3' + \
    ' --binsize=25 --colDispersion=1 --invertOrders=1 --referenceOrderSamplePosition=2300 --graces=1' + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaInstrumentProfileCalibration: ####
def InstrumentProfileCommand(Dirs, product, geomproduct, gainproduct, masterbias, masterflat, mastercomp, masteralign, ipmethod, badpix, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["PROFPLOTFILE"] + ' --datafilename=' + plots["PROFDATAFILE"] + ' --scriptfilename=' + plots["PROFSCRIPTFILE"]
    
    commandline = Dirs.EXE + "operaInstrumentProfileCalibration --outputProf=" + product + \
    ' --geometryfilename=' + geomproduct + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + \
    ' --badpixelmask=' + badpix + ' --mastercomparison=' +  mastercomp + ' --masterfabperot=' + masteralign + ' --gainfilename=' + gainproduct + \
    ' --xSize=' + str(Instmode.IPXSIZE) + ' --ySize=' + str(Instmode.IPYSIZE) + ' --xSampling=' + str(Instmode.IPXSAMPLING) + \
    ' --ySampling=' + str(Instmode.IPYSAMPLING) + ' --referenceLineWidth=' + str(Instmode.REFERENCELINEWIDTH) + \
    ' --binsize=100 --ordernumber=-999 --method=' + str(ipmethod) + ' --tilt=' + str(Instmode.TILTANGLE) + \
    ' --spectralElementHeight=1.0 --maxthreads=4 --minimumlines=5' + \
    ' --LocalMaxFilterWidth=3.0 --DetectionThreshold=0.2 --MinPeakDepth=1.5' + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaExtractionApertureCalibration: ####
def ApertureCommand(Dirs, product, geomproduct, profproduct, orderspacing, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["APERPLOTFILE"] + \
        ' --datafilename=' + plots["APERDATAFILE"] + \
            ' --scriptfilename=' + plots["APERSCRIPTFILE"] + \
                ' --tiltplotfilename=' + plots["APERTILTPLOTFILE"] + \
                    ' --tiltdata1filename=' + plots["APERTILTDATA1FILE"] + \
                        ' --tiltdata2filename=' + plots["APERTILTDATA2FILE"] + \
                            ' --tiltscriptfilename=' + plots["APERTILTSCRIPTFILE"]

    commandline = Dirs.EXE + "operaExtractionApertureCalibration --outputApertureFile=" + product + \
    ' --inputgeom=' + geomproduct + ' --inputprof=' + profproduct + ' --inputorderspacing=' + orderspacing + \
    ' --numberOfBeams=' + str(Instmode.NUMBEROFBEAMS) + ' --gapBetweenBeams=' + str(Instmode.APERGAP) + \
    ' --apertureHeight=' + str(Instmode.APERTUREHEIGHT) + ' --apertureWidth=' + str(Instmode.APERAPERTURE) + \
    ' --constantTilt=' + str(Instmode.CONSTANTTILTFLAG) + \
    ' --backgroundAperture=1.0 --pickImageRow=0 --nRowSamples=10 --xbin=10' + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for Raw Extraction of Comparison spectra: ####
def compRawExtractionCommand(Dirs, product, inputImage, masterbias, masterflat, badpix, gainproduct, geomproduct, profproduct, aperproduct) :
    commandline = Dirs.EXE + 'operaExtraction --outputSpectraFile=' + product + \
    ' --inputImage=' + inputImage + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + ' --badpixelmask=' + badpix + \
    ' --inputGainFile=' + gainproduct + ' --inputGeometryFile=' + geomproduct + \
    ' --inputInstrumentProfileFile=' + profproduct + ' --inputApertureFile=' + aperproduct + \
    ' --spectrumtype=5 --spectrumtypename=RawBeamSpectrum --starplusskymode=0  --maxthreads=4'
    
    return commandline
###########################################

#### Function to generate a command line for Optimal Extraction of spectra: ####
def calibrationExtractionCommand(Dirs, product, inputImage, masterbias, masterflat, badpix, gainproduct, geomproduct, profproduct, aperproduct) :
    commandline = Dirs.EXE + 'operaExtraction --outputSpectraFile=' + product + \
    ' --inputImage=' + inputImage + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + ' --badpixelmask=' + badpix + \
    ' --inputGainFile=' + gainproduct + ' --inputGeometryFile=' + geomproduct + \
    ' --inputInstrumentProfileFile=' + profproduct + ' --inputApertureFile=' + aperproduct + \
    ' --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --backgroundBinsize=300 --sigmaclip=6 ' + \
    ' --removeBackground=0 --iterations=3 --onTargetProfile=1 --usePolynomialFit=0 --starplusskymode=0 --maxthreads=4'
    
    return commandline
###########################################

#### Function to generate a command line for operaWavelengthCalibration: ####
def WavelengthCommand(Dirs, product, geomproduct, compspectrum, Instmode, config, plots) :
    
    plotstring = ' --ordersplotfilename=' + plots["WAVEORDSPLOTFILE"] + ' --specplotfilename=' + plots["WAVESPECPLOTFILE"] + \
    ' --ordersscriptfilename=' + plots["WAVEORDSCRIPTFILE"] + ' --specscriptfilename=' + plots["WAVESPECSCRIPTFILE"] + \
    ' --ordersdatafilename=' + plots["WAVEORDSDATAFILE"] + ' --atlasdatafilename=' + plots["WAVEATLASDATAFILE"] + \
    ' --compdatafilename=' + plots["WAVECOMPDATAFILE"] + ' --linesdatafilename=' + plots["WAVELINESDATAFILE"]
    
    commandline = Dirs.EXE + 'operaWavelengthCalibration --outputWaveFile=' + product + \
    ' --inputGeomFile=' + geomproduct + ' --uncalibrated_spectrum=' + compspectrum +  \
    ' --inputWaveFile=' + config.WAVEFIRSTGUESS + ' --uncalibrated_linewidth=' + str(Instmode.WAVEUNCALLINEWIDTH) + \
    ' --atlas_lines=' + config.THARATLASLINES + ' --atlas_spectrum=' + config.THARATLASSPECTRUM + \
    ' --parseSolution=0 --ParRangeSizeInPerCent=1.0 --NpointsPerPar=3000 --maxNIter=40 --minNumberOfLines=40' +\
    ' --maxorderofpolynomial=4 --dampingFactor=0.85 --initialAcceptableMismatch=1.5 --nsigclip=2.25 ' +\
    ' --normalizeUncalibratedSpectrum=0 --normalizationBinSize=180 --LocalMaxFilterWidth=6' +\
    ' --DetectionThreshold=0.05 --MinPeakDepth=1.0' + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaStitchOrders: ####
def StitchOrdersCommand(Dirs, product, compspectrum, wave) :
    
    commandline = Dirs.EXE + 'operaStitchOrders --outputWaveFile=' + product + \
    ' --inputSpectrum=' + compspectrum + ' --inputWaveFile=' + wave + \
    ' --orderOfReference=37 --DWavelengthRange=0.2 --DWavelengthStep=0.00005 --XCorrelationThreshold=0.1 --sigmaThreshold=2.5'
    
    return commandline
###########################################

#### Function to generate a command line for operaCreateFlatFieldFluxCalibration: ####
def FlatFluxCalibrationCommand(Dirs, product, flatspectrum, Instmode, wave) :
    
    commandline = Dirs.EXE + 'operaCreateFlatFieldFluxCalibration --outputFluxCalibrationFile=' + product + \
    ' --inputMasterFlatSpectrum=' + flatspectrum + ' --wavelengthCalibration=' + wave + \
    ' --wavelengthForNormalization=' + str(Instmode.WAVELENGTHFORNORMALIZATION) + ' --binsize=500'
    
    return commandline
###########################################


#### Function to generate a command line for OBJECT list on screen: ####
def ObjectListCommand(Dirs, Instmode, Readmode, keyword) :
    commandline = Dirs.EXE + 'operaQueryImageInfo --directory=' + Dirs.DATADIR + \
    ' -q "INSTRUME GSLICER EREADSPD OBSTYPE"' + ' INSTRUME="'+ Instmode.INSTRUME + \
    '" GSLICER="'+Instmode.INSTRUMENTMODEKEY+'" EREADSPD="'+Readmode.READOUTSPEED+'" OBSTYPE='+keyword
    
    return commandline
###########################################

#### Function to generate a command line to create OBJECT list into a file: ####
def ObjectListCommandToFile(Dirs, Instmode, Readmode, keyword, listfilename) :
    commandline = ObjectListCommand(Dirs, Instmode, Readmode, keyword) + \
    ' > ' + listfilename
    return commandline
###########################################

#### Function to generate a command line for Optimal Extraction of object spectra: ####
def objectExtractionCommand(Dirs, product, inputImage, masterbias, masterflat, badpix, gainproduct, geomproduct, profproduct, aperproduct, Instmode) :
    commandline = Dirs.EXE + 'operaExtraction --outputSpectraFile=' + product + \
    ' --inputImage=' + inputImage + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + ' --badpixelmask=' + badpix + \
    ' --inputGainFile=' + gainproduct + ' --inputGeometryFile=' + geomproduct + \
    ' --inputInstrumentProfileFile=' + profproduct + ' --inputApertureFile=' + aperproduct + \
    ' --starplusskymode=' + str(Instmode.STARPLUSKYMODEFLAG) + ' ' + Instmode.INVERTSKYFIBERFLAG + \
    ' --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --backgroundBinsize=300 --sigmaclip=6 ' + \
    ' --removeBackground=0 --iterations=3 --onTargetProfile=1 --usePolynomialFit=0 --maxthreads=4'
    
    return commandline
###########################################


#### Function to generate a command line for Telluric Wavelength Correction: ####
def TelluricWaveCommand(Dirs, product, inputSpectrum, wave, flatSpectrum, config, Instmode, plotbool) :
    
    if plotbool :
        plotstring = ' --xcorrsplotfilename=' + "xcorr_tmp.eps" + ' --specplotfilename=' + "spec_tmp.eps" + \
    ' --xcorrscriptfilename=' + "xcorr_tmp.gnu" + ' --specscriptfilename=' + "spec_tmp.gnu" + \
    ' --xcorrdatafilename=' + "xcorr_tmp.dat" + ' --xcorrfitdatafilename=' + "xcorr-fit_tmp.dat" + \
    ' --specdatafilename=' + "spec_tmp.dat"
    else :
        plotstring = ''
    
    if Instmode.STARPLUSKYMODEFLAG != 0 :
        flagstring = ' --StarPlusSky'
    else :
        flagstring = ''
    
    commandline = Dirs.EXE + 'operaTelluricWavelengthCorrection --outputWaveFile=' + product + \
    ' --inputObjectSpectrum=' + inputSpectrum + ' --inputWaveFile=' + wave + \
    ' --telluric_lines=' + config.TELLURICLINES + ' --telluric_spectrum=' + config.TELLURICSPECTRUM + \
    ' --inputWavelengthMaskForTelluric=' + config.TELLURICWAVELENGTHMASK + \
    ' --spectralResolution=' + str(Instmode.SPECTRALRESOLUTION) + \
    ' --radialVelocityRange=' + str(Instmode.RADIALVELOCITYRANGE) + \
    ' --radialVelocityStep=' + str(Instmode.RADIALVELOCITYSTEP) + \
    ' --XCorrelationThreshold=0.1 --normalizationBinsize=110' + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --useFitToFindMaximum' + \
        flagstring + plotstring

    return commandline
##########################################

#### Function to generate a command line for Telluric Wavelength Correction: ####
def RadialVelocityCommand(Dirs, product, inputSpectrum, wave, rvel, flatSpectrum, config, Instmode, headerinfostrg, plotbool) :
    
    #    if plotbool :
    plotstring = ' --xcorrsplotfilename=' + "xcorr_rv.eps" + ' --specplotfilename=' + "spec_rv.eps" + \
    ' --xcorrscriptfilename=' + "xcorr_rv.gnu" + ' --specscriptfilename=' + "spec_rv.gnu" + \
    ' --xcorrdatafilename=' + "xcorr_rv.dat" + ' --xcorrfitdatafilename=' + "xcorr-fit_rv.dat" + \
    ' --specdatafilename=' + "spec_rv.dat"
    #    else :
    #        plotstring = ''
    
    commandline = Dirs.EXE + 'operaRadialVelocity --outputRVFile=' + product + \
    ' --inputObjectSpectrum=' + inputSpectrum + ' --inputWaveFile=' + wave + \
    ' --telluric_lines=' + config.TELLURICLINES + \
    ' --inputWavelengthMask=' + config.RVXCORRWAVELENGTHMASK + \
    ' --inputStellarSpectrum=' + config.SYNTHETICSPECTRUM + \
    ' --spectralResolution=' + str(Instmode.SPECTRALRESOLUTION) + \
    ' --radialVelocitySearchRange=' + str(Instmode.RADIALVELOCITYSEARCHRANGE) + \
    ' --radialVelocitySearchStep=' + str(Instmode.RADIALVELOCITYSEARCHSTEP) + \
    ' --XCorrelationThreshold=0.1 --normalizationBinsize=110' + \
    ' --inputBarycentricCorrection=' + rvel + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --useFitToFindMaximum=1' + \
    ' --StarPlusSky=' + str(Instmode.STARPLUSKYMODEFLAG) + \
    ' --headerData="' + headerinfostrg + '"' + plotstring
    
    return commandline
##########################################


#### Function to generate a command line for Barycentric Wavelength Correction: ####
def HeliocentricWaveCommand(Dirs, product, wave, ra, dec, mjdate, exptime, startHA) :
    
    commandline = Dirs.EXE + 'operaHeliocentricWavelengthCorrection --outputRVelFile=' + product + \
    ' --inputWaveFile=' + wave + ' --observatory_coords="19:49:36 -155:28:18" --observatory_elevation=4207' + \
    ' --object_coords="' + str(ra) + ' ' + str(dec) + '" --MJDTime=' + str(mjdate) + ' --etime=' + str(exptime) + " --ha_start=" + str(startHA)
    
    return commandline
##########################################

#### Function to generate a command line to create a flux calibration spectrum *.fcal: ####
def CreateFcalCommand(Dirs, product, inputSpectrum, stdcaldata, flatSpectrum, config, Instmode, aperture, wave, exptime):
    
    commandline = Dirs.EXE + 'operaCreateFluxCalibration --outputFluxCalibrationFile=' + product + \
    ' --inputUncalibratedSpectrum=' + inputSpectrum + ' --inputCalibratedSpectrum=' + stdcaldata  + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --inputWavelengthMaskForRefContinuum=' + config.ATLASWAVELENGTHMASK + \
    ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK + ' --inputWaveFile=' + wave + \
    ' --inputApertureFile=' + aperture  + ' --wavelengthForNormalization=' + str(Instmode.WAVELENGTHFORNORMALIZATION) + ' --exposureTime=' + str(exptime) + \
    ' --numberOfPointsInUniformSample=150 --numberOfPointsInUniformRefSample=70 --binsize=500'
    
    return commandline
##########################################

#### Function to generate a command line to create a LE flat response file *.s: ####
def CreateFlatResponseCommand(Dirs, product, inputSpectrum, inputImage, stdcaldata, flatSpectrum, config, Instmode, aperture, wave):
    
    commandline = Dirs.EXE + 'operaCreateFlatResponse --outputFlatResponseFile=' + product + \
    ' --inputUncalibratedSpectrum=' + inputSpectrum + ' --inputSpectrumFITSImage=' + inputImage + \
    ' --outputFITS' + ' --inputCalibratedSpectrum=' + stdcaldata  + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --inputWavelengthMaskForRefContinuum=' + config.ATLASWAVELENGTHMASK + \
    ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK + ' --inputWaveFile=' + wave + \
    ' --wavelengthForNormalization=' + str(Instmode.WAVELENGTHFORNORMALIZATION) + \
    ' --numberOfPointsInUniformSample=300 --numberOfPointsInUniformRefSample=100 --binsize=800'
    
    return commandline
##########################################

#### Function to generate a command line to create a master flux calibration spectrum *.fcal: ####
def MasterFcalCommand(Dirs, product, inputFcalFiles, inputRefSpectrum, wave):
    
    inputfileentries = ""
    for file in inputFcalFiles :
        inputfileentries += ' --inputfcal=' + file
    
    commandline = Dirs.EXE + 'operaMasterFluxCalibration --outputfcal=' + product + \
    inputfileentries + ' --inputWaveFile=' + wave + ' --inputReferenceSpectrum=' + inputRefSpectrum + \
    ' --combineMethod=1'
    
    return commandline
##########################################

#### Function to generate a command line for calibrated spectrum *.spc: ####
def SpcModuleCommand(Dirs, product, Instmode, config, inputSpectrum, flatSpectrum, inputFcal, rvelwave, tellwave, wave, objectname, exptime) :
    
    commandline = Dirs.EXE + Instmode.SPCMODULE + ' --outputCalibratedSpectrum=' + product + \
    ' --inputUncalibratedSpectrum=' + inputSpectrum + ' --inputFlatFluxCalibration=' + flatSpectrum  + \
    ' --fluxCalibration=' + inputFcal + ' --flatResponse=' + Dirs.CONFIGDIR + Instmode.GRACESFLATRESPONSE + \
    ' --radialvelocitycorrection=' + rvelwave + ' --telluriccorrection=' + tellwave + ' --wavelengthCalibration=' + wave +\
    ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK + ' ' + Instmode.INVERTSKYFIBERFLAG + \
    ' --object="' + objectname + '" --etime=' + str(exptime) + \
    ' --spectrumtype=17 --numberOfPointsInUniformSample=150 --normalizationBinsize=750 --AbsoluteCalibration=0'
    
    return commandline
##########################################

#### Function to generate a command line for converting from *.spc to Libre-Esprit format: ####
def GenLEFormatsCommand(Dirs, config, product,inputSpectrum, LESPECTRUMTYPE, objectname, fluxtype, wavetype) :
    
    commandline = Dirs.EXE + 'operaGenerateLEFormats --outputLEfilename=' + product + \
    ' --inputOperaSpectrum=' + inputSpectrum + ' --LibreEspritSpectrumType=' + str(LESPECTRUMTYPE) + \
    ' --object="' + objectname + '" --fluxType=' + str(fluxtype) + ' --wavelengthType=' + str(wavetype) + \
    ' --LEorderwavelength=' + config.LEORDERWAVELENGTH
    
    return commandline
##########################################

#### Function to generate the final product in FITS format: ####
def CreateLEFITSProductCommand(Dirs, product,inputUS, inputUN, inputUW, inputNW, objectfile, rvelwave) :
    
    commandline = Dirs.EXE + 'operaCreateProduct --output=' + product + \
    ' --input=' + objectfile + \
    ' --ufile=' + inputUS + \
    ' --nfile=' + inputUN + \
    ' --uwfile=' + inputUW + \
    ' --nwfile=' + inputNW + \
    ' --rvel=' + rvelwave + \
    ' --spectrumtype=21 --centralsnr --compressiontype=21'
    
    return commandline
##########################################

#### Function to generate a command line for raw polarimetry *.p.gz: ####
def PolarCommand(Dirs, product, INPUT1, INPUT2, INPUT3, INPUT4, wave, stokes) :
    
    commandline = Dirs.EXE + 'operaPolar --output=' + product + \
    ' --input1=' + INPUT1 + ' --input2=' + INPUT2 + ' --input3=' + INPUT3 + ' --input4=' + INPUT4 + \
    ' --inputWaveFile=' + wave + ' --stokesparameter=' + str(stokes) + \
    ' --numberofexposures=4 --method=2 --ordernumber=-999 ' 
    
    return commandline
##########################################

#### Function to generate a command line for polarimetry calibrated product *.pol.gz : ####
def CalibratedPolarCommand(Dirs, product, Instmode, config, polar, flatSpectrum, inputFcal, rvelwave, tellwave, wave, objectname, exptime) :
    
    commandline = Dirs.EXE + 'operaPolarimetryCorrection --outputCalibratedSpectrum=' + product + \
    ' --polar=' + polar + ' --inputFlatFluxCalibration=' + flatSpectrum + \
    ' --fluxCalibration=' + inputFcal + ' --flatResponse=' + config.OLAPAFLATRESPONSE + \
    ' --radialvelocitycorrection=' + rvelwave + ' --telluriccorrection=' + tellwave + ' --wavelengthCalibration=' + wave + \
    ' --object="' + objectname + '" --etime=' + exptime + ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK + \
    ' --numberOfPointsInUniformSample=150 --normalizationBinsize=750 --AbsoluteCalibration=0'
    
    return commandline
##########################################

############# Class to encapsulate modes for reduction ####################
# This class contains the information on all available  graces reduction modes
class ReductionModes :
    'Common base class for reduction modes'
    
    # modes = {[modekey]: [modeName, InstMode, ReadMode, Nobjects, Nbiases, Nflats, Ncomps] }
    modes = {}
    files = {}
    
    def __init__(self, dirs, keywords, allowanyreadout, forcecalibration):
        if os.path.exists(dirs.DATADIR) :
            self.displayStats = "---\n"
            self.displayStats += "STATISTICS for DATA in DIR: " + dirs.DATADIR + "\n"
            self.displayStats += "Allow any readout? " + str(bool(allowanyreadout)) + "\n"
            self.displayStats += "Include mode for calibration alone? " + str(forcecalibration) + "\n"
            self.displayStats += "---\n"
            self.displayStats += "InstMode\tReadout\tobject\tbias\tflat\tcomp\tSELECTED?\n"
            self.displayStats += "-------------------------------------------------------------------\n"
            
            for mode in range(1,3) :
                for read in range(1,4) :
                    
                    Instmode = InstMode(mode)
                    Readmode = ReadoutMode(read)
                    
                    modeKey = Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
                    
                    nobjects = 0
                    nbiases = 0
                    nflats = 0
                    ncomps = 0
                    
                    objectcheck = ObjectListCommand(dirs, Instmode, Readmode, keywords.OBJECTKEYWORD)
                    objects = (subprocess.check_output(objectcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(objects) :
                        nobjects = len(objects)
                    
                    biascheck = testBiasListCommand(dirs, Instmode, Readmode, keywords)
                    biases = (subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(biases) :
                        nbiases = len(biases)
                    
                    flatcheck = testCalibrationListCommand(dirs, Instmode, Readmode,  keywords.FLATKEYWORD, allowanyreadout)
                    flats = (subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(flats) :
                        nflats = len(flats)
                    
                    compcheck = testCalibrationListCommand(dirs, Instmode, Readmode,  keywords.COMPKEYWORD, allowanyreadout)            
                    comps = (subprocess.check_output(compcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(comps) :
                        ncomps = len(comps)
                    
                    selected = "NO"
                    
                    if forcecalibration :
                        if nbiases and nflats and ncomps :
                            self.modes[modeKey] = [modeKey, mode, read, nobjects, nbiases, nflats, ncomps]
                            self.files[modeKey] = [objects, biases, flats, comps]
                            selected = "YES"
                    else :
                        if nobjects and nbiases and nflats and ncomps :
                            self.modes[modeKey] = [modeKey, mode, read, nobjects, nbiases, nflats, ncomps]
                            self.files[modeKey] = [objects, biases, flats, comps]
                            selected = "YES"
                    
                    xtab = "\t"
                    self.displayStats += Instmode.INSTRUMENTMODESHORTNAME + xtab + Readmode.READOUTSPEEDSHORTNAME + "\t" + str(nobjects) + "\t" + str(nbiases) + "\t" + str(nflats) + "\t" + str(ncomps) + "\t" + selected + "\n"
        else :
            print "Error: can\'t find data directory: ", dirs.DATADIR
        
        self.displayStats += "-------------------------------------------------------------------\n"
    
    def getNmodes(self) :
        return len(self.modes)
    
    def displayAllModeStats(self) :
        modeItems = self.modes.items()
        for modeItem in modeItems:
            print "MODE SELECTED: " + modeItem[1][0] + "\n" +\
                str(modeItem[1][3]) + " objects  " +  \
                str(modeItem[1][4]) + " biases  " +  \
                str(modeItem[1][5]) + " flats  " + \
                str(modeItem[1][6]) + " comps  "
                    
    def displayModeStats(self,intrumentmode,readoutspeed) :
        modeItems = self.modes.items()
        for modeItem in modeItems:
            if(modeItem[1][1] == intrumentmode and modeItem[1][2] == readoutspeed) :
                print "MODE SELECTED: " + modeItem[1][0] + "\n" +\
                str(modeItem[1][3]) + " objects  " +  \
                str(modeItem[1][4]) + " biases  " +  \
                str(modeItem[1][5]) + " flats  " + \
                str(modeItem[1][6]) + " comps  "
    
    def displayOverallStats(self) :
        print self.displayStats
    
    def getInstReadModes(self) :
        listmodes = []
        modeItems = self.modes.items()
        for modeItem in modeItems:
            listmodes.append([modeItem[1][1],modeItem[1][2]])
        return listmodes
    
    def displayObjectData(self) :
        filesItems = self.files.items()
        for filesItem in filesItems:
            print "-- "
            print "OBJECTS for MODE: " + filesItem[0]
            for file in filesItem[1][0] :
                print file
    
    def displayCalibrationData(self) :
        filesItems = self.files.items()
        for filesItem in filesItems:
            print "-- "            
            print "CALIBRATION DATA for MODE: " + filesItem[0]
            print "-- "            
            print "BIAS: "
            for file in filesItem[1][1] :
                print file
            print "-- "            
            print "FLAT: "
            for file in filesItem[1][2] :
                print file
            print "-- "        
            print "COMPARISON: "
            for file in filesItem[1][3] :
                print file

    def cleanModes(self) :
        for key, item in self.modes.items():
            del self.modes[key]
            del self.files[key]

###############################################################
