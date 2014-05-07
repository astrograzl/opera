#
# espadons build makefile Generator
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/ -DSTATIC=static
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -DSTATIC=static
AM_LDFLAGS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaWIRCamImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -lPixelSet -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lz -lsofa_c -lpthread -lm
# This is for Linux...
LIBS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaWIRCamImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lsofa_c -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS = operaSNR operaWavelengthCalibration \
			operaMasterFlat operaPolar operaBias operaExtraction \
			operaGain operaReductionSet operaGeometryCalibration \
			operaMasterComparison operaMasterBias \
			operaMasterFabPerot operaInstrumentProfileCalibration operaPixelSensitivityMap \
			operaNormalize operaCreateProduct operaExtractSpectralLines operaExtractionApertureCalibration \
			operaStarPlusSky operaStarOnly operaPolarIntensity operaFluxCalibration \
			operaCreateFluxCalibration operaTelluricWavelengthCorrection operaBarycentricWavelengthCorrection \
			operaOrderSpacingCalibration operaEchelleDispersionCalibration operaPolarimetryCorrection \
			operaCalculateSpectralResolution operaMasterFluxCalibration \
			operaObjectInTheSkySetup operaObservingConditionsSetup operaSpectrographSetup \
			operaTelescopeSetup operaNormalizeAcrossOrders operaCreateFlatFieldFluxCalibration \
			operaTelluricCorrection operaStitchOrders operaGenerateLEFormats

#
# wcs support
#				
if LIBWCS
AM_LDFLAGS += -lwcs
LIBS += -lwcs
endif

operaGenerateLEFormats_SOURCES = operaGenerateLEFormats.cpp operaGenerateLEFormats.h

operaStitchOrders_SOURCES = operaStitchOrders.cpp operaStitchOrders.h

operaTelluricCorrection_SOURCES = operaTelluricCorrection.cpp operaTelluricCorrection.h

operaCreateFlatFieldFluxCalibration_SOURCES = operaCreateFlatFieldFluxCalibration.cpp operaCreateFlatFieldFluxCalibration.h

operaNormalizeAcrossOrders_SOURCES = operaNormalizeAcrossOrders.cpp operaNormalizeAcrossOrders.h

operaTelescopeSetup_SOURCES = operaTelescopeSetup.cpp operaTelescopeSetup.h

operaSpectrographSetup_SOURCES = operaSpectrographSetup.cpp operaSpectrographSetup.h

operaObservingConditionsSetup_SOURCES = operaObservingConditionsSetup.cpp operaObservingConditionsSetup.h

operaObjectInTheSkySetup_SOURCES = operaObjectInTheSkySetup.cpp operaObjectInTheSkySetup.h

operaMasterFluxCalibration_SOURCES = operaMasterFluxCalibration.cpp operaMasterFluxCalibration.h

operaCalculateSpectralResolution_SOURCES = operaCalculateSpectralResolution.cpp operaCalculateSpectralResolution.h

operaPolarimetryCorrection_SOURCES = operaPolarimetryCorrection.cpp operaPolarimetryCorrection.h

operaEchelleDispersionCalibration_SOURCES = operaEchelleDispersionCalibration.cpp operaEchelleDispersionCalibration.h

operaOrderSpacingCalibration_SOURCES = operaOrderSpacingCalibration.cpp operaOrderSpacingCalibration.h

operaCreateFluxCalibration_SOURCES = operaCreateFluxCalibration.cpp operaCreateFluxCalibration.h

operaFluxCalibration_SOURCES = operaFluxCalibration.cpp operaFluxCalibration.h

operaExtractSpectralLines_SOURCES = operaExtractSpectralLines.cpp operaExtractSpectralLines.h

operaStarPlusSky_SOURCES = operaStarPlusSky.cpp operaStarPlusSky.h

operaPolarIntensity_SOURCES = operaPolarIntensity.cpp operaPolarIntensity.h

operaStarOnly_SOURCES = operaStarOnly.cpp operaStarOnly.h

operaExtraction_SOURCES = operaExtraction.cpp operaExtraction.h

operaNormalize_SOURCES = operaNormalize.cpp operaNormalize.h

operaCreateProduct_SOURCES = operaCreateProduct.cpp operaCreateProduct.h

operaGeometryCalibration_SOURCES = operaGeometryCalibration.cpp operaGeometryCalibration.h

operaSNR_SOURCES = operaSNR.cpp operaSNR.h

operaBias_SOURCES = operaBias.cpp operaBias.h

operaMasterFlat_SOURCES = operaMasterFlat.cpp operaMasterFlat.h
			
operaMasterComparison_SOURCES = operaMasterComparison.cpp operaMasterComparison.h
			
operaMasterBias_SOURCES = operaMasterBias.cpp operaMasterBias.h
			
operaMasterFabPerot_SOURCES = operaMasterFabPerot.cpp operaMasterFabPerot.h
			
operaPolar_SOURCES = operaPolar.cpp operaPolar.h

operaGain_SOURCES = operaGain.cpp operaGain.h
			
operaReductionSet_SOURCES = operaReductionSet.cpp operaReductionSet.h
			
operaWavelengthCalibration_SOURCES = operaWavelengthCalibration.cpp operaWavelengthCalibration.h
			
operaInstrumentProfileCalibration_SOURCES = operaInstrumentProfileCalibration.cpp operaInstrumentProfileCalibration.h

operaExtractionApertureCalibration_SOURCES = operaExtractionApertureCalibration.cpp operaExtractionApertureCalibration.h

operaPixelSensitivityMap_SOURCES = operaPixelSensitivityMap.cpp operaPixelSensitivityMap.h

operaTelluricWavelengthCorrection_SOURCES = operaTelluricWavelengthCorrection.cpp operaTelluricWavelengthCorrection.h

operaBarycentricWavelengthCorrection_SOURCES = operaBarycentricWavelengthCorrection.cpp operaBarycentricWavelengthCorrection.h

