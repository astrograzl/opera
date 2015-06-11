#
# espadons build makefile Generator
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/
AM_LDFLAGS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -lPixelSet -loperaStats -loperaLib -lArgumentHandler -loperaArgumentHandler -loperaCommonModuleElements -lfftw3 -lgzstream -lcfitsio -lz -lsofa_c -lpthread -lm
# This is for Linux...
LIBS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaStats -loperaLib -lArgumentHandler -loperaArgumentHandler -loperaCommonModuleElements -lfftw3 -lgzstream -lcfitsio -lsofa_c -lz -lpthread -lm

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
			operaStitchOrders operaGenerateLEFormats operaCreateFlatResponse operaHeliocentricWavelengthCorrection

#
# wcs support
#				
if LIBWCS
AM_LDFLAGS += -lwcs
LIBS += -lwcs
endif

operaCreateFlatResponse_SOURCES = operaCreateFlatResponse.cpp

operaGenerateLEFormats_SOURCES = operaGenerateLEFormats.cpp

operaStitchOrders_SOURCES = operaStitchOrders.cpp

operaCreateFlatFieldFluxCalibration_SOURCES = operaCreateFlatFieldFluxCalibration.cpp

operaNormalizeAcrossOrders_SOURCES = operaNormalizeAcrossOrders.cpp operaNormalizeAcrossOrders.h

operaTelescopeSetup_SOURCES = operaTelescopeSetup.cpp operaTelescopeSetup.h

operaSpectrographSetup_SOURCES = operaSpectrographSetup.cpp operaSpectrographSetup.h

operaObservingConditionsSetup_SOURCES = operaObservingConditionsSetup.cpp operaObservingConditionsSetup.h

operaObjectInTheSkySetup_SOURCES = operaObjectInTheSkySetup.cpp operaObjectInTheSkySetup.h

operaMasterFluxCalibration_SOURCES = operaMasterFluxCalibration.cpp

operaCalculateSpectralResolution_SOURCES = operaCalculateSpectralResolution.cpp operaCalculateSpectralResolution.h

operaPolarimetryCorrection_SOURCES = operaPolarimetryCorrection.cpp operaExtendedSpectrumCreation.cpp operaExtendedSpectrumCreation.h

operaEchelleDispersionCalibration_SOURCES = operaEchelleDispersionCalibration.cpp operaEchelleDispersionCalibration.h

operaOrderSpacingCalibration_SOURCES = operaOrderSpacingCalibration.cpp

operaCreateFluxCalibration_SOURCES = operaCreateFluxCalibration.cpp

operaFluxCalibration_SOURCES = operaFluxCalibration.cpp operaFluxCalibration.h

operaExtractSpectralLines_SOURCES = operaExtractSpectralLines.cpp operaExtractSpectralLines.h

operaStarPlusSky_SOURCES = operaStarPlusSky.cpp operaExtendedSpectrumCreation.cpp operaExtendedSpectrumCreation.h

operaPolarIntensity_SOURCES = operaPolarIntensity.cpp operaExtendedSpectrumCreation.cpp operaExtendedSpectrumCreation.h

operaStarOnly_SOURCES = operaStarOnly.cpp operaExtendedSpectrumCreation.cpp operaExtendedSpectrumCreation.h

operaExtraction_SOURCES = operaExtraction.cpp

operaNormalize_SOURCES = operaNormalize.cpp operaNormalize.h

operaCreateProduct_SOURCES = operaCreateProduct.cpp

operaGeometryCalibration_SOURCES = operaGeometryCalibration.cpp

operaSNR_SOURCES = operaSNR.cpp

operaBias_SOURCES = operaBias.cpp

operaMasterFlat_SOURCES = operaMasterFlat.cpp operaMasterCalibration.cpp operaMasterCalibration.h
			
operaMasterComparison_SOURCES = operaMasterComparison.cpp
			
operaMasterBias_SOURCES = operaMasterBias.cpp operaMasterCalibration.cpp operaMasterCalibration.h
			
operaMasterFabPerot_SOURCES = operaMasterFabPerot.cpp operaMasterCalibration.cpp operaMasterCalibration.h
			
operaPolar_SOURCES = operaPolar.cpp

operaGain_SOURCES = operaGain.cpp
			
operaReductionSet_SOURCES = operaReductionSet.cpp operaReductionSet.h
			
operaWavelengthCalibration_SOURCES = operaWavelengthCalibration.cpp operaWavelengthCalibration.h
			
operaInstrumentProfileCalibration_SOURCES = operaInstrumentProfileCalibration.cpp

operaExtractionApertureCalibration_SOURCES = operaExtractionApertureCalibration.cpp

operaPixelSensitivityMap_SOURCES = operaPixelSensitivityMap.cpp operaPixelSensitivityMap.h

operaTelluricWavelengthCorrection_SOURCES = operaTelluricWavelengthCorrection.cpp operaTelluricWavelengthCorrection.h

operaBarycentricWavelengthCorrection_SOURCES = operaBarycentricWavelengthCorrection.cpp operaCentricWavelengthCorrection.cpp operaCentricWavelengthCorrection.h

operaHeliocentricWavelengthCorrection_SOURCES = operaHeliocentricWavelengthCorrection.cpp operaCentricWavelengthCorrection.cpp operaCentricWavelengthCorrection.h
