#
# espadons build makefile Generator for analysis espadons
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/local/lib/ -L/usr/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/local/include/ -L/usr/local/lib/ -L/usr/lib/
AM_LDFLAGS = -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaInstrumentEnvironmentSetup -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -lPixelSet -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lz -lsofa_c -lpthread -lm
# This is for Linux...
LIBS = -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaInstrumentEnvironmentSetup -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lsofa_c -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS = operaBinPolarData operaRadialVelocity operaBinFluxData

operaBinPolarData_SOURCES = operaBinPolarData.cpp

operaRadialVelocity_SOURCES = operaRadialVelocity.cpp operaRadialVelocity.h

operaBinFluxData_SOURCES = operaBinFluxData.cpp operaBinFluxData.h

