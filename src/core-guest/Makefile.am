#########################################################################################
#
# Guest build makefile Generator
#
#########################################################################################
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/  -L../../lib -L/usr/lib/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/ -DSTATIC=static
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/  -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/local/lib/ -DSTATIC=static
AM_LDFLAGS = -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaWIRCamImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -lPixelSet -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lz -lsofa_c -lpthread -lm
# This is for Linux...
LIBS = -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaWIRCamImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lsofa_c -lz -lpthread -lm

#########################################################################################
# this lists the binaries to produce -- add all your modules here
#########################################################################################

bin_PROGRAMS = extraction

extraction_SOURCES = extraction.cpp
