#########################################################################################
# This configuration file lists all the subirectories to be traversed 
# to create Makefiles...
#########################################################################################

AUTOMAKE_OPTIONS = foreign

AM_MAKEFLAGS	= PREFIX=$(prefix)
ACLOCAL_AMFLAGS = -I m4

# add doc to this if you want to build the doxygen documentation
SUBDIRS = \
	src/libraries sofa/20120301_a/c/src src/tools test \
	src/core-espadons src/core-wircam src/core-guest \
	src/analysis-espadons src/analysis-wircam src/analysis-guest \
	scripts scripts/espadons scripts/wircam scripts/guest \
	scripts/sitelle