########################## PNGwriter #########################################
#
#   Website: Main:             http://pngwriter.sourceforge.net/
#            GitHub.com:       https://github.com/pngwriter/pngwriter
#            Sourceforge.net:  http://sourceforge.net/projects/pngwriter/
#
#
#    Author:                    Paul Blackburn https://github.com/individual61
#                               Axel Huebl https://github.com/ax3l
#
#    Email:                     individual61@users.sourceforge.net
#
#    Version:                   0.5.5 (August 2015)
#
#    Description:               Library that allows plotting a 48 bit
#                               PNG image pixel by pixel, which can
#                               then be opened with a graphics program.
#
#    License:                   GNU General Public License
#                               (C) 2002-2015 Paul Blackburn
#                               (C) 2013-2015 Axel Huebl
#
##############################################################################

######################## IMPORTANT ############################################
# This makefile is meant to help beginning programmers learn about simple 
# Makefiles. The compilation uses the object file, and not the library,
# because I belive it to be a more 'hands-on' approach.
###############################################################################


PREFIX=/usr/local
PNGPATH=../png/pngwriter

ifdef P_FREETYPE
FT_ARG_CFLAGS = -DNO_FREETYPE
FT_ARG_LIBS =
else
FT_ARG_CFLAGS = `freetype-config --cflags`
FT_ARG_LIBS = `freetype-config --libs`
endif

CXX=g++

CXXFLAGS= -O3 -Wall -Wno-deprecated $(FT_ARG_CFLAGS)

INC=  -I$(PNGPATH)/src -I$(PREFIX)/include/

LIBS= -L$(PNGPATH)/src -L$(PREFIX)/lib/ -lz -lpngwriter -lpng $(FT_ARG_LIBS)

INSTALL=install

SELF=$(PNGPATH)make.include.linux
EXAMPLES= featurevisualizer rgb

all: $(EXAMPLES)

featurevisualizer: featurevisualizer.cc
	$(CXX) $(CXXFLAGS) $(INC) featurevisualizer.cc -o featurevisualizer $(LIBS) 

rgb: rgb.cc
	$(CXX) $(CXXFLAGS) $(INC) rgb.cc -o rgb $(LIBS) 

clean	:    
	rm -f $(EXAMPLES) *~ 
	rm -f Makefile~
	rm -f .DS_Store *.png
