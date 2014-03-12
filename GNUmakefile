
name := iceblock
G4TARGET := $(name)
G4EXLIB := true

ifndef G4WORKDIR
G4WORKDIR := .
endif


#EXTRALIBS := $(shell $(ROOTSYS)/bin/root-config --libs) -lReflex -lCintex 

.PHONY: all
all: lib bin

######################
# Add Root Clases
######################

#dict.cxx :
#	rootcint -f dict.cxx -c include/Analysis.h include/LinkDef.h

#include $(G4INSTALL)/config/architecture.gmk

#G4NOHIST := true
#CPPFLAGS += $(shell root-config --cflags)
# Here comes dict.cxx
#LDFLAGS  += $(shell root-config --libs) dict.cxx

include $(G4INSTALL)/config/binmake.gmk
