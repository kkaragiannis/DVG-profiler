#/*
# * $Id$ ::2B6A7A!
# * Copyright (c) 2005 Dr. Vahan Simonyan and Dr. Raja Mazumder.
# * This software is protected by U.S. Copyright Law and International
# * Treaties. Unauthorized use, duplication, reverse engineering, any
# * form of redistribution, use in part or as a whole, other than by
# * prior, express, written and signed agreement is subject to penalties.
# * If you have received this file in error, please notify copyright
# * holder and destroy this and any other copies. All rights reserved.
# */

PROGRAMS += dna-diprofiler

PROGRAMS_LIBS = -lviolib $(VLIB_LIB)
PROGRAMS_XLIBS = $(MYSQL_LIB)

dna-diprofiler_SRCS = dna-diprofiler.cpp

OBJECT_TYPES = svc-diprofiler

MENUITEMS = \
/Portal/Post-Alignment/Defective_Viral_Genomes_Detection

ifeq ($(devdir),)
    devdir = $(abspath ../../..)
endif
include $(devdir)/build/makefile.meta

PNOHY: dist
dist::
	$(call dist-service-register,dna-diprofiler,DI Profiler,/Sequencing/)
	$(call dist-app,dna-diprofiler,dna-diprofiler)
