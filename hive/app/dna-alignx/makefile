#/*
# *  ::718604!
# * 
# * Copyright(C) November 20, 2014 U.S. Food and Drug Administration
# * Authors: Dr. Vahan Simonyan (1), Dr. Raja Mazumder (2), et al
# * Affiliation: Food and Drug Administration (1), George Washington University (2)
# * 
# * All rights Reserved.
# * 
# * The MIT License (MIT)
# * 
# * Permission is hereby granted, free of charge, to any person obtaining
# * a copy of this software and associated documentation files (the "Software"),
# * to deal in the Software without restriction, including without limitation
# * the rights to use, copy, modify, merge, publish, distribute, sublicense,
# * and/or sell copies of the Software, and to permit persons to whom the
# * Software is furnished to do so, subject to the following conditions:
# * 
# * The above copyright notice and this permission notice shall be included
# * in all copies or substantial portions of the Software.
# * 
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# * DEALINGS IN THE SOFTWARE.
# */

PROGRAMS = dna-alignx

PROGRAMS_LIBS = -lviolib -ldmlib $(VLIB_LIB)
PROGRAMS_XLIBS = $(MYSQL_LIB) $(LIBZIP_LIB)

dna-alignx_SRCS = dna-alignx.cpp

ifeq ($(devdir),)
    devdir = $(abspath ../../..)
endif
include $(devdir)/build/makefile.meta

DAYS_TO_KEEP ?= 90

# no SPACES around '=' signs in dist-* functions!!
PNOHY: dist
dist::
	$(call dist-service-register,dna-alignx,HIVE-optimized External Alignment,/Sequencing/, \
        maxJobs=4 \
    )
	$(call dist-service-set-config,dna-alignx, \
        indexDirectory='@global_scratch_dir@data/' \
        singleSubjectIndexPersistent='1' \
    )
	$(call dist-app,dna-alignx,dna-alignx)
	$(call dist-src,dna-alignx,dna-alignx.sh.os$(SLIB_PLATFORM),dna-alignx.sh)
	$(call dist-src,dna-alignx,dna-alignx-blast-common.sh.os$(SLIB_PLATFORM),dna-alignx-blast-common.sh)
	$(call dist-src,dna-alignx,dna-alignx-blast.sh.os$(SLIB_PLATFORM),dna-alignx-blast.sh)
	$(call dist-src,dna-alignx,dna-alignx-blastx.sh.os$(SLIB_PLATFORM),dna-alignx-blastx.sh)
	$(call dist-src,dna-alignx,dna-alignx-tblastx.sh.os$(SLIB_PLATFORM),dna-alignx-tblastx.sh)
	$(call dist-src,dna-alignx,dna-alignx-blat.sh.os$(SLIB_PLATFORM),dna-alignx-blat.sh)
	$(call dist-src,dna-alignx,dna-alignx-bowtie.sh.os$(SLIB_PLATFORM),dna-alignx-bowtie.sh)
	$(call dist-src,dna-alignx,dna-alignx-bowtie2.sh.os$(SLIB_PLATFORM),dna-alignx-bowtie2.sh)
	$(call dist-src,dna-alignx,dna-alignx-bwa.sh.os$(SLIB_PLATFORM),dna-alignx-bwa.sh)
	$(call dist-src,dna-alignx,dna-alignx-clustal.sh.os$(SLIB_PLATFORM),dna-alignx-clustal.sh)
	$(call dist-src,dna-alignx,dna-alignx-mafft.sh.os$(SLIB_PLATFORM),dna-alignx-mafft.sh)
	$(call dist-src,dna-alignx,dna-alignx-tophat.sh.os$(SLIB_PLATFORM),dna-alignx-tophat.sh)
