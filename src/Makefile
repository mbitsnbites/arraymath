# -*- Mode: Makefile; tab-width: 4; indent-tabs-mode: t -*-
#-------------------------------------------------------------------------------
# ArrayMath - an array math library
#-------------------------------------------------------------------------------
# Copyright(c) 2013 Marcus Geelnard
#
# This software is provided 'as-is', without any express or implied warranty.
# In no event will the authors be held liable for any damages arising from the
# use of this software.
#
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
#
# 1. The origin of this software must not be misrepresented; you must not claim
#    that you wrote the original software. If you use this software in a
#    product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
#
# 2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
#
# 3. This notice may not be removed or altered from any source distribution.
#-------------------------------------------------------------------------------

.PHONY: all clean

CPP = g++
CPPFLAGS = -c -O3 -ffast-math -W -Wall -Wextra -pedantic -I.
RM = rm -rf

# Automatically generate dependency files
CPPFLAGS += -MMD -MP

# x86 SIMD activation
CPPGLAGS += -DAM_HAS_SSE -DAM_HAS_SSE2 -DAM_HAS_AVX
CPPFLAGS_SSE = $(CPPFLAGS) -msse
CPPFLAGS_SSE2 = $(CPPFLAGS) -msse -msse2
CPPFLAGS_AVX = $(CPPFLAGS) -msse -msse2 -mavx

OUTDIR = out

OBJS = $(OUTDIR)/ArrayMath.o \
       $(OUTDIR)/ArrayMathGeneric.o \
       $(OUTDIR)/ArrayMathSSE.o \
       $(OUTDIR)/ArrayMathSSE2.o \
       $(OUTDIR)/CPUFeatureDetector.o \
       $(OUTDIR)/FilterFactory.o \
       $(OUTDIR)/FilterGeneric.o \
       $(OUTDIR)/RandomGeneric.o \
       $(OUTDIR)/RandomSSE2.o

all: $(OUTDIR)/Benchmark $(OUTDIR)/Test

clean:
	$(RM) $(OBJS) $(OBJS:.o=.d) $(OUTDIR)/Test $(OUTDIR)/Benchmark

$(OUTDIR)/Benchmark: $(OBJS) $(OUTDIR)/Benchmark.o
	$(CPP) -o $@ $^

$(OUTDIR)/Test: $(OBJS) $(OUTDIR)/Test.o
	$(CPP) -o $@ $^

$(OUTDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: common/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: generic/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/ArrayMathSSE.o: x86/ArrayMathSSE.cpp
	$(CPP) $(CPPFLAGS_SSE) -o $@ $<

$(OUTDIR)/ArrayMathSSE2.o: x86/ArrayMathSSE2.cpp
	$(CPP) $(CPPFLAGS_SSE2) -o $@ $<

$(OUTDIR)/RandomSSE2.o: x86/RandomSSE2.cpp
	$(CPP) $(CPPFLAGS_SSE2) -o $@ $<

$(OUTDIR)/%.o: test/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

# These dependency files are generated automatically when building the object
# files.
-include $(OBJS:.o=.d)

