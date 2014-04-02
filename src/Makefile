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

OUTDIR = out
LIBRARY = $(OUTDIR)/libarraymath.a

# C compiler settings
CC = gcc
CFLAGS = --std=c99 -c -O3 -ffast-math -W -Wall -Wextra -pedantic -I. -Ikiss_fft

# C++ compiler settings
CPP = g++
CPPFLAGS = -c -O3 -ffast-math -W -Wall -Wextra -pedantic -I. -Ikiss_fft

# Compiler settings for the test programs
TESTFLAGS = -c -g -W -Wall -Wextra -pedantic -I.

# Automatically generate dependency files
CFLAGS += -MMD -MP
CPPFLAGS += -MMD -MP

# Linker settings
LD = $(CPP)
LFLAGS = -L$(OUTDIR)
LIBS = -larraymath

# Library archiving settings
AR = ar
ARFLAGS = rcs

RM = rm -rf

# x86 SIMD activation
CPPFLAGS += -DAM_HAS_SSE -DAM_HAS_SSE2 -DAM_HAS_SSE4 -DAM_HAS_AVX -DAM_USE_KISS_FFT
CPPFLAGS_SSE = $(CPPFLAGS) -msse
CPPFLAGS_SSE2 = $(CPPFLAGS) -msse2
CPPFLAGS_SSE4 = $(CPPFLAGS) -msse4.1
CPPFLAGS_AVX = $(CPPFLAGS) -mavx

# FFT objects
FFT_OBJS = $(OUTDIR)/kiss_fft.o \
           $(OUTDIR)/kiss_fftr.o \
           $(OUTDIR)/FFTKiss.o

# Core ArrayMath objects
OBJS = $(OUTDIR)/ArrayMath.o \
       $(OUTDIR)/ArrayMathAVX.o \
       $(OUTDIR)/ArrayMathGeneric.o \
       $(OUTDIR)/ArrayMathSSE.o \
       $(OUTDIR)/ArrayMathSSE4.o \
       $(OUTDIR)/CPUFeatureDetector.o \
       $(OUTDIR)/FilterAVX.o \
       $(OUTDIR)/FilterFactory.o \
       $(OUTDIR)/FilterGeneric.o \
       $(OUTDIR)/FFTFactory.o \
       $(OUTDIR)/RandomGeneric.o \
       $(OUTDIR)/RandomSSE2.o \
       $(FFT_OBJS)

BENCHMARK_OBJS = $(OUTDIR)/Benchmark.o

TEST_OBJS = $(OUTDIR)/Test.o \
            $(OUTDIR)/Tester.o \
            $(OUTDIR)/TestArrayMath.o \
            $(OUTDIR)/TestFFTFactory.o \
            $(OUTDIR)/TestFilterFactory.o \
            $(OUTDIR)/TestRandom.o

all: $(LIBRARY) $(OUTDIR)/Benchmark $(OUTDIR)/Test

clean:
	$(RM) $(OBJS) $(OBJS:.o=.d) $(BENCHMARK_OBJS) $(BENCHMARK_OBJS:.o=.d) $(TEST_OBJS) $(TEST_OBJS:.o=.d) $(LIBRARY) $(OUTDIR)/Test $(OUTDIR)/Benchmark

$(LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(OUTDIR)/Benchmark: $(LIBRARY) $(BENCHMARK_OBJS)
	$(CPP) $(LFLAGS) -o $@ $(BENCHMARK_OBJS) $(LIBS)

$(OUTDIR)/Test: $(LIBRARY) $(TEST_OBJS)
	$(CPP) $(LFLAGS) -o $@ $(TEST_OBJS) $(LIBS)

$(OUTDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: common/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: generic/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/ArrayMathSSE.o: x86/ArrayMathSSE.cpp
	$(CPP) $(CPPFLAGS_SSE) -o $@ $<

$(OUTDIR)/ArrayMathSSE4.o: x86/ArrayMathSSE4.cpp
	$(CPP) $(CPPFLAGS_SSE4) -o $@ $<

$(OUTDIR)/ArrayMathAVX.o: x86/ArrayMathAVX.cpp
	$(CPP) $(CPPFLAGS_AVX) -o $@ $<

$(OUTDIR)/FilterAVX.o: x86/FilterAVX.cpp
	$(CPP) $(CPPFLAGS_AVX) -o $@ $<

$(OUTDIR)/RandomSSE2.o: x86/RandomSSE2.cpp
	$(CPP) $(CPPFLAGS_SSE2) -o $@ $<

$(OUTDIR)/%.o: test/%.cpp
	$(CPP) $(TESTFLAGS) -o $@ $<

# Kiss FFT compilation rules
$(OUTDIR)/%.o: kiss_fft/%.c
	$(CC) $(CFLAGS) -o $@ $<

$(OUTDIR)/%.o: kiss_fft/tools/%.c
	$(CC) $(CFLAGS) -o $@ $<


# These dependency files are generated automatically when building the object
# files.
-include $(OBJS:.o=.d)

