# -*- Mode: Makefile; tab-width: 4; indent-tabs-mode: t -*-
#-------------------------------------------------------------------------------
# ArrayMath - an array math library
#-------------------------------------------------------------------------------
# Copyright(c) 2013-2014 Marcus Geelnard
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

# CPU architecture selection. Supported values: "x86" and "arm".
ARCH ?= x86

# FFT back end selection. Supported values: "kiss_fft".
FFT ?= kiss_fft

OUTDIR = out_$(ARCH)

# C compiler settings.
CC = gcc
CFLAGS = --std=c99 -c -O3 -ffast-math -W -Wall -Wextra -pedantic -I.

# C++ compiler settings.
CPP = g++
CPPFLAGS = -c -O3 -ffast-math -W -Wall -Wextra -pedantic -I.

# Compiler settings for the test programs.
TESTFLAGS = -c -g -W -Wall -Wextra -pedantic -I.

# Automatically generate dependency files.
CFLAGS += -MMD -MP
CPPFLAGS += -MMD -MP

# Linker settings.
LD = $(CPP)
LFLAGS = -L$(OUTDIR)
LIBS = -larraymath

# Library archiving settings.
AR = ar
ARFLAGS = rcs

RM = rm -rf

ifeq ($(ARCH), x86)
# x86 SIMD activation.
CPPFLAGS += -DAM_HAS_SSE -DAM_HAS_SSE2 -DAM_HAS_SSE4 -DAM_HAS_AVX

# x86 specific SIMD flags.
CPPFLAGS_SSE = $(CPPFLAGS) -msse
CPPFLAGS_SSE2 = $(CPPFLAGS) -msse2
CPPFLAGS_SSE4 = $(CPPFLAGS) -msse4.1
CPPFLAGS_AVX = $(CPPFLAGS) -mavx
endif

ifeq ($(ARCH), arm)
# Use cross compiling tool chain.
CC = arm-linux-gnueabihf-gcc
CPP = arm-linux-gnueabihf-g++

# ARM architecture flags.
CPPFLAGS += -march=armv7-a -mtune=generic-armv7-a

# ARM SIMD activation.
CPPFLAGS += -DAM_HAS_NEON

# ARM specific SIMD flags.
CPPFLAGS_NEON = $(CPPFLAGS) -mfloat-abi=hard -mfpu=neon
endif

ifeq ($(FFT), kiss_fft)
# Kiss FFT activation.
CFLAGS   += -Ikiss_fft
CPPFLAGS += -DAM_USE_KISS_FFT -Ikiss_fft

# Kiss FFT object files.
FFT_OBJS = $(OUTDIR)/kiss_fft.o \
           $(OUTDIR)/kiss_fftr.o \
           $(OUTDIR)/FFTKiss.o
endif

# Core ArrayMath objects.
OBJS = $(OUTDIR)/ArrayMath.o \
       $(OUTDIR)/ArrayMathGeneric.o \
       $(OUTDIR)/CPUFeatureDetector.o \
       $(OUTDIR)/FilterFactory.o \
       $(OUTDIR)/FilterGeneric.o \
       $(OUTDIR)/FFTFactory.o \
       $(OUTDIR)/RandomGeneric.o \
       $(FFT_OBJS)

ifeq ($(ARCH), x86)
OBJS += $(OUTDIR)/ArrayMathAVX.o \
        $(OUTDIR)/ArrayMathSSE.o \
        $(OUTDIR)/ArrayMathSSE4.o \
        $(OUTDIR)/FilterAVX.o \
        $(OUTDIR)/RandomSSE2.o
endif

ifeq ($(ARCH), arm)
OBJS += $(OUTDIR)/ArrayMathARM.o \
        $(OUTDIR)/ArrayMathNEON.o \
        $(OUTDIR)/RandomNEON.o
endif

BENCHMARK_OBJS = $(OUTDIR)/Benchmark.o

TEST_OBJS = $(OUTDIR)/Test.o \
            $(OUTDIR)/Tester.o \
            $(OUTDIR)/TestArrayMath.o \
            $(OUTDIR)/TestFFTFactory.o \
            $(OUTDIR)/TestFilterFactory.o \
            $(OUTDIR)/TestRandom.o

LIBRARY = $(OUTDIR)/libarraymath.a

all: $(LIBRARY) $(OUTDIR)/Benchmark $(OUTDIR)/Test

clean:
	$(RM) $(OBJS) $(OBJS:.o=.d) $(BENCHMARK_OBJS) $(BENCHMARK_OBJS:.o=.d) $(TEST_OBJS) $(TEST_OBJS:.o=.d) $(LIBRARY) $(OUTDIR)/Test $(OUTDIR)/Benchmark

# Static library build rule.
$(LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

# Test programs build rules.
$(OUTDIR)/Benchmark: $(LIBRARY) $(BENCHMARK_OBJS)
	$(CPP) $(LFLAGS) -o $@ $(BENCHMARK_OBJS) $(LIBS)

$(OUTDIR)/Test: $(LIBRARY) $(TEST_OBJS)
	$(CPP) $(LFLAGS) -o $@ $(TEST_OBJS) $(LIBS)

# Generic sources.
$(OUTDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: common/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: generic/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/%.o: test/%.cpp
	$(CPP) $(TESTFLAGS) -o $@ $<

ifeq ($(ARCH), x86)
# x86 specific sources.
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
endif

ifeq ($(ARCH), arm)
# ARM specific sources.
$(OUTDIR)/ArrayMathARM.o: arm/ArrayMathARM.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<

$(OUTDIR)/ArrayMathNEON.o: arm/ArrayMathNEON.cpp
	$(CPP) $(CPPFLAGS_NEON) -o $@ $<

$(OUTDIR)/RandomNEON.o: arm/RandomNEON.cpp
	$(CPP) $(CPPFLAGS_NEON) -o $@ $<
endif

ifeq ($(FFT), kiss_fft)
# Kiss FFT sources.
$(OUTDIR)/%.o: kiss_fft/%.c
	$(CC) $(CFLAGS) -o $@ $<

$(OUTDIR)/%.o: kiss_fft/tools/%.c
	$(CC) $(CFLAGS) -o $@ $<
endif

# These dependency files are generated automatically when building the object
# files.
-include $(OBJS:.o=.d)

