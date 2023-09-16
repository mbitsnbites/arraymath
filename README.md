## This repo has moved to: https://gitlab.com/mbitsnbites/arraymath

# ArrayMath

ArrayMath is a portable low level C++ library for performing primitive math
operations on arrays.

The idea is to provide a simple interface to highly optimized math kernels.


## Usage

To build the static library, a test program and a benchmark program, run `make`
in the `src/` folder (requires GCC).

It is also possible to cross compile for ARM (`make ARCH=arm`), which can be
useful for testing in a QEMU environment, for instance.

An experimental MSVC makefile (`nmake /F Makefile.msvc`) can be used for
building the test and benchmark programs from command line.

Provided that you've built the static library or separate object files
that you link to your project, you will find the necessary information in
the doxygen documentation (run `doxygen` in the `doc/` folder to generate
the HTML documentation).


## SIMD acceleration

A key feature of ArrayMath is that it has been optimized for different SIMD
architectures. There are currently x86 (SSE, AVX, etc) and ARM (NEON) optimized
versions of many routines.

ArrayMath will select the most appropriate version of each routine in the
library at runtime, depending on the CPU it is running on.


## ArrayMath license

ArrayMath is released under the zlib/libpng license:

```
 Copyright (c) 2013-2014 Marcus Geelnard

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from the
 use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not claim
    that you wrote the original software. If you use this software in a
    product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
```

## Third party licenses

### Kiss FFT

ArrayMath uses [Kiss FFT](http://kissfft.sourceforge.net/) by Mark Borgerding,
which is released under a BSD license:

```
 Copyright (c) 2003-2010 Mark Borgerding

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to
      endorse or promote products derived from this software without specific
      prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### sse\_mathfun & neon\_mathfun

ArrayMath uses [sse\_mathfun and neon\_mathfun](http://gruntthepeon.free.fr/ssemath/)
by Julien Pommier, which are released under the zlib/libpng license:

```
 Copyright (C) 2007 / 2011 Julien Pommier

 This software is provided 'as-is', without any express or implied
 warranty.  In no event will the authors be held liable for any damages
 arising from the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.
 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.
 3. This notice may not be removed or altered from any source distribution.
```

