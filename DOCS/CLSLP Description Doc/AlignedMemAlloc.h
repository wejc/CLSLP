///////////////////////////////////////////////////////////////////////////////////////////////////
//  FilePath:  CLSLP_20241021_122444/DOCS/CLSLP Description Doc/AlignedMemAlloc.h
//  Tag:       CLSLP_20241021_122444
//  TimeDate:  20241021_122708
//  Email:     WEJC@WEJC.COM
//  Author:    William Earl Jones
///////////////////////////////////////////////////////////////////////////////////////////////////
//                           Copyright (c) 2024 William Earl Jones                               //
//                                                                                               //
//                            ---Standard BSD 3-Clause License---                                //
//                                                                                               //
//  Redistribution and use in source and binary forms, with or without modification, are         //
//  permitted provided that the following conditions are met:                                    //
//                                                                                               //
//    1. Redistributions of source code must retain the above copyright notice, this list of     //
//       conditions and the following disclaimer.                                                //
//                                                                                               //
//    2. Redistributions in binary form must reproduce the above copyright notice, this list of  //
//       conditions and the following disclaimer in the documentation and/or other materials     //
//       provided with the distribution.                                                         //
//                                                                                               //
//    3. Neither the name of the copyright holder nor the names of its contributors may be used  //
//       to endorse or promote products derived from this software without specific prior        //
//       written permission.                                                                     //
//                                                                                               //
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS  //
//  OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF              //
//  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE   //
//  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,    //
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF           //
//  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       //
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR     //
//  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, //
//  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                           //
/////////////////////////////////////////////////////////////////////////////////////////<wejc>////

#pragma once

#include  <cassert>
#include  <complex>
#include  <iostream>
#include  <memory>
#include  <vector>
#include  <cstdlib>
#include  <fftw3.h>

using namespace std;


// This allocator allows vectors to use FFTW's allocation routine which allows proper
//   memory alignment for SSE processing
template<typename Tdata>
class  SSE_Allocator : public std::allocator<Tdata>
{
public:
  template  <typename U>

  struct rebind { typedef SSE_Allocator<U> other; };

  Tdata* allocate(size_t N)
  {
    const auto  AlignBytes    = 16U;
    const auto  AlignShifts   = 4U;
    const auto  AllocSize0    = N * sizeof(Tdata);
    const auto  AllocSizeMod  = ((AllocSize0+(AlignBytes-1)) >> AlignShifts) << AlignShifts;

    return static_cast<Tdata*>(std::aligned_alloc( AlignBytes, AllocSizeMod ));
  }

  void deallocation(Tdata* data, std::size_t size)
  {
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
    std::free(data);
  }
};


// This allocator allows vectors to use FFTW's allocation routine which allows proper
//   memory alignment for AVX processing
template<typename Tdata>
class  AVX_Allocator : public std::allocator<Tdata>
{
public:
  template  <typename U>
  struct rebind { typedef AVX_Allocator<U> other; };

  Tdata* allocate(size_t N)
  {
    const auto  AlignBytes    = 32U;
    const auto  AlignShifts   = 5U;
    const auto  AllocSize0    = N * sizeof(Tdata);
    const auto  AllocSizeMod  = ((AllocSize0+(AlignBytes-1)) >> AlignShifts) << AlignShifts;

    return static_cast<Tdata*>(std::aligned_alloc( AlignBytes, AllocSizeMod ));
  }

  void deallocation(Tdata* data, std::size_t size)
  {
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
    std::free(data);
  }
};


// AVector properly memory aligned for SSE
template<typename T>
using SSE_Vector = vector<T,SSE_Allocator<T>>;


// AVector properly memory aligned for AVX
template<typename T>
using AVX_Vector = vector<T,AVX_Allocator<T>>;


