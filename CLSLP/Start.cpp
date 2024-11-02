///////////////////////////////////////////////////////////////////////////////////////////////////
//  FilePath:  CLSLP_20241021_122444/CLSLP/Start.cpp
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

#include  <fftw3.h>
#include  <iostream>
#include  <fstream>
#include  <chrono>
#include  "CLSLP.h"


double  BlackmanHarrisWindow( const double K, const double N );


int  main( int /* NumArgs */, char** /* ArgList */ )
{
  const uint32_t           FftSize      = 1U << 12;
  const uint32_t           NumFiltCoefs = 2047U;
  vector<complex<double>>  OutCoefs(NumFiltCoefs);

  std::ofstream  OutFiltCoefs("FilterCoefs.dat");
  std::ofstream  OutSpectral ("Spectral.dat");
  std::ofstream  OutPhase    ("Phase.dat");

  std::cout << ">>[FilterSynthesis: Example Filters]" << std::endl;

  {
    // CLSLP LPFIR Filter Synthesis
    CLSLP  Filt;

    // This is an amplitude, not dB, that multiplies A0 and A1 in the Filt.Add() table below
    Filt.SetGain(1);    // Not dB, optional call --- Gain assumed 1 otherwise

    // ‘AddSymmetric’ is for Real-Valued Filter Coefficients    ([0,.5] mod 1.) Normalized Hz
    // ‘Add’          is for Complex-Valued Filter Coefficients ([0,1.) mod 1.) Normalized Hz
    //          Type                F0     F1     A0     A1   AuxMult
    Filt.Add( CLSLP::eLin_AbsErr,  0.00,  0.10,  0.00,  0.00,  1.00 );
    Filt.Add( CLSLP::eLin_AbsErr,  0.10,  0.20,  1.00,   1e5,  1.00 );
    Filt.Add( CLSLP::eLin_AbsErr,  0.20,  0.30,  0.00,  0.00,  1.00 );
    Filt.Add( CLSLP::eExp_AbsErr,  0.30,  0.40,   1e5,  1.00,  1.00 );
    Filt.Add( CLSLP::eLin_AbsErr,  0.40,  0.50,  0.00,  0.00,  1.00 );
    Filt.Add( CLSLP::eLin_AbsErr,  0.50,  0.60,  0.00,  0.00,  1.00 );
    Filt.Add( CLSLP::eLin_RelErr,  0.60,  0.70,  1.00,   1e5,  1.00 );
    Filt.Add( CLSLP::eLin_AbsErr,  0.70,  0.80,  0.00,  0.00,  1.00 );
    Filt.Add( CLSLP::eExp_RelErr,  0.80,  0.90,   1e5,  1.00,  1.00 );
    Filt.Add( CLSLP::eLin_AbsErr,  0.90,  1.00,  0.00,  0.00,  1.00 );

    // Generate the filter coefficients
    const auto Clock_0 = chrono::high_resolution_clock::now();

    // Calculate the FIR Filter Coeffcients
    Filt.GenFilter( NumFiltCoefs, OutCoefs );

    const auto Clock_1 = chrono::high_resolution_clock::now();

    for( auto K=0U ; K < NumFiltCoefs ; ++K )
    {
      OutFiltCoefs << OutCoefs[K] << std::endl;
    }

    const double DeltaTime =
        double(chrono::duration_cast <chrono::microseconds> (Clock_1 - Clock_0).count());

    cout << "--FIR Filter Coefficient Generation: NumCoefs:"
         << NumFiltCoefs << "  [TotalSynthPeriod:"
         << DeltaTime << " uS  PerCoefSynthPeriod:"
         << (DeltaTime/NumFiltCoefs) << " uS/Coef]"
         << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // PLOT GENERATION FOLLOWS // Not needed for CLSLP LPFIR Filter Synthesis //
    ////////////////////////////////////////////////////////////////////////////

    cout << ">> GENERATING PLOTS" << std::endl;

    // Take DFT of these coefficents
    {
      fftw_complex  *In=nullptr, *Out=nullptr;
      fftw_plan     Plan=nullptr;

      In  = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * FftSize );
      Out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * FftSize );

      // Plan the FFT
      Plan = fftw_plan_dft_1d( FftSize, In, Out, FFTW_FORWARD,  FFTW_ESTIMATE );

      // Fill input data.
      //   No need to iterate on [0,FftSize).  Only (0,NumFiltCoefs) non-zero.
      for( auto K=0U ; K < NumFiltCoefs ; ++K )
      {
        const double  BHW = BlackmanHarrisWindow( K, NumFiltCoefs );

        In[K][0] = BHW * real(OutCoefs[K]);
        In[K][1] = BHW * imag(OutCoefs[K]);
      }

      for( auto K=NumFiltCoefs ; K < FftSize ; ++K )
      {
        In[K][0] = 0;
        In[K][1] = 0;
      }

      // Perform the FFT
      fftw_execute(Plan);

      for( int32_t  K = 0 ; K < (signed) FftSize ; ++K )
      {
        const double  A = ( (double) K / FftSize );
        double        Mag = 0;

        // Amplitude
        if( K >= 0 )
        {
           Mag = 10.*log10(Out[K][0]*Out[K][0] + Out[K][1]*Out[K][1]);
        }
        else
        {
          assert( ( K + (signed) FftSize ) >= 0 );

          const uint32_t  K2 = K + FftSize;

          Mag = 10.*log10(Out[K2][0]*Out[K2][0] + Out[K2][1]*Out[K2][1]);
        }

        OutSpectral <<  A  << "  " << Mag << std::endl;

        // Phase Response
        double  Phase = 0;

        // Phase
        if( K >= 0 )
        {
          Phase = atan2( Out[K][1], Out[K][0] );
        }
        else
        {
          Phase = atan2( Out[K+FftSize-1U][1], Out[K+FftSize-1U][0] );
        }

        OutPhase <<  A  << "  " << Phase << std::endl;
      }

      // Clean Up
      fftw_destroy_plan(Plan);
      fftw_free(In);
      fftw_free(Out);
    }
  }

  std::cout << "<<[FilterSynthesis]" << std::endl;

  return  0;
}


double  BlackmanHarrisWindow( const double Num, const double Den )
{
  const double  Mult = M_PI * Num / Den;

  const double  BHD =   0.35875
                      - 0.48829 * cos(2*Mult)
                      + 0.14128 * cos(4*Mult)
                      - 0.01168 * cos(6*Mult);

  return  BHD;
}

