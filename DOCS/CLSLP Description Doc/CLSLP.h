///////////////////////////////////////////////////////////////////////////////////////////////////
//  FilePath:  CLSLP_20241021_122444/DOCS/CLSLP Description Doc/CLSLP.h
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

#include  <complex>
#include  <vector>
#include  <sstream>
#include  <list>
#include  <iomanip>
#include  <cassert>
#include  <cstdint>
#include  <gsl/gsl_sf.h>

using namespace std;


/*/////////////////////////////////////////////////////////////////////////////////////////////////
  class  CLSLP

  The CLSLP C++ class is a tool that generates the coefficients for a linear-phase FIR
  filter filters. These filters are specified in line segments of frequency and amplitude
  (i.e. freq:[F0,F1], amp: [A0,A1]).  Two types of segments: linear and exponential.  Linear
  segment’s amplitude are specified by an equation of the form B*f+A where B and A are picked
  so that the line-segment end points, freq:[F0,F1] and amp: [A0,A1], are satisfied.  The other
  type of segment is the exponential segment defined by exp(B*f+A).  In this case the end points
  are specified by the user as freq:[F0,F1], amp: [loge(A0),loge(A1)].  Please read the PDF that
  accompanies this software for examples of use and mathematical background.  These algorithms
  were developed and published by Amin G. Jaffer and William E. Jones in the 1994 timeframe.

  The Diagnostics member function gives an example of each segment type and the code associated
  with FIR filter coefficient generation and evaluation.
/////////////////////////////////////////////////////////////////////////////////////////////////*/

class  CLSLP
{
public:

  // Some constants
  static const double  TPI;
  static const double  PI2;

  // Segment types
  //  The interval is defined by [F0,F1] in frequency and [A0,A1] in amplitude for Linear Segments
  //  The interval is defined by [F0,F1] in frequency and [loge(A0),loge(A1)] in amplitude for
  //      Exponential Segments.

  //  Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
  typedef enum  eSegType
  {
    eExp_RelErr,   // exp(B*f+A) segment, relative error
    eExp_AbsErr,   // exp(B*f+A) segment, absolute error
    eLin_RelErr,   // B*f+A segment, relative error
    eLin_AbsErr    // B*f+A segment, absolute error
  } SegType;

public:

  // Exponential to a imaginary power
  static inline complex<double> ExpI( const double ImagPart )
  {
    return  complex<double>( cos(ImagPart), sin(ImagPart) );
  }

  // The E1 exponential integral
  static inline complex<double> ExpInt_1_I( const double ImagArg )
  {
    return  complex<double>( -gsl_sf_Ci(fabs(ImagArg)), -.5*M_PI+gsl_sf_Si(ImagArg) );
  }

  // The E2 exponential integral
  static inline complex<double> ExpInt_2_I( const double Arg )
  {
    return  ( ExpI(-Arg) - MultI( Arg * ExpInt_1_I( Arg ) ) );
  }

  // Multiply by i
  static inline complex<double> MultI( const complex<double> Arg )
  {
    return  complex<double>( -imag(Arg), real(Arg) );
  }

  //>> START class Seg
  class  Seg
  {
  public:

    // Linear Segments are: B*f+A.  The interval is defined by [F0,F1] in frequency and [A0,A1] in amplitude.
    // Exponential Segments are: exp(B*f+A).  The interval is defined by [F0,F1] in frequency and
    //                                        [loge(A0),loge(A1)] in amplitude.
    //  Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
    SegType  Type;
    double   F0, F1, A0, A1, W;

    Seg( const SegType Type_Param,
         const double  F0_Param, const double  F1_Param,
         const double  A0_Param, const double  A1_Param,
         const double  Weight_Param=1. )
    {
      Type = Type_Param;
      F0   = F0_Param;
      F1   = F1_Param;
      A0   = A0_Param;
      A1   = A1_Param;
      W    = Weight_Param;
    }

    // This call modifies the Q matrix and R vectors to account for this segment.  It in turn calls
    //    specific generator functions based on the type of segment desired.
    //   Linear Segments are: B*f+A.  The interval is defined by [F0,F1] in frequency and [A0,A1] in amplitude.
    //   Exponential Segments are: exp(B*f+A).  The interval is defined by [F0,F1] in frequency and
    //                                          [loge(A0),loge(A1)] in amplitude.
    //  Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
    void   Gen_Segment( vector<complex<double>>&  Q,
                        vector<complex<double>>&  R,
                        const double              Gain=1.,
                        const double              Weight=1. );

    // Display a segment definition
    string   Display( const double SampleRate=1 ) const
    {
      ostringstream  oS;

      oS << "SEG[";

      switch(Type)
      {
      case eExp_RelErr:
        oS << "Exp_Rel";
      break;

      case eExp_AbsErr:
        oS << "Exp_Abs";
      break;

      case eLin_RelErr:
        oS << "Lin_Rel";
      break;

      case eLin_AbsErr:
        oS << "Abs_Lin";
      break;

      default:
        oS << "Error";
        assert(false);
      break;
      }

      oS  << "] "
          << setiosflags(ios::fixed) << setprecision(6) << setw(12) << (this->F0*SampleRate)
          << ": "
          << setiosflags(ios::fixed) << setprecision(8) << setw(12) << this->A0
          << "  <--->  "
          << setiosflags(ios::fixed) << setprecision(6) << setw(12) << (this->F1*SampleRate)
          << ": "
          << setiosflags(ios::fixed) << setprecision(8) << setw(12) << this->A1
          << "  W:"
          << setiosflags(ios::fixed) << setprecision(8) << setw(12) << this->W
          << "]";

      return  oS.str();
    }

protected:

    // exp(B*f+A) segment, relative error
    // -Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
    void  Gen_ExpRel_Segment( vector<std::complex<double>>&  Q,
                              vector<std::complex<double>>&  R,
                              const double                   Gain,
                              const double                   Weight );

    // exp(B*f+A) segment, relative error
    // -Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
    void  Gen_ExpAbs_Segment( vector<std::complex<double>>&  Q,
                              vector<std::complex<double>>&  R,
                              const double                   Gain,
                              const double                   Weight );

    // B*f+A segment, relative error
    // -Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
    void  Gen_LinRel_Segment( vector<std::complex<double>>&  Q,
                              vector<std::complex<double>>&  R,
                              const double                   Gain,
                              const double                   Weight );

    // B*f+A segment, absolute error
    // -Note that frequencies are specified in the interval [0,1).  1/2 is PRF/2 at all sample rates.
    void  Gen_LinAbs_Segment( vector<std::complex<double>>&  Q,
                              vector<std::complex<double>>&  R,
                              const double                   Gain,
                              const double                   Weight );
  };
  //<< END class  Seg

  // Linear system solution of 'matrix(Q) vector(Coefs) = vector(R)' via generalized Levenson's recusion
  bool  Levenson( vector<complex<double>>&  pR,
                  vector<complex<double>>&  g,
                  vector<complex<double>>&  X  );

  complex<double>  InnerProdR( const uint32_t                  Start,
                               const uint32_t                  StopPlus1,
                               const vector<complex<double>>&  InData1,
                               const vector<complex<double>>&  InData2    ) const;

  complex<double>  InnerProdC( const uint32_t                  Start,
                               const uint32_t                  StopPlus1,
                               const vector<complex<double>>&  InData1,
                               const vector<complex<double>>&  InData2    ) const;

private:

  double  Gain;   // Default amplitude (1) unless set by '.SetAmplitude(xx)'

  // The segments list
  list<Seg>  Segs;

public:

  CLSLP(void);

  void SetGain( const double GainParam=1 )
  {
    this->Gain = GainParam;
  }

  double GetGain(void) const
  {
    return  this->Gain;
  }

  // Add a new filter segment
  //   Linear Segments are: B*f+A.  The interval is defined by [F0,F1] in frequency and [A0,A1] in amplitude.
  //   Exponential Segments are: exp(B*f+A).  The interval is defined by [F0,F1] in frequency and
  //                                          [loge(A0),loge(A1)] in amplitude.
  //  Note that frequencies are specified in the interval [0,1).  1/2 is SampleRate/2 here.
  void  Add( const SegType Type,
             const double  F0, const double  F1,
             const double  A0, const double  A1,
             const double  Weight=1. )
  {
    assert( F0 < F1 );

    // Add: Design for a assymetric filter.  This type of filter will process the real and imaginary channels
    //      with different coefficients.  [To generate a Symmetric, real-valued, filter, specifiy desired
    //      filter responses in a conjugate-symmetric pattern about SampleRate/2.]
    // ===  by the design.
     Segs.push_back( Seg( Type, F0, F1, A0*this->Gain, A1*this->Gain, Weight ) );
  }

  // Add a new filter segment symmetric between the positive andnegative frequencies
  //   Linear Segments are: B*f+A.  The interval is defined by [F0,F1] in frequency and [A0,A1] in amplitude.
  //   Exponential Segments are: exp(B*f+A).  The interval is defined by [F0,F1] in frequency and
  //                                          [loge(A0),loge(A1)] in amplitude.
  //  Note that frequencies are specified in the interval [0,1).  1/2 is SampleRate/2 at all sample rates.
  void  AddSymmetric( const SegType Type,
                      const double  F0, const double  F1,
                      const double  A0, const double  A1,
                      const double  Weight=1. )
  {
    assert( F0 <= .5 );
    assert( F0 <= .5 );
    assert( F0 <  F1 );

    // AddSymmetric: Filter response is conjugate-symmetric about PRF/2 and linear-phase response is assumed
    // ============  by the fiter design algorithm.
    Segs.push_back( Seg( Type,   F0,   F1, A0*this->Gain, A1*this->Gain, Weight ) );
    Segs.push_back( Seg( Type, 1-F1, 1-F0, A1*this->Gain, A0*this->Gain, Weight ) );
  }

  // Generate the filter coefficients out into Out.  This involves the solution of the linear system
  //    'matrix(Q) vector(Coefs) = vector(R)'.
  bool  GenFilter( const unsigned int FilterLen, vector<complex<double>>& Out );

  // Get the number of filter segments
  unsigned int  GetNumSegments(void) const
  {
    return  Segs.size();
  }

  // Clear all filter segments
  void  Clear(void)
  {
    Segs.clear();
  }

  // Get segment list reference
  list<Seg>&  GetSegs(void)
  {
    return  this->Segs;
  }
};


