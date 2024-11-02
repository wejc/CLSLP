///////////////////////////////////////////////////////////////////////////////////////////////////
//  FilePath:  CLSLP_20241021_122444/OTHER_CLSLP_EXAMPLES/Filter_6_Complex_Type 1024 pt/CLSLP.cpp
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

#include  <cassert>
#include  <cmath>
#include  <fstream>
#include  <iostream>
#include  "CLSLP.h"

using namespace std;


// Constants
const double  CLSLP::TPI = 2. * M_PI;
const double  CLSLP::PI2 = M_PI * M_PI;


// Generate the filter coefficients into Out
bool  CLSLP::GenFilter( const unsigned int        FilterLen,
                        vector<complex<double>>&  Out       )
{
  // Is zero length
  if( FilterLen == 0 )  return false;

  vector<complex<double>>  Q(FilterLen);
  vector<complex<double>>  R(FilterLen);

  // Over The Frequency Segments
  auto  I = Segs.begin();

  for(  ; I != Segs.end() ; ++I )
  {
    I->Gen_Segment( Q, R );
  }

  // Perform Levenson-Durbin Recursion To Solve The System
  const bool  Solution = LevensonDurbin( Q, R, Out );

  return  Solution;
}


// Fill the Q-matrix and the R-vector for linear system solution to filter coefficients
void  CLSLP::Seg::Gen_Segment( vector<complex<double>>&  Q,
                               vector<complex<double>>&  R,
                               const double              Gain,
                               const double              Weight )
{
  // Pick the frequency interval segment (line) type
  switch(this->Type)
  {
    case eExp_RelErr:
      this->Gen_ExpRel_Segment( Q, R, Gain, Weight );
    break;

    case eExp_AbsErr:
      this->Gen_ExpAbs_Segment( Q, R, Gain, Weight );
    break;

    case eLin_RelErr:
      this->Gen_LinRel_Segment( Q, R, Gain, Weight );
    break;

    case eLin_AbsErr:
      this->Gen_LinAbs_Segment( Q, R, Gain, Weight );
    break;

    default:
      assert(false);   // Set exception in DEBUG mode
    break;
  }
}


// exp(B*f+A) segment, relative error
void  CLSLP::Seg::Gen_ExpRel_Segment( vector<complex<double>>&  Q,
                                      vector<complex<double>>&  R,
                                      const double              Gain,
                                      const double              Weight  )
{
  const double  Min_Segment = 1e-32;   // Need > 0 here
  const double  MidPoint    = .5 * ( Q.size() - 1 );
  const double  L_A0 = Gain * this->A0;
  const double  L_A1 = Gain * this->A1;
  double        lA1  = 0;
  double        lA0  = 0;
  complex<double>  Ret { 0 };

  if( L_A0 > Min_Segment )
  {
    lA0 = log(L_A0);
  }
  else
  {
    lA0 = log(Min_Segment);
  }

  if( L_A1 > Min_Segment )
  {
    lA1 = log(L_A1);
  }
  else
  {
    lA1 = log(Min_Segment);
  }

  // Put in standard form
  const double  sB  = (lA1  - lA0) / (this->F1 - this->F0);
  const double  sA  = lA0 - sB  * this->F0;

  // Form The Q Toeplitz Matrix (Vector Here)
  for( auto L=0U ; L < Q.size() ; ++L )
  {
    if( !((L==0) && (sB == 0)) )
    {
      Ret   = ExpI(TPI*L*this->F1) * exp(-2.*(sB*this->F1+sA));
      Ret  -= ExpI(TPI*L*this->F0) * exp(-2.*(sB*this->F0+sA));
      Q[L] +=  Weight * Ret / complex<double>( -2.*sB, TPI*L );
    }
    else
    {
      Q[L] +=  Weight * exp(-2.*sA) * (this->F1 - this->F0);
    }
  }

  // Form The R Vector
  for( auto L=0U ; L < R.size() ; ++L )
  {
    if( !((L==MidPoint) && (sB == 0)) )
    {
      Ret  = ExpI(-TPI*this->F1*(L-MidPoint)) * exp(-(sB*this->F1+sA));
      Ret -= ExpI(-TPI*this->F0*(L-MidPoint)) * exp(-(sB*this->F0+sA));
      R[R.size()-1-L] +=  Weight * Ret / complex<double>( -sB, -TPI*(L-MidPoint));
    }
    else
    {
      R[R.size()-1-L] +=  Weight * exp(sA) * (this->F1 - this->F0);
    }
  }
}


// exp(B*f+A) segment, absolute error
void  CLSLP::Seg::Gen_ExpAbs_Segment( vector<complex<double>>&  Q,
                                      vector<complex<double>>&  R,
                                      const double              Gain,
                                      const double              Weight  )
{
  const double     Min_Segment = 1e-32;   // Need > 0 here
  const double     MidPoint = .5 * ( Q.size() - 1 );
  const double     L_A0 = Gain * this->A0;
  const double     L_A1 = Gain * this->A1;
  double           lA1  = 0;
  double           lA0  = 0;
  complex<double>  Ret { 0 };

  if( L_A0 > Min_Segment )
  {
    lA0 = log(L_A0);
  }
  else
  {
    lA0 = log(Min_Segment);
  }

  if( L_A1 > Min_Segment )
  {
    lA1 = log(L_A1);
  }
  else
  {
    lA1 = log(Min_Segment);
  }

  // Put in standard form
  const double  sB  = (lA1  - lA0) / (this->F1 - this->F0);
  const double  sA  = lA0 - sB  * this->F0;

  // Form The Q Toeplitz Matrix (Vector Here)
  for( auto L=0U ; L < Q.size() ; ++L )
  {
    // See if we can assume a zero denominator
    if( L != 0 )
    { // L != 0
      Ret  = -ExpI( TPI*L*this->F1 );
      Ret -= -ExpI( TPI*L*this->F0 );

      Q[L]  +=  Weight * MultI( Ret ) / (TPI*L);
    }
    else
    { // L == 0
      Q[L]  +=  Weight * ( this->F1 - this->F0 );
    }
  }

  // Form The R Vector
  for( auto L=0U ; L < R.size() ; ++L )
  {
    if( !((L == MidPoint) && (fabs(sB) <= Min_Segment)) )
    { // Not midpoint
      Ret  =  ExpI(-TPI*this->F1*(L-MidPoint)) * exp(sB*this->F1+sA);
      Ret -=  ExpI(-TPI*this->F0*(L-MidPoint)) * exp(sB*this->F0+sA);

      R[R.size()-1-L] += Weight * Ret / complex<double>( sB, -TPI*(L-MidPoint) );
    }
    else
    { // Midpoint
      R[R.size()-1-L] += Weight * ( this->F1 - this->F0 ) * exp(sA);
    }
  }
}


// B*f+A segment, relative error
void  CLSLP::Seg::Gen_LinRel_Segment( vector<complex<double>>&  Q,
                                      vector<complex<double>>&  R,
                                      const double              Gain,
                                      const double              Weight  )
{
  const double     MidPoint = .5 * ( Q.size() - 1 );
  const double     L_A0 = Gain * this->A0;
  const double     L_A1 = Gain * this->A1;
  double           Phase0 = 0;
  double           Phase1 = 0;
  complex<double>  Ret { 0 };

  // Put in standard form
  const double  sB = (L_A1  - L_A0) / (this->F1 - this->F0);
  const double  sA = L_A0 - sB  * this->F0;

  // Form The Q Toeplitz Matrix (Vector Here)
  for( auto L=0U ; L < Q.size() ; ++L )
  {
    // Handle special cases
    if( L != 0 )
    { // L != 0
      if( sB != 0 )
      { // B != 0  &&  L != 0
        Phase1 = -TPI*L*(sB*this->F1+sA) / sB;
        Phase0 = -TPI*L*(sB*this->F0+sA) / sB;
        Ret    =  -ExpInt_2_I(Phase1) / (sB*(sB*this->F1+sA));
        Ret   -=  -ExpInt_2_I(Phase0) / (sB*(sB*this->F0+sA));
        Q[L]  +=  Weight * Ret * ExpI(-TPI*sA*L/sB);
      }
      else
      { // B == 0  &&  L != 0
        Ret   =  -MultI( ExpI(TPI*L*this->F1) );
        Ret  -=  -MultI( ExpI(TPI*L*this->F0) );
        Q[L] +=  Weight * Ret / (TPI*sA*sA*L);
      }
    }
    else
    { // L == 0  &&  B != 0
      if( sB != 0 )
      { // B!= 0
        Ret   =  -1. / (sB*(sB*this->F1+sA));
        Ret  -=  -1. / (sB*(sB*this->F0+sA));
        Q[L] +=  Weight * Ret;
      }
      else
      { // L == 0  &&  B == 0
        Q[L] +=  Weight * (this->F1 - this->F0) / (sA*sA);
      }
    }
  }

  // Form The R Vector (Vector Here)
  for( auto L=0U ; L < R.size() ; ++L )
  {
    // Handle special cases
    if( L != MidPoint )
    { // L != MidPoint
      if( sB != 0 )
      { // B != 0  &&  L != MidPoint
        Phase1 = TPI * (sB*this->F1+sA)*(L-MidPoint) / sB;
        Phase0 = TPI * (sB*this->F0+sA)*(L-MidPoint) / sB;
        Ret    =  -ExpInt_1_I(Phase1);
        Ret   -=  -ExpInt_1_I(Phase0);
        R[R.size()-1-L]  +=  Weight * Ret * ExpI(TPI*sA*(L-MidPoint)/sB) / sB;
      }
      else
      { // B == 0  &&  L != MidPoint
        Ret    =  MultI( ExpI(-TPI*this->F1*(L-MidPoint)));
        Ret   -=  MultI( ExpI(-TPI*this->F0*(L-MidPoint)));
        R[R.size()-1-L]  +=  Weight * Ret / (TPI*sA*(L-MidPoint));
      }
    }
    else
    { // L == MidPoint  &&  B != 0
      if( sB != 0 )
      { // B != 0
        Ret  =  log(sB*this->F1+sA);
        Ret -=  log(sB*this->F0+sA);
        R[R.size()-1-L]  +=  Weight * Ret / sB;
      }
      else
      { // L == MidPoint  &&  B == 0
        R[R.size()-1-L]  +=  Weight * (this->F1 - this->F0) / sA;
      }
    }
  }
}


// B*f+A segment, absolute error
void  CLSLP::Seg::Gen_LinAbs_Segment( vector<complex<double>>&  Q,
                                      vector<complex<double>>&  R,
                                      const double              Gain,
                                      const double              Weight  )
{
  const double  MinQR_Denom = 1e-6;
  const double  MidPoint = .5 * ( Q.size() - 1 );
  const double  L_A0   = Gain * this->A0;
  const double  L_A1   = Gain * this->A1;
        double  Phase0 = 0;
        double  Phase1 = 0;
        double  Denom  = 0;

  // Put in standard form
  const double  sB  = (L_A1 - L_A0) / (this->F1 - this->F0);
  const double  sA  = L_A0 - sB * this->F0;

  // Form The Q Toeplitz Matrix (Vector Here)
  for( auto L=0U ; L < Q.size() ; ++L )
  {
    Denom = -TPI * L;

    // See if we can assume a zero denominator
    if( abs(Denom) >= MinQR_Denom )
    { // Denominator not zero
      Phase0 = TPI * this->F0 * L;
      Phase1 = TPI * this->F1 * L;

      // Update Q
      Q[L] += Weight * ( complex<double>( cos(Phase1), sin(Phase1) ) -
                         complex<double>( cos(Phase0), sin(Phase0) ) ) /
                         complex<double>( 0.,  -Denom  );
    }
    else
    { // Denominator zero
      // Update Q
      Q[L] += Weight * (this->F1 - this->F0);
    }
  }

  // Form The R Vector
  for( auto L=0U ; L < R.size() ; ++L )
  {
    Denom = 4. * PI2 * pow( (double) L - MidPoint, 2 );

    if( abs(Denom) >= MinQR_Denom )
    { // Denominator not zero
      // Update R
      const complex<double>  A1( 0.,  L*TPI*(sB*this->F1 + sA) );
      const complex<double>  A0( 0.,  L*TPI*(sB*this->F0 + sA) );
      const complex<double>  B1( sB, -TPI*MidPoint*(sB*this->F1 + sA) );
      const complex<double>  B0( sB, -TPI*MidPoint*(sB*this->F0 + sA) );

      Phase0 = -TPI * F0 * ( L - MidPoint );
      Phase1 = -TPI * F1 * ( L - MidPoint );

      const complex<double>  X= ((A1 + B1) * complex<double>( cos(Phase1), sin(Phase1) )) -
                                ((A0 + B0) * complex<double>( cos(Phase0), sin(Phase0) )) ;

      // Update R
      R[R.size()-1-L] += Weight * X / Denom;
   }
    else
    { // Denominator zero
      // Update R
      R[R.size()-1-L] += Weight * ( ( .5*sB*this->F1*this->F1 + sA*this->F1 ) -
                                    ( .5*sB*this->F0*this->F0 + sA*this->F0 ) );
   }
  }
}


// Must Have At Least 2 Dimensions
bool  CLSLP::LevensonDurbin( const vector<complex<double>>& Tm_in,
                             const vector<complex<double>>& Yv_in,
                                   vector<complex<double>>& Xv_out  )
{
  vector<complex<double>>  F(Xv_out.size());
  vector<complex<double>>  Fnew(Xv_out.size());
  complex<double>          Csum=0, Delta=0, Alpha=0, Xcorr=0, Ex=0;
  const double             MinDivisor = 1e-12;
  unsigned int             K = 0;

  // Size and data restrictions
  if( Tm_in.size()  < 2 )  return  false;
  if( Yv_in.size()  < 2 )  return  false;
  if( Xv_out.size() < 2 )  return  false;
  if( norm(Tm_in[0]) < (MinDivisor*MinDivisor) )  return  false;

  // Step 0
  Xv_out[0] = Yv_in[0]/Tm_in[0];
  F[0]      = -Tm_in[1]/Tm_in[0];

  for( auto Step=1U ; Step < Xv_out.size() ; ++Step )
  {
    // Calculate Delta
    for( K=0U, Csum=0 ; K < Step ; ++K )
    {
      Csum += Tm_in[K+1] * conj(F[K]);
    }

    Delta = Csum + Tm_in[0];

    if( norm(Delta) < (MinDivisor*MinDivisor) )  return  false;

    // Calculate Ex
    for( K=0U, Ex=0 ; K < Step ; ++K )
    {
      Ex += Tm_in[K+1] * Xv_out[Step-1-K];
    }

    Xcorr = ( Yv_in[Step] - Ex ) / Delta;

    // Update Xv_out
    for( auto L=0U ; L<Step ; L++ )
    {
      Xv_out[L] += Xcorr * conj(F[Step-1-L]);
    }

    // Output new X
    Xv_out[Step] = Xcorr;

    if( (Step+1) < Xv_out.size() )
    {
      for( K=0U, Csum=0 ; K < Step ; ++K )
      {
        Csum += Tm_in[K+1] * F[Step-1-K];
      }

      // Calculate Alpha=1/(1-Eb Ef*)
      Alpha = -( Tm_in[Step+1] + Csum ) / Delta;

      for( auto L=0U ; L<Step ; ++L )
      {
        Fnew[L] = F[L] + Alpha * conj(F[Step-1-L]);
      }

      // Fnew ---> F
      std::copy( &Fnew[0], &Fnew[Step], &F[0] );

      // New last F-Vector data
      F[Step] = Alpha;
    }
  }

  // Good return
  return  true;
}


