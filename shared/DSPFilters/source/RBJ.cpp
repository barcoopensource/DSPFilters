/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#include "DspFilters/Common.h"
#include "DspFilters/RBJ.h"

namespace Dsp {

namespace RBJ {

void LowPass::setup (double sampleRate,
                     double cutoffFrequency,
                     double q)
{
  double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  double AL = sn / (2 * q);
  double b0 = (1 - cs) / 2;
  double b1 =  1 - cs;
  double b2 = (1 - cs) / 2;
  double a0 =  1 + AL;
  double a1 = -2 * cs;
  double a2 =  1 - AL;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void HighPass::setup (double sampleRate,
                      double cutoffFrequency,
                      double q)
{
  double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  double AL = sn / ( 2 * q );
  double b0 =  (1 + cs) / 2;
  double b1 = -(1 + cs);
  double b2 =  (1 + cs) / 2;
  double a0 =  1 + AL;
  double a1 = -2 * cs;
  double a2 =  1 - AL;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void BandPass1::setup (double sampleRate,
                       double centerFrequency,
                       double bandWidth)
{
  double w0 = 2 * doublePi * centerFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  double AL = sn / ( 2 * bandWidth );
  // There was a bug in BandPass2::setup. If you run into any issue do check if bandwidth is
  // interpreted correctly, maybe we have to change this like for BandPass2 and BandStop
  double b0 = bandWidth * AL;// sn / 2;
  double b1 =  0;
  double b2 = -bandWidth * AL;//-sn / 2;
  double a0 =  1 + AL;
  double a1 = -2 * cs;
  double a2 =  1 - AL;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void BandPass2::setup (double sampleRate,
                       double centerFrequency,
                       double bandWidth)
{
  double w0 = 2 * doublePi * centerFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  // double AL = sn / ( 2 * bandWidth ); ->this was the original calculation, but it is
  // misleading since bandWidth is interpreted as quality
  double AL = sn * sinh( doubleLn2/2 * bandWidth * w0/sn );
  double b0 =  AL;
  double b1 =  0;
  double b2 = -AL;
  double a0 =  1 + AL;
  double a1 = -2 * cs;
  double a2 =  1 - AL;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void BandStop::setup (double sampleRate,
                      double centerFrequency,
                      double bandWidth)
{
  double w0 = 2 * doublePi * centerFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  // double AL = sn / ( 2 * bandWidth ); ->this was the original calculation, but it is
  // misleading since bandWidth is interpreted as quality
  double AL = sn * sinh( doubleLn2/2 * bandWidth * w0/sn );
  double b0 =  1;
  double b1 = -2 * cs;
  double b2 =  1;
  double a0 =  1 + AL;
  double a1 = -2 * cs;
  double a2 =  1 - AL;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void LowShelf::setup (double sampleRate,
                      double cutoffFrequency,
                      double gainDb,
                      double shelfSlope)
{
  double A  = pow (10, gainDb/40);
  double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  double AL = sn / 2 * ::std::sqrt ((A + 1/A) * (1/shelfSlope - 1) + 2);
  double sq = 2 * sqrt(A) * AL;
  double b0 =    A*( (A+1) - (A-1)*cs + sq );
  double b1 =  2*A*( (A-1) - (A+1)*cs );
  double b2 =    A*( (A+1) - (A-1)*cs - sq );
  double a0 =        (A+1) + (A-1)*cs + sq;
  double a1 =   -2*( (A-1) + (A+1)*cs );
  double a2 =        (A+1) + (A-1)*cs - sq;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void HighShelf::setup (double sampleRate,
                       double cutoffFrequency,
                       double gainDb,
                       double shelfSlope)
{
  double A  = pow (10, gainDb/40);
  double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  double AL = sn / 2 * ::std::sqrt ((A + 1/A) * (1/shelfSlope - 1) + 2);
  double sq = 2 * sqrt(A) * AL;
  double b0 =    A*( (A+1) + (A-1)*cs + sq );
  double b1 = -2*A*( (A-1) + (A+1)*cs );
  double b2 =    A*( (A+1) + (A-1)*cs - sq );
  double a0 =        (A+1) - (A-1)*cs + sq;
  double a1 =    2*( (A-1) - (A+1)*cs );
  double a2 =        (A+1) - (A-1)*cs - sq;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

// HighShelf2: Anymix addition (based on animix svn r11957 "anymix: using smoothed dsp filter for DDEQ" by MARDO)
void HighShelf2::setup (double sampleRate,
                        double cutoffFrequency,
                        double gainLin,
                        double Q)
{
  const double A = std::max(0.0, gainLin);
  const double aminus1 = A - 1.0;
  const double aplus1 = A + 1.0;
  const double omega = (doublePi * 2.0 * std::max(cutoffFrequency, 2.0)) / sampleRate;
  const double coso = std::cos(omega);
  const double beta = std::sin(omega) * std::sqrt(A) / Q;
  const double aminus1TimesCoso = aminus1 * coso;

  setCoefficients (aplus1 - aminus1TimesCoso + beta,
                   2.0 * (aminus1 - aplus1 * coso),
                   aplus1 - aminus1TimesCoso - beta,
                   A * (aplus1 + aminus1TimesCoso + beta),
                   A * -2.0 * (aminus1 + aplus1 * coso),
                   A * (aplus1 + aminus1TimesCoso - beta));
}

void BandShelf::setup (double sampleRate,
                       double centerFrequency,
                       double gainDb,
                       double bandWidth)
{
  double A  = pow (10, gainDb/40);
  double w0 = 2 * doublePi * centerFrequency / sampleRate;
  double cs = cos(w0);
  double sn = sin(w0);
  double AL = sn * sinh( doubleLn2/2 * bandWidth * w0/sn );
  assert (!Dsp::is_nan (AL));
  double b0 =  1 + AL * A;
  double b1 = -2 * cs;
  double b2 =  1 - AL * A;
  double a0 =  1 + AL / A;
  double a1 = -2 * cs;
  double a2 =  1 - AL / A;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

void AllPass::setup (double sampleRate,
                     double phaseFrequency,
                     double q)
{
  double w0 = 2 * doublePi * phaseFrequency / sampleRate;
  double cs = cos (w0);
  double sn = sin (w0);
  double AL = sn / ( 2 * q );
  double b0 =  1 - AL;
  double b1 = -2 * cs;
  double b2 =  1 + AL;
  double a0 =  1 + AL;
  double a1 = -2 * cs;
  double a2 =  1 - AL;
  setCoefficients (a0, a1, a2, b0, b1, b2);
}

}

}
