/* ============================================================================
 * Copyright 2021 The University of Utah
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of BlueQuartz Software, the US Air Force, the University of Utah nor the names of its contributors may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the following contracts:
 *
 *
 * This code contained herein is based upon work supported by the following grants:
 *    DOE Office of Nuclear Energy's Nuclear Energy University Program Grant No.: DE-NE0008799
 *    DOD Office of Economic Adjustment Grant No.: ST1605-19-03
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include "SIMPLib/Math/SIMPLibMath.h"

namespace EigenstrainsHelper
{
namespace SIMPLMath = SIMPLib::Constants;

/**
 *@brief Small structure/class to wrap the 81 elements of a 4D tensor and access them via i,j,k,l indices
 */
template <class T, size_t I, size_t J, size_t K, size_t L>
struct Tensor4D
{
  static_assert(I > 0 && J > 0 && K > 0 && L > 0);
  static constexpr size_t Size = I * J * K * L;
  T& operator()(size_t i, size_t j, size_t k, size_t l) noexcept
  {
    size_t idx = i * J * K * L + j * K * L + k * L + l;
    return data[idx];
  }
  constexpr const T& operator()(size_t i, size_t j, size_t k, size_t l) const noexcept
  {
    size_t idx = i * J * K * L + j * K * L + k * L + l;
    return data[idx];
  }
  void init(T value)
  {
    data.fill(value);
  }
  std::array<T, Size> data;
};

using Tensor4DType = Tensor4D<double, 3, 3, 3, 3>;

// -----------------------------------------------------------------------------
// Calculates 32-point Gaussian quadrature of the input function
// -----------------------------------------------------------------------------
template <typename Functor>
double gauss_integration(Functor function, double lowerBound, double upperBound)
{
  // Hard coded for 32-point integration
  // clang-format off
  double points[32] = { 0.0483076656877383162348126,  0.1444719615827964934851864,  0.2392873622521370745446032,  0.3318686022821276497799168,
                        0.4213512761306353453641194,  0.5068999089322293900237475,  0.5877157572407623290407455,  0.6630442669302152009751152,
                        0.7321821187402896803874267,  0.7944837959679424069630973,  0.8493676137325699701336930,  0.8963211557660521239653072,
                        0.9349060759377396891709191,  0.9647622555875064307738119,  0.9856115115452683354001750,  0.9972638618494815635449811,
                       -0.0483076656877383162348126, -0.1444719615827964934851864, -0.2392873622521370745446032, -0.3318686022821276497799168,
                       -0.4213512761306353453641194, -0.5068999089322293900237475, -0.5877157572407623290407455, -0.6630442669302152009751152,
                       -0.7321821187402896803874267, -0.7944837959679424069630973, -0.8493676137325699701336930, -0.8963211557660521239653072,
                       -0.9349060759377396891709191, -0.9647622555875064307738119, -0.9856115115452683354001750, -0.9972638618494815635449811};

  double weights[32] = {0.0965400885147278005667648,  0.0956387200792748594190820,  0.0938443990808045656391802,  0.0911738786957638847128686,
                        0.0876520930044038111427715,  0.0833119242269467552221991,  0.0781938957870703064717409,  0.0723457941088485062253994,
                        0.0658222227763618468376501,  0.0586840934785355471452836,  0.0509980592623761761961632,  0.0428358980222266806568786,
                        0.0342738629130214331026877,  0.0253920653092620594557526,  0.0162743947309056706051706,  0.0070186100094700966004071,
                        0.0965400885147278005667648,  0.0956387200792748594190820,  0.0938443990808045656391802,  0.0911738786957638847128686,
                        0.0876520930044038111427715,  0.0833119242269467552221991,  0.0781938957870703064717409,  0.0723457941088485062253994,
                        0.0658222227763618468376501,  0.0586840934785355471452836,  0.0509980592623761761961632,  0.0428358980222266806568786,
                        0.0342738629130214331026877,  0.0253920653092620594557526,  0.0162743947309056706051706,  0.0070186100094700966004071};
  // clang-format on

  double evaluation = 0;
  double sum = 0;
  double range = (upperBound - lowerBound) / 2;
  double mid = (upperBound + lowerBound) / 2;

  for(size_t i = 0; i < 32; i++)
  {
    sum += weights[i] * function(range * points[i] + mid);
  }

  evaluation = range * sum;
  return evaluation;
}

// -----------------------------------------------------------------------------
// Calculates the fourth-rank Eshelby tensor using Poisson's ratio and ellipsoid information
// Equations referenced come from Chap. 2 in "Micromechanics of Defects in Solids" 2nd Ed. by Toshio Mura, 1987
// -----------------------------------------------------------------------------
Tensor4DType find_eshelby(double a, double b, double c, double nu, bool ellipsoidal)
{
  Tensor4DType eshelbyTensor;
  eshelbyTensor.init(0.0);

  // Ellipsoidal solution criteria
  if((a - b > 1e-5 || b - c > 1e-5) && ellipsoidal)
  {
    double IVector[3] = {0};
    double IArray[3][3] = {0};
    double aa = std::pow(a, 2);
    double bb = std::pow(b, 2);
    double cc = std::pow(c, 2);
    double axesSq[3] = {aa, bb, cc};

    if(a - b < 1e-5 && b - c > 1e-5) // Oblate spheroid a = b > c (Eq. 11.28)
    {
      IVector[0] = IVector[1] = ((2 * SIMPLMath::k_PiD * aa * c) / std::pow((aa - cc), 1.5)) * (std::acos(c / a) - (c / a) * std::sqrt(1 - cc / aa));
      IVector[2] = 4 * SIMPLMath::k_PiD - 2 * IVector[0];

      IArray[0][2] = IArray[2][0] = IArray[1][2] = IArray[2][1] = (IVector[0] - IVector[2]) / (cc - aa);
      IArray[0][0] = IArray[1][1] = IArray[0][1] = IArray[1][0] = SIMPLMath::k_PiD / aa - IArray[0][2] / 4;
      IArray[2][2] = ((4 * SIMPLMath::k_PiD) / cc - 2 * IArray[0][2]) / 3;
    }
    else if(a - b > 1e-5 && b - c < 1e-5) // Prolate spheroid a > b = c (Eq. 11.29)
    {
      IVector[1] = IVector[2] = ((2 * SIMPLMath::k_PiD * a * cc) / std::pow((aa - cc), 1.5)) * ((a / c) * std::sqrt(aa / cc - 1) - std::acosh(a / c));
      IVector[0] = 4 * SIMPLMath::k_PiD - 2 * IVector[1];

      IArray[0][1] = IArray[1][0] = IArray[0][2] = IArray[2][0] = (IVector[1] - IVector[0]) / (aa - bb);
      IArray[1][1] = IArray[2][2] = IArray[1][2] = IArray[2][1] = SIMPLMath::k_PiD / bb - IArray[0][1] / 4;
      IArray[0][0] = ((4 * SIMPLMath::k_PiD) / aa - 2 * IArray[0][1]) / 3;
    }
    else // Ellipsoid
    {
      // Functional form of elliptic integrals of first and second kind (Eq. 11.18)
      double theta = std::asin(std::sqrt(1 - cc / aa));
      double k = std::sqrt((aa - bb) / (aa - cc));
      auto F = [k](double w) { return 1 / std::sqrt(1 - std::pow(k, 2) * std::pow(std::sin(w), 2)); };
      auto E = [k](double w) { return std::sqrt(1 - std::pow(k, 2) * std::pow(std::sin(w), 2)); };

      // Calculate elliptic integrals w/ 32-point Gauss quadrature integration
      double FIntegral = gauss_integration(F, 0, theta);
      double EIntegral = gauss_integration(E, 0, theta);

      // This would be the boost implementation of 30-point integration
      // boost::math::quadrature::gauss<double, 30> integrator;
      // double FIntegral = integrator.integrate(F, 0, theta);
      // double EIntegral = integrator.integrate(E, 0, theta);

      // Solve I1, I2, I3 (Eq. 11.17 and Formula 1 in Eq. 11.19)
      IVector[0] = ((4.0 * SIMPLMath::k_PiD * a * b * c) / ((aa - bb) * std::sqrt(aa - cc))) * (FIntegral - EIntegral);
      IVector[2] = ((4.0 * SIMPLMath::k_PiD * a * b * c) / ((bb - cc) * std::sqrt(aa - cc))) * ((b * std::sqrt(aa - cc)) / (a * c) - EIntegral);
      IVector[1] = 4.0 * SIMPLMath::k_PiD - IVector[0] - IVector[2];

      // Solve for I off-diagonal terms (Formula 4 in Eq. 11.19)
      for(size_t m = 0; m < 3; m++)
      {
        for(size_t n = 0; n < 3; n++)
        {
          if(m != n)
          {
            IArray[m][n] = (IVector[n] - IVector[m]) / (axesSq[m] - axesSq[n]);
          }
        }
      }

      // Solve for I diagonal terms (Formula 2 in Eq. 11.19)
      IArray[0][0] = (4.0 * SIMPLMath::k_PiD / axesSq[0] - IArray[0][1] - IArray[0][2]) / 3;
      IArray[1][1] = (4.0 * SIMPLMath::k_PiD / axesSq[1] - IArray[1][0] - IArray[1][2]) / 3;
      IArray[2][2] = (4.0 * SIMPLMath::k_PiD / axesSq[2] - IArray[2][0] - IArray[2][1]) / 3;
    }

    // Ellipsoid cyclic permutation (Eq. 11.16)
    for(size_t i = 0; i < 3; i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        for(size_t k = 0; k < 3; k++)
        {
          for(size_t l = 0; l < 3; l++)
          {
            if(i == j && k == l)
            {
              if(j == k)
              {
                // Case for first diagonal terms (S11-S33)
                eshelbyTensor(i, j, k, l) = 3.0 / (8.0 * SIMPLMath::k_PiD * (1.0 - nu)) * axesSq[k] * IArray[i][k] + ((1.0 - 2.0 * nu) / (8.0 * SIMPLMath::k_PiD * (1.0 - nu))) * IVector[i];
              }
              else
              {
                // Case for off-diagonal (S12, S21, S13, S31, S23, S32)
                eshelbyTensor(i, j, k, l) = 1.0 / (8.0 * SIMPLMath::k_PiD * (1.0 - nu)) * axesSq[k] * IArray[i][k] - ((1.0 - 2.0 * nu) / (8.0 * SIMPLMath::k_PiD * (1.0 - nu))) * IVector[i];
              }
            }
            else if((i == k && j == l) || (i == l && j == k))
            {
              // Case for second diagonal terms (S44-S66)
              eshelbyTensor(i, j, k, l) =
                  (axesSq[i] + axesSq[j]) / (16.0 * SIMPLMath::k_PiD * (1.0 - nu)) * IArray[i][j] + ((1.0 - 2.0 * nu) / (16.0 * SIMPLMath::k_PiD * (1.0 - nu))) * (IVector[i] + IVector[j]);
            }
            else
            {
              eshelbyTensor(i, j, k, l) = 0;
            }
          }
        }
      }
    }
  }
  else
  {
    // Spherical cyclic permutation (Eq. 11.21)
    for(size_t i = 0; i < 3; i++)
    {
      for(size_t j = 0; j < 3; j++)
      {
        for(size_t k = 0; k < 3; k++)
        {
          for(size_t l = 0; l < 3; l++)
          {
            if(i == j && k == l)
            {
              if(j == k)
              {
                // Case for first diagonal terms (S11-S33)
                eshelbyTensor(i, j, k, l) = ((7.0 - 5.0 * nu) / (15.0 * (1.0 - nu)));
              }
              else
              {
                // Case for off-diagonal (S12, S21, S13, S31, S23, S32)
                eshelbyTensor(i, j, k, l) = ((5.0 * nu - 1.0) / (15.0 * (1.0 - nu)));
              }
            }
            else if((i == k && j == l) || (i == l && j == k))
            {
              // Case for second diagonal terms (S44-S66)
              eshelbyTensor(i, j, k, l) = ((4.0 - 5.0 * nu) / (15.0 * (1.0 - nu)));
            }
            else
            {
              eshelbyTensor(i, j, k, l) = 0;
            }
          }
        }
      }
    }
  }
  return eshelbyTensor;
}
} // namespace EigenstrainsHelper
