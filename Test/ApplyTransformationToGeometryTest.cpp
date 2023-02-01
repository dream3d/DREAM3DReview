/* ============================================================================
 * Copyright (c) 2019-2019 BlueQuartz Software, LLC
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
 * Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
 * contributors may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-15-D-5231
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/*
TEST ONLY CHECKS FEATURES ADDED FOR IMAGE TRANSFORMATION BASED ON RotateSampleRefFrameTest
*/
#include <QtCore/QFile>

#include <Eigen/Dense>

#include <iostream>

class ApplyTransformationToGeometryTest
{
private:
public:
  ApplyTransformationToGeometryTest() = default;
  ~ApplyTransformationToGeometryTest() = default;
  ApplyTransformationToGeometryTest(const ApplyTransformationToGeometryTest&) = delete;            // Copy Constructor Not Implemented
  ApplyTransformationToGeometryTest(ApplyTransformationToGeometryTest&&) = delete;                 // Move Constructor Not Implemented
  ApplyTransformationToGeometryTest& operator=(const ApplyTransformationToGeometryTest&) = delete; // Copy Assignment Not Implemented
  ApplyTransformationToGeometryTest& operator=(ApplyTransformationToGeometryTest&&) = delete;      // Move Assignment Not Implemented

  // -----------------------------------------------------------------------------
  //
  // -----------------------------------------------------------------------------
  void operator()()
  {
    std::cout << "----Start ApplyTransformationToGeometryTest----\n";

    int err = EXIT_SUCCESS;
    // clang-format off
    using Matrix4f = Eigen::Matrix<float, 4, 4, Eigen::RowMajor>;
    Matrix4f A;
     A << 5.0F, 0.0F, 0.0F, 0.0F,
          0.0F, 1.0F, 0.0F, 0.0F,
          0.0F, 0.0F, 1.0F, 0.0F,
          0.0F, 0.0F, 0.0F, 1.0F;
    Matrix4f B;
    B << 0.866025403784439F, -0.5, 0.0F, 0.0F,
    0.5, 0.866025403784439F, 0.0F, 0.0F,
    0.0F, 0.0F, 1.0F, 0.0F,
    0.0F, 0.0F, 0.0F, 1.0F;
    
    
    Matrix4f C;
     C << 1.0F, 0.0F, 0.0F, 2.0F,
          0.0F, 1.0F, 0.0F, 2.0F,
          0.0F, 0.0F, 1.0F, 0.0F,
          0.0F, 0.0F, 0.0F, 1.0F;


  auto transformationMatrix = C * B * A;
std::cout << transformationMatrix << std::endl;

/*
 4.33013     -0.5        0        2
     2.5 0.866025        0        2
       0        0        1        0
       0        0        0        1
*/

    std::cout << "----End ApplyTransformationToGeometryTest----\n";





  }
};
