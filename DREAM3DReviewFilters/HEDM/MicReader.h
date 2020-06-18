/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2010, Dr. Michael A. Groeber (US Air Force Research Laboratories
 * All rights reserved.
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include <string>

#include "EbsdLib/EbsdLib.h"
#include "EbsdLib/IO/EbsdReader.h"
#include "EbsdLib/Core/EbsdSetGetMacros.h"
#include "EbsdLib/Core/EbsdLibConstants.h"

#include "MicConstants.h"
#include "MicHeaderEntry.h"
#include "MicPhase.h"

/**
 * @class MicReader MicReader.h EbsdLib/HEDM/MicReader.h
 * @brief This class is a self contained HEDM OIM .Mic file reader and will read a
 * single .Mic file and store all the data in column centric pointers.
 * @author Michael A. Jackson for BlueQuartz Software
 * @date Mar 1, 2010
 * @version 1.0
 */
class MicReader : public EbsdReader
{
public:
  MicReader();
  ~MicReader() override;

  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, InfileBasename, Mic::InfileBasename)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, InfileSerialLength, Mic::InfileSerialLength)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, OutfileBasename, Mic::OutfileBasename)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, OutfileSerialLength, Mic::OutfileSerialLength)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, OutStructureBasename, Mic::OutStructureBasename)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, BCPeakDetectorOffset, Mic::BCPeakDetectorOffset)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, InFileType, Mic::InFileType)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, OutfileExtension, Mic::OutfileExtension)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, InfileExtesnion, Mic::InfileExtesnion)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, BeamEnergyWidth, Mic::BeamEnergyWidth)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, BeamDirection, Mic::BeamDirection)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, BeamDeflectionChiLaue, Mic::BeamDeflectionChiLaue)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, BeamHeight, Mic::BeamHeight)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, BeamEnergy, Mic::BeamEnergy)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, DetectorFilename, Mic::DetectorFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, OptimizationConstrainFilename, Mic::OptimizationConstrainFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, EtaLimit, Mic::EtaLimit)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, SampleFilename, Mic::SampleFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, StructureFilename, Mic::StructureFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, RotationRangeFilename, Mic::RotationRangeFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, FundamentalZoneFilename, Mic::FundamentalZoneFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, SampleSymmetry, Mic::SampleSymmetry)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MinAmplitudeFraction, Mic::MinAmplitudeFraction)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, MaxQ, Mic::MaxQ)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MaxInitSideLength, Mic::MaxInitSideLength)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MinSideLength, Mic::MinSideLength)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, LocalOrientationGridRadius, Mic::LocalOrientationGridRadius)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MinLocalResolution, Mic::MinLocalResolution)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MaxLocalResolution, Mic::MaxLocalResolution)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MaxAcceptedCost, Mic::MaxAcceptedCost)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MaxConvergenceCost, Mic::MaxConvergenceCost)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, MaxMCSteps, Mic::MaxMCSteps)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MCRadiusScaleFactor, Mic::MCRadiusScaleFactor)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, SuccessiveRestarts, Mic::SuccessiveRestarts)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, SecondsBetweenSave, Mic::SecondsBetweenSave)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, NumParameterOptimizationSteps, Mic::NumParameterOptimizationSteps)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, NumElementToOptimizePerPE, Mic::NumElementToOptimizePerPE)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, OptimizationFilename, Mic::OptimizationFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, DetectionLimitFilename, Mic::DetectionLimitFilename)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, ParameterMCInitTemperature, Mic::ParameterMCInitTemperature)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, OrientationSearchMethod, Mic::OrientationSearchMethod)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, CoolingFraction, Mic::CoolingFraction)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, ThermalizeFraction, Mic::ThermalizeFraction)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, ParameterRefinements, Mic::ParameterRefinements)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, NumDetectors, Mic::NumDetectors)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, DetectorSpacing, Mic::DetectorSpacing)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, DetectorSpacingDeviation, Mic::DetectorSpacingDeviation)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, DetectorOrientationDeviationInEuler, Mic::DetectorOrientationDeviationInEuler)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, DetectorOrientationDeviationInSO3, Mic::DetectorOrientationDeviationInSO3)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, ParamMCMaxLocalRestarts, Mic::ParamMCMaxLocalRestarts)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, ParamMCMaxGlobalRestarts, Mic::ParamMCMaxGlobalRestarts)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, ParamMCNumGlobalSearchElements, Mic::ParamMCNumGlobalSearchElements)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, ConstrainedOptimization, Mic::ConstrainedOptimization)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, SearchVolumeReductionFactor, Mic::SearchVolumeReductionFactor)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, FileNumStart, Mic::FileNumStart)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, FileNumEnd, Mic::FileNumEnd)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, SampleLocation, Mic::SampleLocation)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, SampleOrientation, Mic::SampleOrientation)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, EnableStrain, Mic::EnableStrain)
  EBSDHEADER_INSTANCE_PROPERTY(MicStringHeaderEntry, std::string, SampleCenter, Mic::SampleCenter)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, SampleRadius, Mic::SampleRadius)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MaxDeepeningHitRatio, Mic::MaxDeepeningHitRatio)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, ConsistencyError, Mic::ConsistencyError)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, BraggFilterTolerance, Mic::BraggFilterTolerance)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, MinAccelerationThreshold, Mic::MinAccelerationThreshold)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, MaxDiscreteCandidates, Mic::MaxDiscreteCandidates)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, XDim, Mic::XDim)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<int>, int, YDim, Mic::YDim)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, XRes, Mic::XRes)
  EBSDHEADER_INSTANCE_PROPERTY(MicHeaderEntry<float>, float, YRes, Mic::YRes)

  EBSD_INSTANCE_PROPERTY(QVector<MicPhase::Pointer>, PhaseVector)

  EBSD_POINTER_PROPERTY(Euler1, Euler1, float)
  EBSD_POINTER_PROPERTY(Euler2, Euler2, float)
  EBSD_POINTER_PROPERTY(Euler3, Euler3, float)
  EBSD_POINTER_PROPERTY(X, X, float)
  EBSD_POINTER_PROPERTY(Y, Y, float)
  EBSD_POINTER_PROPERTY(Confidence, Conf, float)
  EBSD_POINTER_PROPERTY(Phase, Phase, int)
  EBSD_POINTER_PROPERTY(Level, Level, int)
  EBSD_POINTER_PROPERTY(Up, Up, int)

  /**
   * @brief Returns the pointer to the data for a given feature
   * @param featureName The name of the feature to return the pointer to.
   */
  void* getPointerByName(const std::string& featureName) override;

  /**
   * @brief Returns an enumeration value that depicts the numerical
   * primitive type that the data is stored as (Int, Float, etc).
   * @param featureName The name of the feature.
   */
  EbsdLib::NumericTypes::Type getPointerType(const std::string& featureName) override;

  /**
   * @brief Reads the complete HEDM .Mic file.
   * @return 1 on success
   */
  int readFile() override;

  /**
   * @brief Reads ONLY the header portion of the HEDM .Mic file
   * @return 1 on success
   */
  int readHeaderOnly() override;

  int getXDimension() override;
  int getYDimension() override;
  void setXDimension(int xD) override;
  void setYDimension(int yD) override;
  float getXStep();
  float getYStep();

protected:
  /** @brief 'free's the allocated memory and sets the pointer to nullptr
   */
  void deletePointers();

private:
  MicPhase::Pointer m_CurrentPhase;

  int xDim;
  int yDim;
  float xRes;
  float yRes;

  int readMicFile();

  int readDatFile();

  /** @brief Parses the value from a single line of the header section of the HEDM .Mic file
   * @param line The line to parse
   */
  void parseHeaderLine(std::string& line);

  /** @brief Parses the data from a line of data from the HEDM .Mic file
   * @param line The line of data to parse
   */
  void parseDataLine(std::string& line, size_t i);

  /**
   * @brief initPointers
   * @param numElements
   */
  void initPointers(size_t numElements);

public:
  MicReader(const MicReader&) = delete;            // Copy Constructor Not Implemented
  MicReader(MicReader&&) = delete;                 // Move Constructor Not Implemented
  MicReader& operator=(const MicReader&) = delete; // Copy Assignment Not Implemented
  MicReader& operator=(MicReader&&) = delete;      // Move Assignment Not Implemented
};
