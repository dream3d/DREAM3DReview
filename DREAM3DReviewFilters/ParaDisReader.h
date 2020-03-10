/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
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

#include <memory>

#include <QtCore/QFile>
#include <QtCore/QString>
#include <vector>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/CoreFilters/FileReader.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/MeshStructs.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class ParaDisReader ParaDisReader.h DREAM3DLib/IO/ParaDisReader.h
 * @brief
 * @author mjackson
 * @date Sep 28, 2011
 * @version $Revision$
 */
class DREAM3DReview_EXPORT ParaDisReader : public FileReader
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(ParaDisReader SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(ParaDisReader)
  PYB11_FILTER_NEW_MACRO(ParaDisReader)
  PYB11_FILTER_PARAMETER(DataArrayPath, EdgeDataContainerName)
  PYB11_FILTER_PARAMETER(QString, VertexAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, EdgeAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, InputFile)
  PYB11_FILTER_PARAMETER(float, BurgersVector)
  PYB11_FILTER_PARAMETER(QString, NumberOfArmsArrayName)
  PYB11_FILTER_PARAMETER(QString, NodeConstraintsArrayName)
  PYB11_FILTER_PARAMETER(QString, BurgersVectorsArrayName)
  PYB11_FILTER_PARAMETER(QString, SlipPlaneNormalsArrayName)
  PYB11_FILTER_PARAMETER(QString, DomainBoundsArrayName)
  PYB11_PROPERTY(DataArrayPath EdgeDataContainerName READ getEdgeDataContainerName WRITE setEdgeDataContainerName)
  PYB11_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)
  PYB11_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)
  PYB11_PROPERTY(QString InputFile READ getInputFile WRITE setInputFile)
  PYB11_PROPERTY(float BurgersVector READ getBurgersVector WRITE setBurgersVector)
  PYB11_PROPERTY(QString NumberOfArmsArrayName READ getNumberOfArmsArrayName WRITE setNumberOfArmsArrayName)
  PYB11_PROPERTY(QString NodeConstraintsArrayName READ getNodeConstraintsArrayName WRITE setNodeConstraintsArrayName)
  PYB11_PROPERTY(QString BurgersVectorsArrayName READ getBurgersVectorsArrayName WRITE setBurgersVectorsArrayName)
  PYB11_PROPERTY(QString SlipPlaneNormalsArrayName READ getSlipPlaneNormalsArrayName WRITE setSlipPlaneNormalsArrayName)
  PYB11_PROPERTY(QString DomainBoundsArrayName READ getDomainBoundsArrayName WRITE setDomainBoundsArrayName)
#endif

public:
  using Self = ParaDisReader;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ParaDisReader> New();

  /**
   * @brief Returns the name of the class for ParaDisReader
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ParaDisReader
   */
  static QString ClassName();

  virtual ~ParaDisReader();
  /**
   * @brief Setter property for EdgeDataContainerName
   */
  void setEdgeDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for EdgeDataContainerName
   * @return Value of EdgeDataContainerName
   */
  DataArrayPath getEdgeDataContainerName() const;

  Q_PROPERTY(DataArrayPath EdgeDataContainerName READ getEdgeDataContainerName WRITE setEdgeDataContainerName)
  /**
   * @brief Setter property for VertexAttributeMatrixName
   */
  void setVertexAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for VertexAttributeMatrixName
   * @return Value of VertexAttributeMatrixName
   */
  QString getVertexAttributeMatrixName() const;

  Q_PROPERTY(QString VertexAttributeMatrixName READ getVertexAttributeMatrixName WRITE setVertexAttributeMatrixName)
  /**
   * @brief Setter property for EdgeAttributeMatrixName
   */
  void setEdgeAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for EdgeAttributeMatrixName
   * @return Value of EdgeAttributeMatrixName
   */
  QString getEdgeAttributeMatrixName() const;

  Q_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)

  /**
   * @brief Setter property for InputFile
   */
  void setInputFile(const QString& value);
  /**
   * @brief Getter property for InputFile
   * @return Value of InputFile
   */
  QString getInputFile() const;

  Q_PROPERTY(QString InputFile READ getInputFile WRITE setInputFile)
  /**
   * @brief Setter property for BurgersVector
   */
  void setBurgersVector(float value);
  /**
   * @brief Getter property for BurgersVector
   * @return Value of BurgersVector
   */
  float getBurgersVector() const;

  Q_PROPERTY(float BurgersVector READ getBurgersVector WRITE setBurgersVector)

  /**
   * @brief Setter property for NumberOfArmsArrayName
   */
  void setNumberOfArmsArrayName(const QString& value);
  /**
   * @brief Getter property for NumberOfArmsArrayName
   * @return Value of NumberOfArmsArrayName
   */
  QString getNumberOfArmsArrayName() const;

  Q_PROPERTY(QString NumberOfArmsArrayName READ getNumberOfArmsArrayName WRITE setNumberOfArmsArrayName)

  /**
   * @brief Setter property for NodeConstraintsArrayName
   */
  void setNodeConstraintsArrayName(const QString& value);
  /**
   * @brief Getter property for NodeConstraintsArrayName
   * @return Value of NodeConstraintsArrayName
   */
  QString getNodeConstraintsArrayName() const;

  Q_PROPERTY(QString NodeConstraintsArrayName READ getNodeConstraintsArrayName WRITE setNodeConstraintsArrayName)

  /**
   * @brief Setter property for BurgersVectorsArrayName
   */
  void setBurgersVectorsArrayName(const QString& value);
  /**
   * @brief Getter property for BurgersVectorsArrayName
   * @return Value of BurgersVectorsArrayName
   */
  QString getBurgersVectorsArrayName() const;

  Q_PROPERTY(QString BurgersVectorsArrayName READ getBurgersVectorsArrayName WRITE setBurgersVectorsArrayName)

  /**
   * @brief Setter property for SlipPlaneNormalsArrayName
   */
  void setSlipPlaneNormalsArrayName(const QString& value);
  /**
   * @brief Getter property for SlipPlaneNormalsArrayName
   * @return Value of SlipPlaneNormalsArrayName
   */
  QString getSlipPlaneNormalsArrayName() const;

  Q_PROPERTY(QString SlipPlaneNormalsArrayName READ getSlipPlaneNormalsArrayName WRITE setSlipPlaneNormalsArrayName)

  /**
   * @brief Setter property for DomainBoundsArrayName
   */
  void setDomainBoundsArrayName(const QString& value);
  /**
   * @brief Getter property for DomainBoundsArrayName
   * @return Value of DomainBoundsArrayName
   */
  QString getDomainBoundsArrayName() const;

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief This method will read the options from a file
   * @param reader The reader that is used to read the options from a file
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  ParaDisReader();

  /**
   * @brief readHeader
   * @return
   */
  virtual int readHeader();

  /**
   * @brief readFile
   * @return
   */
  virtual int readFile();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief updateVertexInstancePointers
   */
  void updateVertexInstancePointers();

  /**
   * @brief updateEdgeInstancePointers
   */
  void updateEdgeInstancePointers();

private:
  std::weak_ptr<DataArray<int32_t>> m_NumberOfArmsPtr;
  int32_t* m_NumberOfArms = nullptr;
  std::weak_ptr<DataArray<int32_t>> m_NodeConstraintsPtr;
  int32_t* m_NodeConstraints = nullptr;
  std::weak_ptr<DataArray<float>> m_BurgersVectorsPtr;
  float* m_BurgersVectors = nullptr;
  std::weak_ptr<DataArray<float>> m_SlipPlaneNormalsPtr;
  float* m_SlipPlaneNormals = nullptr;
  std::weak_ptr<DataArray<float>> m_DomainBoundsPtr;
  float* m_DomainBounds = nullptr;

  DataArrayPath m_EdgeDataContainerName = {};
  QString m_VertexAttributeMatrixName = {};
  QString m_EdgeAttributeMatrixName = {};
  QString m_InputFile = {};
  float m_BurgersVector = {};
  QString m_NumberOfArmsArrayName = {};
  QString m_NodeConstraintsArrayName = {};
  QString m_BurgersVectorsArrayName = {};
  QString m_SlipPlaneNormalsArrayName = {};
  QString m_DomainBoundsArrayName = {};

  QFile m_InStream;

  int m_NumVerts;
  int m_NumEdges;

  int m_FileVersion;

public:
  ParaDisReader(const ParaDisReader&) = delete;            // Copy Constructor Not Implemented
  ParaDisReader(ParaDisReader&&) = delete;                 // Move Constructor Not Implemented
  ParaDisReader& operator=(const ParaDisReader&) = delete; // Copy Assignment Not Implemented
  ParaDisReader& operator=(ParaDisReader&&) = delete;      // Move Assignment Not Implemented
};
