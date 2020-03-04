/* ============================================================================
 * Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2012 Dr. Michael A. Groeber (US Air Force Research Laboratories)
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

#include <QtCore/QDateTime>
#include <QtCore/QScopedPointer>
#include <QtCore/QString>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataArrays/StringDataArray.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewFilters/HEDM/MicPhase.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"

class MicReader;

struct Mic_Private_Data
{
  SizeVec3Type dims;
  FloatVec3Type resolution;
  FloatVec3Type origin;
  QVector<MicPhase::Pointer> phases;
};

enum MIC_READ_FLAG
{
  MIC_FULL_FILE,
  MIC_HEADER_ONLY
};

// our PIMPL private class
class ReadMicDataPrivate;

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @class ReadMicData ReadMicData.h /FilterCategoryFilters/ReadMicData.h
 * @brief
 * @author
 * @date
 * @version 1.0
 */
class DREAM3DReview_EXPORT ReadMicData : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  PYB11_CREATE_BINDINGS(ReadMicData SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(ReadMicData)
  PYB11_FILTER_NEW_MACRO(ReadMicData)
  PYB11_FILTER_PARAMETER(DataArrayPath, DataContainerName)
  PYB11_FILTER_PARAMETER(QString, CellEnsembleAttributeMatrixName)
  PYB11_FILTER_PARAMETER(QString, CellAttributeMatrixName)
  PYB11_FILTER_PARAMETER(bool, FileWasRead)
  PYB11_FILTER_PARAMETER(QString, InputFile)
  PYB11_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)
  PYB11_PROPERTY(QString CellEnsembleAttributeMatrixName READ getCellEnsembleAttributeMatrixName WRITE setCellEnsembleAttributeMatrixName)
  PYB11_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)
  PYB11_PROPERTY(bool FileWasRead READ getFileWasRead WRITE setFileWasRead)
  PYB11_PROPERTY(QString InputFile READ getInputFile WRITE setInputFile)
#endif

  Q_DECLARE_PRIVATE(ReadMicData)

public:
  using Self = ReadMicData;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ReadMicData> New();

  /**
   * @brief Returns the name of the class for ReadMicData
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ReadMicData
   */
  static QString ClassName();

  ~ReadMicData() override;

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  DataArrayPath getDataContainerName() const;

  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /**
   * @brief Setter property for CellEnsembleAttributeMatrixName
   */
  void setCellEnsembleAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for CellEnsembleAttributeMatrixName
   * @return Value of CellEnsembleAttributeMatrixName
   */
  QString getCellEnsembleAttributeMatrixName() const;

  Q_PROPERTY(QString CellEnsembleAttributeMatrixName READ getCellEnsembleAttributeMatrixName WRITE setCellEnsembleAttributeMatrixName)

  /**
   * @brief Setter property for CellAttributeMatrixName
   */
  void setCellAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for CellAttributeMatrixName
   * @return Value of CellAttributeMatrixName
   */
  QString getCellAttributeMatrixName() const;

  Q_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)

  /**
   * @brief Setter property for FileWasRead
   */
  void setFileWasRead(bool value);
  /**
   * @brief Getter property for FileWasRead
   * @return Value of FileWasRead
   */
  bool getFileWasRead() const;

  Q_PROPERTY(bool FileWasRead READ getFileWasRead)

  /**
   * @brief Setter property for PhaseNameArrayName
   */
  void setPhaseNameArrayName(const QString& value);
  /**
   * @brief Getter property for PhaseNameArrayName
   * @return Value of PhaseNameArrayName
   */
  QString getPhaseNameArrayName() const;

  /**
   * @brief Setter property for MaterialNameArrayName
   */
  void setMaterialNameArrayName(const QString& value);
  /**
   * @brief Getter property for MaterialNameArrayName
   * @return Value of MaterialNameArrayName
   */
  QString getMaterialNameArrayName() const;

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
   * @brief Setter property for CellEulerAnglesArrayName
   */
  void setCellEulerAnglesArrayName(const QString& value);
  /**
   * @brief Getter property for CellEulerAnglesArrayName
   * @return Value of CellEulerAnglesArrayName
   */
  QString getCellEulerAnglesArrayName() const;

  /**
   * @brief Setter property for CellPhasesArrayName
   */
  void setCellPhasesArrayName(const QString& value);
  /**
   * @brief Getter property for CellPhasesArrayName
   * @return Value of CellPhasesArrayName
   */
  QString getCellPhasesArrayName() const;

  /**
   * @brief Setter property for CrystalStructuresArrayName
   */
  void setCrystalStructuresArrayName(const QString& value);
  /**
   * @brief Getter property for CrystalStructuresArrayName
   * @return Value of CrystalStructuresArrayName
   */
  QString getCrystalStructuresArrayName() const;

  /**
   * @brief Setter property for LatticeConstantsArrayName
   */
  void setLatticeConstantsArrayName(const QString& value);
  /**
   * @brief Getter property for LatticeConstantsArrayName
   * @return Value of LatticeConstantsArrayName
   */
  QString getLatticeConstantsArrayName() const;

  /**
   * @brief This returns the group that the filter belonds to. You can select
   * a different group if you want. The string returned here will be displayed
   * in the GUI for the filter
   */
  QString getCompiledLibraryName() const override;
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;
  QString getGroupName() const override;

  /**
   * @brief This returns the group that the filter belonds to. You can select
   * a different group if you want. The string returned here will be displayed
   * in the GUI for the filter
   */
  QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

  /**
   * @brief This returns a string that is displayed in the GUI. It should be readable
   * and understandable by humans.
   */
  QString getHumanLabel() const override;

  /**
   * @brief This method will instantiate all the end user settable options/parameters
   * for this filter
   */
  void setupFilterParameters() override;

  /**
   * @brief This method will read the options from a file
   * @param reader The reader that is used to read the options from a file
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief This function runs some sanity checks on the DataContainer and inputs
   * in an attempt to ensure the filter can process the inputs.
   */
  void preflight() override;

  /* These are non-exposed to the user through the GUI. Manual Pipelines are OK to set them */
  /**
   * @brief Setter property for RefFrameZDir
   */
  void setRefFrameZDir(uint32_t value);
  /**
   * @brief Getter property for RefFrameZDir
   * @return Value of RefFrameZDir
   */
  uint32_t getRefFrameZDir() const;

  /**
   * @brief Setter property for Manufacturer
   */
  void setManufacturer(const Ebsd::OEM& value);
  /**
   * @brief Getter property for Manufacturer
   * @return Value of Manufacturer
   */
  Ebsd::OEM getManufacturer() const;

  /**
   * @brief Setter property for InputFile_Cache
   */
  void setInputFile_Cache(const QString& value);
  /**
   * @brief Getter property for InputFile_Cache
   * @return Value of InputFile_Cache
   */
  QString getInputFile_Cache() const;

  /**
   * @brief Setter property for TimeStamp_Cache
   */
  void setTimeStamp_Cache(const QDateTime& value);
  /**
   * @brief Getter property for TimeStamp_Cache
   * @return Value of TimeStamp_Cache
   */
  QDateTime getTimeStamp_Cache() const;

  /**
   * @brief Setter property for Data
   */
  void setData(const Mic_Private_Data& value);
  /**
   * @brief Getter property for Data
   * @return Value of Data
   */
  Mic_Private_Data getData() const;

  Q_PROPERTY(Mic_Private_Data Data READ getData WRITE setData)

public slots:
  void flushCache();

protected:
  ReadMicData();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

  /**
   * @brief readMicFile This reads the HEDM .mic file and puts the data into the Voxel Data container
   */
  void readMicFile();

  /**
   * @brief This method reads the values for the phase type, crystal structure
   * and precipitate fractions from the EBSD file.
   * @param reader The EbsdReader instance
   * @return Zero/Positive on Success - Negative on error.
   */
  int loadMaterialInfo(MicReader* reader);

  /**
   * @brief populateMicData
   * @param reader
   * @param m
   * @param dims
   */
  void populateMicData(MicReader* reader, DataContainer::Pointer m, std::vector<size_t> dims, MIC_READ_FLAG = MIC_FULL_FILE);

private:
  std::weak_ptr<DataArray<int32_t>> m_CellPhasesPtr;
  int32_t* m_CellPhases = nullptr;
  std::weak_ptr<DataArray<float>> m_CellEulerAnglesPtr;
  float* m_CellEulerAngles = nullptr;
  std::weak_ptr<DataArray<uint32_t>> m_CrystalStructuresPtr;
  uint32_t* m_CrystalStructures = nullptr;
  std::weak_ptr<DataArray<float>> m_LatticeConstantsPtr;
  float* m_LatticeConstants = nullptr;

  DataArrayPath m_DataContainerName = {};
  QString m_CellEnsembleAttributeMatrixName = {};
  QString m_CellAttributeMatrixName = {};
  bool m_FileWasRead = {};
  QString m_PhaseNameArrayName = {};
  QString m_MaterialNameArrayName = {};
  QString m_InputFile = {};
  QString m_CellEulerAnglesArrayName = {};
  QString m_CellPhasesArrayName = {};
  QString m_CrystalStructuresArrayName = {};
  QString m_LatticeConstantsArrayName = {};
  uint32_t m_RefFrameZDir = {};
  Ebsd::OEM m_Manufacturer = {};

  QScopedPointer<ReadMicDataPrivate> const d_ptr;

public:
  ReadMicData(const ReadMicData&) = delete;            // Copy Constructor Not Implemented
  ReadMicData(ReadMicData&&) = delete;                 // Move Constructor Not Implemented
  ReadMicData& operator=(const ReadMicData&) = delete; // Copy Assignment Not Implemented
  ReadMicData& operator=(ReadMicData&&) = delete;      // Move Assignment Not Implemented
};

Q_DECLARE_METATYPE(Mic_Private_Data)
