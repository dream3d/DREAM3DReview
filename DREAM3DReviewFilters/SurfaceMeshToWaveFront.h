/*
 * Your License or Copyright can go here
 */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/Filtering/AbstractFilter.h"


#include "DREAM3DReview/DREAM3DReviewPlugin.h"
    
/**
 * @brief The SurfaceMeshToWaveFront class. See [Filter documentation](@ref surfacemeshtowavefront) for details.
 */
class DREAM3DReview_EXPORT SurfaceMeshToWaveFront : public AbstractFilter
{
  Q_OBJECT

#ifdef SIMPL_ENABLE_PYTHON
  // clang-format off
  PYB11_CREATE_BINDINGS(SurfaceMeshToWaveFront SUPERCLASS AbstractFilter)
  PYB11_FILTER_PARAMETER(QString, OutputWaveFrontFile)
//  PYB11_FILTER_PARAMETER(DataArrayPath, SurfaceMeshFaceLabelsArrayPath)
  PYB11_FILTER_PARAMETER(DataArrayPath, SurfaceMeshNodeTypeArrayPath)
//  PYB11_FILTER_PARAMETER(QVector<DataArrayPath>, SelectedFaceArrays)
//  PYB11_FILTER_PARAMETER(QVector<DataArrayPath>, SelectedVertexArrays)

  // clang-format on
#endif

public:

  using Self = SurfaceMeshToWaveFront;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<SurfaceMeshToWaveFront> New();

  /**
   * @brief Setter property for OutputWaveFrontFile
   */
  void setOutputWaveFrontFile(const QString& value);
  /**
   * @brief Getter property for OutputWaveFrontFile
   * @return Value of OutputWaveFrontFile
   */
  QString getOutputWaveFrontFile() const;

  Q_PROPERTY(QString OutputWaveFrontFile READ getOutputWaveFrontFile WRITE setOutputWaveFrontFile)

  //  /**
  //   * @brief Setter property for DataContainerName
  //   */
  //  void setDataContainerPath(const DataArrayPath& value);
  //  /**
  //   * @brief Getter property for DataContainerPath
  //   * @return Value of DataContainerPath
  //   */
  //  DataArrayPath getDataContainerPath() const;

  //  Q_PROPERTY(DataArrayPath DataContainerPath READ getDataContainerPath WRITE setDataContainerPath)

  //  /**
  //   * @brief Setter property for SurfaceMeshFaceLabelsArrayPath
  //   */
  //  void setSurfaceMeshFaceLabelsArrayPath(const DataArrayPath& value);
  //  /**
  //   * @brief Getter property for SurfaceMeshFaceLabelsArrayPath
  //   * @return Value of SurfaceMeshFaceLabelsArrayPath
  //   */
  //  DataArrayPath getSurfaceMeshFaceLabelsArrayPath() const;

  //  Q_PROPERTY(DataArrayPath SurfaceMeshFaceLabelsArrayPath READ getSurfaceMeshFaceLabelsArrayPath WRITE setSurfaceMeshFaceLabelsArrayPath)

  /**
   * @brief Setter property for SurfaceMeshNodeTypeArrayPath
   */
  void setSurfaceMeshNodeTypeArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SurfaceMeshNodeTypeArrayPath
   * @return Value of SurfaceMeshNodeTypeArrayPath
   */
  DataArrayPath getSurfaceMeshNodeTypeArrayPath() const;

  Q_PROPERTY(DataArrayPath SurfaceMeshNodeTypeArrayPath READ getSurfaceMeshNodeTypeArrayPath WRITE setSurfaceMeshNodeTypeArrayPath)

  //  /**
  //   * @brief Setter property for SelectedFaceArrays
  //   */
  //  void setSelectedFaceArrays(const QVector<DataArrayPath>& value);
  //  /**
  //   * @brief Getter property for SelectedFaceArrays
  //   * @return Value of SelectedFaceArrays
  //   */
  //  QVector<DataArrayPath> getSelectedFaceArrays() const;

  //  Q_PROPERTY(QVector<DataArrayPath> SelectedFaceArrays READ getSelectedFaceArrays WRITE setSelectedFaceArrays)

  //  /**
  //   * @brief Setter property for SelectedVertexArrays
  //   */
  //  void setSelectedVertexArrays(const QVector<DataArrayPath>& value);
  //  /**
  //   * @brief Getter property for SelectedVertexArrays
  //   * @return Value of SelectedVertexArrays
  //   */
  //  QVector<DataArrayPath> getSelectedVertexArrays() const;

  //  Q_PROPERTY(QVector<DataArrayPath> SelectedVertexArrays READ getSelectedVertexArrays WRITE setSelectedVertexArrays)

  /**
   * @brief Returns the name of the class for SurfaceMeshToWaveFront
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for SurfaceMeshToWaveFront
   */
  static QString ClassName();

  ~SurfaceMeshToWaveFront() override;
  


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
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief preflight Reimplemented from @see AbstractFilter class
   */
  void preflight() override;

signals:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter * filter);

  /**
   * @brief parametersChanged Emitted when any Filter parameter is changed internally
   */
  void parametersChanged();

  /**
   * @brief preflightAboutToExecute Emitted just before calling dataCheck()
   */
  void preflightAboutToExecute();

  /**
   * @brief preflightExecuted Emitted just after calling dataCheck()
   */
  void preflightExecuted();

protected:
  SurfaceMeshToWaveFront();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  //  std::weak_ptr<DataArray<int32_t>> m_SurfaceMeshFaceLabelsPtr;
  //  int32_t* m_SurfaceMeshFaceLabels = nullptr;
  std::weak_ptr<DataArray<int8_t>> m_SurfaceMeshNodeTypePtr;
  int8_t* m_SurfaceMeshNodeType = nullptr;
  //  DataArrayPath m_DataContainerPath;

  QString m_OutputWaveFrontFile = {};
  //  DataArrayPath m_SurfaceMeshFaceLabelsArrayPath = {};
  DataArrayPath m_SurfaceMeshNodeTypeArrayPath = {};
  //  QVector<DataArrayPath> m_SelectedFaceArrays = {};
  //  QVector<DataArrayPath> m_SelectedVertexArrays = {};

public:
  SurfaceMeshToWaveFront(const SurfaceMeshToWaveFront&) = delete;            // Copy Constructor Not Implemented
  SurfaceMeshToWaveFront& operator=(const SurfaceMeshToWaveFront&) = delete; // Copy Assignment Not Implemented
  SurfaceMeshToWaveFront(SurfaceMeshToWaveFront &&) = delete;                // Move Constructor Not Implemented
  SurfaceMeshToWaveFront& operator=(SurfaceMeshToWaveFront&&) = delete;      // Move Assignment Not Implemented
};

