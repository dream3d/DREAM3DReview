/* ============================================================================
* Software developed by US federal government employees (including military personnel) 
* as part of their official duties is not subject to copyright protection and is 
* considered “public domain” (see 17 USC Section 105). Public domain software can be used 
* by anyone for any purpose, and cannot be released under a copyright license 
* (including typical open source software licenses).
* 
* This source code file was originally written by United States DoD employees. The
* original source code files are released into the Public Domain.
* 
* Subsequent changes to the codes by others may elect to add a copyright and license
* for those changes.
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _findlayerstatistics_h_
#define _findlayerstatistics_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

class IDataArray;
using IDataArrayWkPtrType = std::weak_ptr<IDataArray>;

/**
 * @brief The FindLayerStatistics class. See [Filter documentation](@ref findlayerstatistics) for details.
 */
class FindLayerStatistics : public AbstractFilter
{
  Q_OBJECT
public:
  using Self = FindLayerStatistics;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindLayerStatistics> New();

  /**
   * @brief Returns the name of the class for FindLayerStatistics
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindLayerStatistics
   */
  static QString ClassName();

  ~FindLayerStatistics() override;

  /**
   * @brief Setter property for SelectedArrayPath
   */
  void setSelectedArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SelectedArrayPath
   * @return Value of SelectedArrayPath
   */
  DataArrayPath getSelectedArrayPath() const;

  Q_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)

  /**
   * @brief Setter property for Plane
   */
  void setPlane(unsigned int value);
  /**
   * @brief Getter property for Plane
   * @return Value of Plane
   */
  unsigned int getPlane() const;

  Q_PROPERTY(unsigned int Plane READ getPlane WRITE setPlane)

  /**
   * @brief Setter property for LayerIDsArrayName
   */
  void setLayerIDsArrayName(const QString& value);
  /**
   * @brief Getter property for LayerIDsArrayName
   * @return Value of LayerIDsArrayName
   */
  QString getLayerIDsArrayName() const;

  Q_PROPERTY(QString LayerIDsArrayName READ getLayerIDsArrayName WRITE setLayerIDsArrayName)

  /**
   * @brief Setter property for LayerAttributeMatrixName
   */
  void setLayerAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for LayerAttributeMatrixName
   * @return Value of LayerAttributeMatrixName
   */
  QString getLayerAttributeMatrixName() const;

  Q_PROPERTY(QString LayerAttributeMatrixName READ getLayerAttributeMatrixName WRITE setLayerAttributeMatrixName)

  /**
   * @brief Setter property for LayerMinArrayName
   */
  void setLayerMinArrayName(const QString& value);
  /**
   * @brief Getter property for LayerMinArrayName
   * @return Value of LayerMinArrayName
   */
  QString getLayerMinArrayName() const;

  Q_PROPERTY(QString LayerMinArrayName READ getLayerMinArrayName WRITE setLayerMinArrayName)

  /**
   * @brief Setter property for LayerMaxArrayName
   */
  void setLayerMaxArrayName(const QString& value);
  /**
   * @brief Getter property for LayerMaxArrayName
   * @return Value of LayerMaxArrayName
   */
  QString getLayerMaxArrayName() const;

  Q_PROPERTY(QString LayerMaxArrayName READ getLayerMaxArrayName WRITE setLayerMaxArrayName)

  /**
   * @brief Setter property for LayerAvgArrayName
   */
  void setLayerAvgArrayName(const QString& value);
  /**
   * @brief Getter property for LayerAvgArrayName
   * @return Value of LayerAvgArrayName
   */
  QString getLayerAvgArrayName() const;

  Q_PROPERTY(QString LayerAvgArrayName READ getLayerAvgArrayName WRITE setLayerAvgArrayName)

  /**
   * @brief Setter property for LayerStdArrayName
   */
  void setLayerStdArrayName(const QString& value);
  /**
   * @brief Getter property for LayerStdArrayName
   * @return Value of LayerStdArrayName
   */
  QString getLayerStdArrayName() const;

  Q_PROPERTY(QString LayerStdArrayName READ getLayerStdArrayName WRITE setLayerStdArrayName)

  /**
   * @brief Setter property for LayerVarArrayName
   */
  void setLayerVarArrayName(const QString& value);
  /**
   * @brief Getter property for LayerVarArrayName
   * @return Value of LayerVarArrayName
   */
  QString getLayerVarArrayName() const;

  Q_PROPERTY(QString LayerVarArrayName READ getLayerVarArrayName WRITE setLayerVarArrayName)

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
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
  * @brief getUuid Return the unique identifier for this filter.
  * @return A QUuid object.
  */
  QUuid getUuid() const override;

protected:
  FindLayerStatistics();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  IDataArrayWkPtrType m_InDataPtr;
  void* m_InData = nullptr;

  std::weak_ptr<DataArray<int32_t>> m_LayerIDsPtr;
  int32_t* m_LayerIDs = nullptr;
  std::weak_ptr<DataArray<float>> m_LayerMinPtr;
  float* m_LayerMin = nullptr;
  std::weak_ptr<DataArray<float>> m_LayerMaxPtr;
  float* m_LayerMax = nullptr;
  std::weak_ptr<DataArray<float>> m_LayerAvgPtr;
  float* m_LayerAvg = nullptr;
  std::weak_ptr<DataArray<float>> m_LayerStdPtr;
  float* m_LayerStd = nullptr;
  std::weak_ptr<DataArray<float>> m_LayerVarPtr;
  float* m_LayerVar = nullptr;

  DataArrayPath m_SelectedArrayPath = {};
  unsigned int m_Plane = {};
  QString m_LayerIDsArrayName = {};
  QString m_LayerAttributeMatrixName = {};
  QString m_LayerMinArrayName = {};
  QString m_LayerMaxArrayName = {};
  QString m_LayerAvgArrayName = {};
  QString m_LayerStdArrayName = {};
  QString m_LayerVarArrayName = {};

  FindLayerStatistics(const FindLayerStatistics&) = delete; // Copy Constructor Not Implemented
  FindLayerStatistics(FindLayerStatistics&&) = delete;      // Move Constructor Not Implemented
  void operator=(const FindLayerStatistics&) = delete;      // Operator '=' Not Implemented
};

#endif /* _findlayerstatistics_h_ */
