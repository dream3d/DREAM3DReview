/* ============================================================================
 * Software developed by US federal government employees (including military personnel)
 * as part of their official duties is not subject to copyright protection and is
 * considered "public domain" (see 17 USC Section 105). Public domain software can be used
 * by anyone for any purpose, and cannot be released under a copyright license
 * (including typical open source software licenses).
 *
 * This source code file was originally written by United States DoD employees. The
 * original source code files are released into the Public Domain.
 *
 * Subsequent changes to the codes by others may elect to add a copyright and license
 * for those changes.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include <memory>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

class IDataArray;
using IDataArrayWkPtrType = std::weak_ptr<IDataArray>;

/**
 * @brief The FindNeighborListStatistics class. See [Filter documentation](@ref FindNeighborListStatistics) for details.
 */
class DREAM3DReview_EXPORT FindNeighborListStatistics : public AbstractFilter
{
  Q_OBJECT

  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(FindNeighborListStatistics SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(FindNeighborListStatistics)
  PYB11_FILTER_NEW_MACRO(FindNeighborListStatistics)

  PYB11_PROPERTY(bool FindLength READ getFindLength WRITE setFindLength)
  PYB11_PROPERTY(bool FindMin READ getFindMin WRITE setFindMin)
  PYB11_PROPERTY(bool FindMax READ getFindMax WRITE setFindMax)
  PYB11_PROPERTY(bool FindMean READ getFindMean WRITE setFindMean)
  PYB11_PROPERTY(bool FindMedian READ getFindMedian WRITE setFindMedian)
  PYB11_PROPERTY(bool FindStdDeviation READ getFindStdDeviation WRITE setFindStdDeviation)
  PYB11_PROPERTY(bool FindSummation READ getFindSummation WRITE setFindSummation)

  PYB11_PROPERTY(DataArrayPath DestinationAttributeMatrix READ getDestinationAttributeMatrix WRITE setDestinationAttributeMatrix)

  PYB11_PROPERTY(QString LengthArrayName READ getLengthArrayName WRITE setLengthArrayName)
  PYB11_PROPERTY(QString MinimumArrayName READ getMinimumArrayName WRITE setMinimumArrayName)
  PYB11_PROPERTY(QString MaximumArrayName READ getMaximumArrayName WRITE setMaximumArrayName)
  PYB11_PROPERTY(QString MeanArrayName READ getMeanArrayName WRITE setMeanArrayName)
  PYB11_PROPERTY(QString MedianArrayName READ getMedianArrayName WRITE setMedianArrayName)
  PYB11_PROPERTY(QString StdDeviationArrayName READ getStdDeviationArrayName WRITE setStdDeviationArrayName)
  PYB11_PROPERTY(QString SummationArrayName READ getSummationArrayName WRITE setSummationArrayName)

  PYB11_PROPERTY(DataArrayPath SelectedArrayPath READ getSelectedArrayPath WRITE setSelectedArrayPath)

  PYB11_END_BINDINGS()
  // End Python bindings declarations

public:
  using Self = FindNeighborListStatistics;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindNeighborListStatistics> New();

  /**
   * @brief Returns the name of the class for FindNeighborListStatistics
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindNeighborListStatistics
   */
  static QString ClassName();

  ~FindNeighborListStatistics() override;

  /**
   * @brief Setter property for FindLength
   */
  void setFindLength(bool value);
  /**
   * @brief Getter property for FindLength
   * @return Value of FindLength
   */
  bool getFindLength() const;
  Q_PROPERTY(bool FindLength READ getFindLength WRITE setFindLength)

  /**
   * @brief Setter property for FindMin
   */
  void setFindMin(bool value);
  /**
   * @brief Getter property for FindMin
   * @return Value of FindMin
   */
  bool getFindMin() const;
  Q_PROPERTY(bool FindMin READ getFindMin WRITE setFindMin)

  /**
   * @brief Setter property for FindMax
   */
  void setFindMax(bool value);
  /**
   * @brief Getter property for FindMax
   * @return Value of FindMax
   */
  bool getFindMax() const;
  Q_PROPERTY(bool FindMax READ getFindMax WRITE setFindMax)

  /**
   * @brief Setter property for FindMean
   */
  void setFindMean(bool value);
  /**
   * @brief Getter property for FindMean
   * @return Value of FindMean
   */
  bool getFindMean() const;
  Q_PROPERTY(bool FindMean READ getFindMean WRITE setFindMean)

  /**
   * @brief Setter property for FindMedian
   */
  void setFindMedian(bool value);
  /**
   * @brief Getter property for FindMedian
   * @return Value of FindMedian
   */
  bool getFindMedian() const;
  Q_PROPERTY(bool FindMedian READ getFindMedian WRITE setFindMedian)

  /**
   * @brief Setter property for FindStdDeviation
   */
  void setFindStdDeviation(bool value);
  /**
   * @brief Getter property for FindStdDeviation
   * @return Value of FindStdDeviation
   */
  bool getFindStdDeviation() const;
  Q_PROPERTY(bool FindStdDeviation READ getFindStdDeviation WRITE setFindStdDeviation)

  /**
   * @brief Setter property for FindSummation
   */
  void setFindSummation(bool value);
  /**
   * @brief Getter property for FindSummation
   * @return Value of FindSummation
   */
  bool getFindSummation() const;
  Q_PROPERTY(bool FindSummation READ getFindSummation WRITE setFindSummation)

  /**
   * @brief Setter property for DestinationAttributeMatrix
   */
  void setDestinationAttributeMatrix(const DataArrayPath& value);
  /**
   * @brief Getter property for DestinationAttributeMatrix
   * @return Value of DestinationAttributeMatrix
   */
  DataArrayPath getDestinationAttributeMatrix() const;
  Q_PROPERTY(DataArrayPath DestinationAttributeMatrix READ getDestinationAttributeMatrix WRITE setDestinationAttributeMatrix)

  /**
   * @brief Setter property for LengthArrayName
   */
  void setLengthArrayName(const QString& value);
  /**
   * @brief Getter property for LengthArrayName
   * @return Value of LengthArrayName
   */
  QString getLengthArrayName() const;
  Q_PROPERTY(QString LengthArrayName READ getLengthArrayName WRITE setLengthArrayName)

  /**
   * @brief Setter property for MinimumArrayName
   */
  void setMinimumArrayName(const QString& value);
  /**
   * @brief Getter property for MinimumArrayName
   * @return Value of MinimumArrayName
   */
  QString getMinimumArrayName() const;
  Q_PROPERTY(QString MinimumArrayName READ getMinimumArrayName WRITE setMinimumArrayName)

  /**
   * @brief Setter property for MaximumArrayName
   */
  void setMaximumArrayName(const QString& value);
  /**
   * @brief Getter property for MaximumArrayName
   * @return Value of MaximumArrayName
   */
  QString getMaximumArrayName() const;
  Q_PROPERTY(QString MaximumArrayName READ getMaximumArrayName WRITE setMaximumArrayName)

  /**
   * @brief Setter property for MeanArrayName
   */
  void setMeanArrayName(const QString& value);
  /**
   * @brief Getter property for MeanArrayName
   * @return Value of MeanArrayName
   */
  QString getMeanArrayName() const;
  Q_PROPERTY(QString MeanArrayName READ getMeanArrayName WRITE setMeanArrayName)

  /**
   * @brief Setter property for MedianArrayName
   */
  void setMedianArrayName(const QString& value);
  /**
   * @brief Getter property for MedianArrayName
   * @return Value of MedianArrayName
   */
  QString getMedianArrayName() const;
  Q_PROPERTY(QString MedianArrayName READ getMedianArrayName WRITE setMedianArrayName)

  /**
   * @brief Setter property for StdDeviationArrayName
   */
  void setStdDeviationArrayName(const QString& value);
  /**
   * @brief Getter property for StdDeviationArrayName
   * @return Value of StdDeviationArrayName
   */
  QString getStdDeviationArrayName() const;
  Q_PROPERTY(QString StdDeviationArrayName READ getStdDeviationArrayName WRITE setStdDeviationArrayName)

  /**
   * @brief Setter property for SummationArrayName
   */
  void setSummationArrayName(const QString& value);
  /**
   * @brief Getter property for SummationArrayName
   * @return Value of SummationArrayName
   */
  QString getSummationArrayName() const;
  Q_PROPERTY(QString SummationArrayName READ getSummationArrayName WRITE setSummationArrayName)

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

protected:
  /**
   * @brief createCompatibleArrays Creates the output statistics arrays with compatible types to the
   * incoming NeighborList
   * @param dataArrayPaths Vector used to store paths for later tuple checking
   */
  template <typename T>
  void createCompatibleArrays(QVector<DataArrayPath>& dataArrayPaths);

  FindNeighborListStatistics();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  IDataArrayWkPtrType m_MinimumPtr;
  IDataArrayWkPtrType m_MaximumPtr;
  IDataArrayWkPtrType m_InputArrayPtr;

  std::weak_ptr<SizeTArrayType> m_LengthPtr;
  std::weak_ptr<FloatArrayType> m_MeanPtr;
  std::weak_ptr<FloatArrayType> m_MedianPtr;
  std::weak_ptr<FloatArrayType> m_StandardDeviationPtr;
  std::weak_ptr<FloatArrayType> m_SummationPtr;

  DataArrayPath m_DestinationAttributeMatrix = {"", "", ""};

  bool m_FindLength = false;
  bool m_FindMin = false;
  bool m_FindMax = false;
  bool m_FindMean = false;
  bool m_FindMedian = false;
  bool m_FindStdDeviation = false;
  bool m_FindSummation = false;

  QString m_LengthArrayName = {"Length"};
  QString m_MinimumArrayName = {"Minimum"};
  QString m_MaximumArrayName = {"Maximum"};
  QString m_MeanArrayName = {"Mean"};
  QString m_MedianArrayName = {"Median"};
  QString m_StdDeviationArrayName = {"StandardDeviation"};
  QString m_SummationArrayName = {"Summation"};

  DataArrayPath m_SelectedArrayPath = {};

public:
  FindNeighborListStatistics(const FindNeighborListStatistics&) = delete;            // Copy Constructor Not Implemented
  FindNeighborListStatistics(FindNeighborListStatistics&&) = delete;                 // Move Constructor Not Implemented
  FindNeighborListStatistics& operator=(const FindNeighborListStatistics&) = delete; // Copy Assignment Not Implemented
  FindNeighborListStatistics& operator=(FindNeighborListStatistics&&) = delete;      // Move Assignment Not Implemented
};
