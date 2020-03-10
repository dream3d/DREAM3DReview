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
#ifndef _finddynamicarraystatistics_h_
#define _finddynamicarraystatistics_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

class IDataArray;
using IDataArrayWkPtrType = std::weak_ptr<IDataArray>;

/**
 * @brief The FindDynamicArrayStatistics class. See [Filter documentation](@ref finddynamicarraystatistics) for details.
 */
class FindDynamicArrayStatistics : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = FindDynamicArrayStatistics;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<FindDynamicArrayStatistics> New();

  /**
   * @brief Returns the name of the class for FindDynamicArrayStatistics
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for FindDynamicArrayStatistics
   */
  static QString ClassName();

  ~FindDynamicArrayStatistics() override;

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
   * @brief Setter property for LengthArrayPath
   */
  void setLengthArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for LengthArrayPath
   * @return Value of LengthArrayPath
   */
  DataArrayPath getLengthArrayPath() const;

  Q_PROPERTY(DataArrayPath LengthArrayPath READ getLengthArrayPath WRITE setLengthArrayPath)

  /**
   * @brief Setter property for MinimumArrayPath
   */
  void setMinimumArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MinimumArrayPath
   * @return Value of MinimumArrayPath
   */
  DataArrayPath getMinimumArrayPath() const;

  Q_PROPERTY(DataArrayPath MinimumArrayPath READ getMinimumArrayPath WRITE setMinimumArrayPath)

  /**
   * @brief Setter property for MaximumArrayPath
   */
  void setMaximumArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MaximumArrayPath
   * @return Value of MaximumArrayPath
   */
  DataArrayPath getMaximumArrayPath() const;

  Q_PROPERTY(DataArrayPath MaximumArrayPath READ getMaximumArrayPath WRITE setMaximumArrayPath)

  /**
   * @brief Setter property for MeanArrayPath
   */
  void setMeanArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MeanArrayPath
   * @return Value of MeanArrayPath
   */
  DataArrayPath getMeanArrayPath() const;

  Q_PROPERTY(DataArrayPath MeanArrayPath READ getMeanArrayPath WRITE setMeanArrayPath)

  /**
   * @brief Setter property for MedianArrayPath
   */
  void setMedianArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MedianArrayPath
   * @return Value of MedianArrayPath
   */
  DataArrayPath getMedianArrayPath() const;

  Q_PROPERTY(DataArrayPath MedianArrayPath READ getMedianArrayPath WRITE setMedianArrayPath)

  /**
   * @brief Setter property for StdDeviationArrayPath
   */
  void setStdDeviationArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for StdDeviationArrayPath
   * @return Value of StdDeviationArrayPath
   */
  DataArrayPath getStdDeviationArrayPath() const;

  Q_PROPERTY(DataArrayPath StdDeviationArrayPath READ getStdDeviationArrayPath WRITE setStdDeviationArrayPath)

  /**
   * @brief Setter property for SummationArrayPath
   */
  void setSummationArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for SummationArrayPath
   * @return Value of SummationArrayPath
   */
  DataArrayPath getSummationArrayPath() const;

  Q_PROPERTY(DataArrayPath SummationArrayPath READ getSummationArrayPath WRITE setSummationArrayPath)

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

signals:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter* filter);

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
  /**
   * @brief createCompatibleArrays Creates the output statistics arrays with compatible types to the
   * incoming NeighborList
   * @param dataArrayPaths Vector used to store paths for later tuple checking
   */
  template <typename T> void createCompatibleArrays(QVector<DataArrayPath>& dataArrayPaths)
  {
    std::vector<size_t> cDims(1, 1);

    if(m_FindMin)
    {
      m_MinimumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMinimumArrayPath(), 0,
                                                                                                            cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMinimumArrayPath());
      }
    }

    if(m_FindMax)
    {
      m_MaximumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMaximumArrayPath(), 0,
                                                                                                            cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMaximumArrayPath());
      }
    }

    if(m_FindMean)
    {
      m_MeanPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMeanArrayPath(), 0,
                                                                                                         cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMeanArrayPath());
      }
    }

    if(m_FindMedian)
    {
      m_MedianPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getMedianArrayPath(), 0,
                                                                                                           cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getMedianArrayPath());
      }
    }

    if(m_FindStdDeviation)
    {
      m_StandardDeviationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(
          this, getStdDeviationArrayPath(), 0, cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getStdDeviationArrayPath());
      }
    }

    if(m_FindSummation)
    {
      m_SummationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, getSummationArrayPath(), 0,
                                                                                                              cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
      if(getErrorCode() >= 0)
      {
        dataArrayPaths.push_back(getSummationArrayPath());
      }
    }
  }

  FindDynamicArrayStatistics();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  IDataArrayWkPtrType m_MinimumPtr;
  void* m_Minimum = nullptr;
  IDataArrayWkPtrType m_MaximumPtr;
  void* m_Maximum = nullptr;
  IDataArrayWkPtrType m_MeanPtr;
  void* m_Mean = nullptr;
  IDataArrayWkPtrType m_MedianPtr;
  void* m_Median = nullptr;
  IDataArrayWkPtrType m_StandardDeviationPtr;
  void* m_StandardDeviation = nullptr;
  IDataArrayWkPtrType m_SummationPtr;
  void* m_Summation = nullptr;
  IDataArrayWkPtrType m_InputArrayPtr;
  void* m_InputArray = nullptr;

  std::weak_ptr<DataArray<int32_t>> m_LengthPtr;
  int32_t* m_Length = nullptr;

  bool m_FindLength = {};
  bool m_FindMin = {};
  bool m_FindMax = {};
  bool m_FindMean = {};
  bool m_FindMedian = {};
  bool m_FindStdDeviation = {};
  bool m_FindSummation = {};
  DataArrayPath m_LengthArrayPath = {};
  DataArrayPath m_MinimumArrayPath = {};
  DataArrayPath m_MaximumArrayPath = {};
  DataArrayPath m_MeanArrayPath = {};
  DataArrayPath m_MedianArrayPath = {};
  DataArrayPath m_StdDeviationArrayPath = {};
  DataArrayPath m_SummationArrayPath = {};
  DataArrayPath m_SelectedArrayPath = {};

  FindDynamicArrayStatistics(const FindDynamicArrayStatistics&) = delete; // Copy Constructor Not Implemented
  FindDynamicArrayStatistics(FindDynamicArrayStatistics&&) = delete;      // Move Constructor Not Implemented
  void operator=(const FindDynamicArrayStatistics&) = delete;             // Operator '=' Not Implemented
};

#endif /* _FindDynamicArrayStatistics_H_ */
