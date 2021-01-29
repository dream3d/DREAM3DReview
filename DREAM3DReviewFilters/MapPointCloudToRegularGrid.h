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
#include "SIMPLib/FilterParameters/IntVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/IGeometry.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The MapPointCloudToRegularGrid class. See [Filter documentation](@ref mappointcloudtoregulargrid) for details.
 */
class DREAM3DReview_EXPORT MapPointCloudToRegularGrid : public AbstractFilter
{
  Q_OBJECT

  PYB11_BEGIN_BINDINGS(MapPointCloudToRegularGrid SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(MapPointCloudToRegularGrid)
  PYB11_FILTER_NEW_MACRO(MapPointCloudToRegularGrid)
  PYB11_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)
  PYB11_PROPERTY(DataArrayPath CreatedImageDataContainerName READ getCreatedImageDataContainerName WRITE setCreatedImageDataContainerName)
  PYB11_PROPERTY(DataArrayPath ImageDataContainerPath READ getImageDataContainerPath WRITE setImageDataContainerPath)
  PYB11_PROPERTY(DataArrayPath VoxelIndicesArrayPath READ getVoxelIndicesArrayPath WRITE setVoxelIndicesArrayPath)
  PYB11_PROPERTY(IntVec3Type GridDimensions READ getGridDimensions WRITE setGridDimensions)
  PYB11_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)
  PYB11_PROPERTY(int SamplingGridType READ getSamplingGridType WRITE setSamplingGridType)
  PYB11_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)
  PYB11_END_BINDINGS()

public:
  using Self = MapPointCloudToRegularGrid;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<MapPointCloudToRegularGrid> New();

  /**
   * @brief Returns the name of the class for MapPointCloudToRegularGrid
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for MapPointCloudToRegularGrid
   */
  static QString ClassName();

  ~MapPointCloudToRegularGrid() override;

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
   * @brief Setter property for CreatedImageDataContainerName
   */
  void setCreatedImageDataContainerName(const DataArrayPath& value);
  /**
   * @brief Getter property for CreatedImageDataContainerName
   * @return Value of CreatedImageDataContainerName
   */
  DataArrayPath getCreatedImageDataContainerName() const;
  Q_PROPERTY(DataArrayPath CreatedImageDataContainerName READ getCreatedImageDataContainerName WRITE setCreatedImageDataContainerName)

  /**
   * @brief Setter property for ImageDataContainerPath
   */
  void setImageDataContainerPath(const DataArrayPath& value);
  /**
   * @brief Getter property for ImageDataContainerPath
   * @return Value of ImageDataContainerPath
   */
  DataArrayPath getImageDataContainerPath() const;
  Q_PROPERTY(DataArrayPath ImageDataContainerPath READ getImageDataContainerPath WRITE setImageDataContainerPath)

  /**
   * @brief Setter property for VoxelIndicesArrayPath
   */
  void setVoxelIndicesArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for VoxelIndicesArrayPath
   * @return Value of VoxelIndicesArrayPath
   */
  DataArrayPath getVoxelIndicesArrayPath() const;
  Q_PROPERTY(DataArrayPath VoxelIndicesArrayPath READ getVoxelIndicesArrayPath WRITE setVoxelIndicesArrayPath)

  /**
   * @brief Setter property for GridDimensions
   */
  void setGridDimensions(const IntVec3Type& value);
  /**
   * @brief Getter property for GridDimensions
   * @return Value of GridDimensions
   */
  IntVec3Type getGridDimensions() const;
  Q_PROPERTY(IntVec3Type GridDimensions READ getGridDimensions WRITE setGridDimensions)

  /**
   * @brief Setter property for UseMask
   */
  void setUseMask(bool value);
  /**
   * @brief Getter property for UseMask
   * @return Value of UseMask
   */
  bool getUseMask() const;
  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  /**
   * @brief Setter property for SamplingGridType
   */
  void setSamplingGridType(int value);
  /**
   * @brief Getter property for SamplingGridType
   * @return Value of SamplingGridType
   */
  int getSamplingGridType() const;
  Q_PROPERTY(int SamplingGridType READ getSamplingGridType WRITE setSamplingGridType)

  /**
   * @brief Setter property for MaskArrayPath
   */
  void setMaskArrayPath(const DataArrayPath& value);
  /**
   * @brief Getter property for MaskArrayPath
   * @return Value of MaskArrayPath
   */
  DataArrayPath getMaskArrayPath() const;
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

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
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

protected:
  MapPointCloudToRegularGrid();

  /**
   * @brief createRegularGrid Creates the structured rectilinear grid to interpolate
   * the point cloud onto
   */
  void createRegularGrid();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<MeshIndexType>> m_VoxelIndicesPtr;
  MeshIndexType* m_VoxelIndices = nullptr;
  std::weak_ptr<DataArray<bool>> m_MaskPtr;
  bool* m_Mask = nullptr;

  DataArrayPath m_DataContainerName = {"", "", ""};
  DataArrayPath m_CreatedImageDataContainerName = {"ImageDataContainer", "", ""};
  DataArrayPath m_ImageDataContainerPath = {"", "", ""};
  DataArrayPath m_VoxelIndicesArrayPath = {"", "", "VoxelIndices"};
  IntVec3Type m_GridDimensions = {10, 10, 10};
  bool m_UseMask = {false};
  int m_SamplingGridType = {0};
  DataArrayPath m_MaskArrayPath = {"", "", ""};

  std::vector<float> m_MeshMinExtents;
  std::vector<float> m_MeshMaxExtents;

public:
  MapPointCloudToRegularGrid(const MapPointCloudToRegularGrid&) = delete;            // Copy Constructor Not Implemented
  MapPointCloudToRegularGrid(MapPointCloudToRegularGrid&&) = delete;                 // Move Constructor Not Implemented
  MapPointCloudToRegularGrid& operator=(const MapPointCloudToRegularGrid&) = delete; // Copy Assignment Not Implemented
  MapPointCloudToRegularGrid& operator=(MapPointCloudToRegularGrid&&) = delete;      // Move Assignment Not Implemented
};
