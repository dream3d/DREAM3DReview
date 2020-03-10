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
#ifndef _interpolatemeshtoregulargrid_h_
#define _interpolatemeshtoregulargrid_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The InterpolateMeshToRegularGrid class. See [Filter documentation](@ref interpolatemeshtoregulargrid) for details.
 */
class InterpolateMeshToRegularGrid : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = InterpolateMeshToRegularGrid;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<InterpolateMeshToRegularGrid> New();

  /**
   * @brief Returns the name of the class for InterpolateMeshToRegularGrid
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for InterpolateMeshToRegularGrid
   */
  static QString ClassName();

  ~InterpolateMeshToRegularGrid() override;

  /**
   * @brief Setter property for SelectedDataContainerName
   */
  void setSelectedDataContainerName(const QString& value);
  /**
   * @brief Getter property for SelectedDataContainerName
   * @return Value of SelectedDataContainerName
   */
  QString getSelectedDataContainerName() const;

  Q_PROPERTY(QString SelectedDataContainerName READ getSelectedDataContainerName WRITE setSelectedDataContainerName)

  /**
   * @brief Setter property for InterpolatedDataContainerName
   */
  void setInterpolatedDataContainerName(const QString& value);
  /**
   * @brief Getter property for InterpolatedDataContainerName
   * @return Value of InterpolatedDataContainerName
   */
  QString getInterpolatedDataContainerName() const;

  Q_PROPERTY(QString InterpolatedDataContainerName READ getInterpolatedDataContainerName WRITE setInterpolatedDataContainerName)

  /**
   * @brief Setter property for InterpolatedAttributeMatrixName
   */
  void setInterpolatedAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for InterpolatedAttributeMatrixName
   * @return Value of InterpolatedAttributeMatrixName
   */
  QString getInterpolatedAttributeMatrixName() const;

  Q_PROPERTY(QString InterpolatedAttributeMatrixName READ getInterpolatedAttributeMatrixName WRITE setInterpolatedAttributeMatrixName)

  /**
   * @brief Setter property for ScaleOrSpecifyNumCells
   */
  void setScaleOrSpecifyNumCells(int value);
  /**
   * @brief Getter property for ScaleOrSpecifyNumCells
   * @return Value of ScaleOrSpecifyNumCells
   */
  int getScaleOrSpecifyNumCells() const;

  Q_PROPERTY(int ScaleOrSpecifyNumCells READ getScaleOrSpecifyNumCells WRITE setScaleOrSpecifyNumCells)

  /**
   * @brief Setter property for SetXDimension
   */
  void setSetXDimension(int value);
  /**
   * @brief Getter property for SetXDimension
   * @return Value of SetXDimension
   */
  int getSetXDimension() const;

  Q_PROPERTY(int SetXDimension READ getSetXDimension WRITE setSetXDimension)

  /**
   * @brief Setter property for SetYDimension
   */
  void setSetYDimension(int value);
  /**
   * @brief Getter property for SetYDimension
   * @return Value of SetYDimension
   */
  int getSetYDimension() const;

  Q_PROPERTY(int SetYDimension READ getSetYDimension WRITE setSetYDimension)

  /**
   * @brief Setter property for ScaleFactorNumCells
   */
  void setScaleFactorNumCells(int value);
  /**
   * @brief Getter property for ScaleFactorNumCells
   * @return Value of ScaleFactorNumCells
   */
  int getScaleFactorNumCells() const;

  Q_PROPERTY(int ScaleFactorNumCells READ getScaleFactorNumCells WRITE setScaleFactorNumCells)

  /**
   * @brief Setter property for OutsideMeshIdentifier
   */
  void setOutsideMeshIdentifier(int value);
  /**
   * @brief Getter property for OutsideMeshIdentifier
   * @return Value of OutsideMeshIdentifier
   */
  int getOutsideMeshIdentifier() const;

  Q_PROPERTY(int OutsideMeshIdentifier READ getOutsideMeshIdentifier WRITE setOutsideMeshIdentifier)

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
  InterpolateMeshToRegularGrid();

  /**
   * @brief createRegularGrid Creates the structured rectilinear grid to interpolate the mesh
   * data onto
   */
  void createRegularGrid();

  /**
   * @brief interpolateMeshPointsToGrid Performs the interpolation between the selected mesh and
   * the created grid
   * @param inputMesh Source mesh
   * @param interpolatedGrid Destination grid
   */
  void interpolateMeshPointsToGrid(DataContainer::Pointer inputMesh, DataContainer::Pointer interpolatedGrid);
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  QString m_SelectedDataContainerName = {};
  QString m_InterpolatedDataContainerName = {};
  QString m_InterpolatedAttributeMatrixName = {};
  int m_ScaleOrSpecifyNumCells = {};
  int m_SetXDimension = {};
  int m_SetYDimension = {};
  int m_ScaleFactorNumCells = {};
  int m_OutsideMeshIdentifier = {};

  QList<QString> m_AttrMatList;
  QMap<QString, QList<QString>> m_AttrArrayMap;
  std::vector<float> m_MeshMinExtents;
  std::vector<float> m_MeshMaxExtents;
  std::vector<bool> m_InsideMesh;
  std::vector<size_t> m_InterpolatedIndex;

  InterpolateMeshToRegularGrid(const InterpolateMeshToRegularGrid&) = delete; // Copy Constructor Not Implemented
  InterpolateMeshToRegularGrid(InterpolateMeshToRegularGrid&&) = delete;      // Move Constructor Not Implemented
  void operator=(const InterpolateMeshToRegularGrid&) = delete;               // Operator '=' Not Implemented
};

#endif /* _InterpolateMeshToRegularGrid_H_ */
