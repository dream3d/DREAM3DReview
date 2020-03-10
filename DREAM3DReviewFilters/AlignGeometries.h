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
#ifndef _aligngeometries_h_
#define _aligngeometries_h_

#include <memory>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The AlignGeometries class. See [Filter documentation](@ref aligngeometries) for details.
 */
class AlignGeometries : public AbstractFilter
{
  Q_OBJECT

public:
  using Self = AlignGeometries;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<AlignGeometries> New();

  /**
   * @brief Returns the name of the class for AlignGeometries
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for AlignGeometries
   */
  static QString ClassName();

  ~AlignGeometries() override;

  /**
   * @brief Setter property for MovingGeometry
   */
  void setMovingGeometry(const DataArrayPath& value);
  /**
   * @brief Getter property for MovingGeometry
   * @return Value of MovingGeometry
   */
  DataArrayPath getMovingGeometry() const;

  Q_PROPERTY(DataArrayPath MovingGeometry READ getMovingGeometry WRITE setMovingGeometry)

  /**
   * @brief Setter property for TargetGeometry
   */
  void setTargetGeometry(const DataArrayPath& value);
  /**
   * @brief Getter property for TargetGeometry
   * @return Value of TargetGeometry
   */
  DataArrayPath getTargetGeometry() const;

  Q_PROPERTY(DataArrayPath TargetGeometry READ getTargetGeometry WRITE setTargetGeometry)

  /**
   * @brief Setter property for AlignmentType
   */
  void setAlignmentType(int value);
  /**
   * @brief Getter property for AlignmentType
   * @return Value of AlignmentType
   */
  int getAlignmentType() const;

  Q_PROPERTY(int AlignmentType READ getAlignmentType WRITE setAlignmentType)

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
  AlignGeometries();

  /**
  * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
  */
  void dataCheck() override;

  /**
  * @brief Initializes all the private instance variables.
  */
  void initialize();

private:
  DataArrayPath m_MovingGeometry = {};
  DataArrayPath m_TargetGeometry = {};
  int m_AlignmentType = {};

  AlignGeometries(const AlignGeometries&); // Copy Constructor Not Implemented
  AlignGeometries(AlignGeometries&&);      // Move Constructor Not Implemented
  void operator=(const AlignGeometries&);  // Operator '=' Not Implemented
};

#endif /* _aligngeometries_h_ */
