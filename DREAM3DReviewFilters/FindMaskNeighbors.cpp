/* ============================================================================
 * Copyright (c) 2009-2016 BlueQuartz Software, LLC
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
 * Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
 * contributors may be used to endorse or promote products derived from this software
 * without specific prior written permission.
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
 * The code contained herein was partially funded by the followig contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "FindMaskNeighbors.h"

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/Geometry/ImageGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
// enum createdPathID : RenameDataPath::DataID_t
//{
//  DataArrayID30 = 30,
//  DataArrayID31 = 31,
//};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindMaskNeighbors::FindMaskNeighbors()
: m_FeatureIdsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::FeatureIds)
, m_MaskArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellAttributeMatrixName, SIMPL::CellData::Mask)
, m_MaskNeighborsArrayPath(SIMPL::Defaults::ImageDataContainerName, SIMPL::Defaults::CellFeatureAttributeMatrixName, "MaskNeighbors")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindMaskNeighbors::~FindMaskNeighbors() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::RequiredArray));
  {
    DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature IDs", FeatureIdsArrayPath, FilterParameter::RequiredArray, FindMaskNeighbors, req));

    DataArraySelectionFilterParameter::RequirementType req2 = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask Array", MaskArrayPath, FilterParameter::RequiredArray, FindMaskNeighbors, req2));
  }
  parameters.push_back(SeparatorFilterParameter::New("Cell Feature Data", FilterParameter::CreatedArray));
  {
    DataArrayCreationFilterParameter::RequirementType req = DataArrayCreationFilterParameter::CreateRequirement(AttributeMatrix::Category::Feature);
    req.dcGeometryTypes = IGeometry::Types(1, IGeometry::Type::Image);
    parameters.push_back(SIMPL_NEW_DA_CREATION_FP("Mask Neighbors", MaskNeighborsArrayPath, FilterParameter::CreatedArray, FindMaskNeighbors, req));
  }
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setMaskNeighborsArrayPath(reader->readDataArrayPath("MaskNeighborsArrayPath", getMaskNeighborsArrayPath()));
  setFeatureIdsArrayPath(reader->readDataArrayPath("FeatureIdsArrayPath", getFeatureIdsArrayPath()));
  setMaskArrayPath(reader->readDataArrayPath("MaskArrayPath", getMaskArrayPath()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getFeatureIdsArrayPath().getDataContainerName());

  QVector<size_t> cDims(1, 1);
  m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getFeatureIdsArrayPath(),
                                                                                                        cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_FeatureIdsPtr.lock())                                                                         /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_MaskPtr =
      getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if(nullptr != m_MaskPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_Mask = m_MaskPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */

  m_MaskNeighborsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<bool>, AbstractFilter, bool>(this, getMaskNeighborsArrayPath(), false, cDims, "");
  if(nullptr != m_MaskNeighborsPtr.lock()) /* Validate the Weak Pointer wraps a non-nullptr pointer to a DataArray<T> object */
  {
    m_MaskNeighbors = m_MaskNeighborsPtr.lock()->getPointer(0);
  } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheck();
  emit preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool FindMaskNeighbors::validNeighbor(int64_t dims[3], int64_t neighborhood[78], size_t index, int64_t x, int64_t y, int64_t z)
{
  int64_t modX = x + neighborhood[3 * index + 0];
  int64_t modY = y + neighborhood[3 * index + 1];
  int64_t modZ = z + neighborhood[3 * index + 2];

  return (modX >= 0 && modX < dims[0]) && (modY >= 0 && modY < dims[1]) && (modZ >= 0 && modZ < dims[2]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
bool FindMaskNeighbors::validNeighbor2D(int64_t dims[2], int64_t neighborhood[16], size_t index, int64_t x, int64_t y)
{
  int64_t modX = x + neighborhood[2 * index + 0];
  int64_t modY = y + neighborhood[2 * index + 1];

  return (modX >= 0 && modX < dims[0]) && (modY >= 0 && modY < dims[1]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::find_surfacefeatures()
{
  int64_t neighborhood[78] = { 1,  0, 0,  -1, 0, 0, 0, 1, 0,  0, -1, 0, 0, 0,  1,  0, 0, -1, 1, 1, 0,  -1, 1,  0, 1, -1, 0,  -1, -1, 0, 1,  0, 1,  1,  0,  -1, -1, 0,  1,
                            -1, 0, -1, 0,  1, 1, 0, 1, -1, 0, -1, 1, 0, -1, -1, 1, 1, 1,  1, 1, -1, 1,  -1, 1, 1, -1, -1, -1, 1,  1, -1, 1, -1, -1, -1, 1,  -1, -1, -1 };
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());

  int64_t xPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getXPoints());
  int64_t yPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getYPoints());
  int64_t zPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getZPoints());

  int64_t dims[3] = {xPoints, yPoints, zPoints};


  int64_t zStride = 0, yStride = 0;
  for(int64_t i = 0; i < zPoints; i++)
  {
    zStride = i * xPoints * yPoints;
    for(int64_t j = 0; j < yPoints; j++)
    {
      yStride = j * xPoints;
      for(int64_t k = 0; k < xPoints; k++)
      {
        int32_t gnum = m_FeatureIds[zStride + yStride + k];
        if(!m_MaskNeighbors[gnum])
        {
          if (!m_Mask[yStride + k])
          {

            std::vector<int> validNeighborIndicies;

            for (size_t n = 0; n < 8; n++)
            {
              if (validNeighbor(dims, neighborhood, n, k, j, i))
              {
                int64_t neighbor = ((i + neighborhood[3 * n + 2]) * dims[1] * dims[0]) + ((j + neighborhood[3 * n + 1]) * dims[0]) + (k + neighborhood[3 * n + 0]);
                if (m_Mask[neighbor])
                {
                  m_MaskNeighbors[gnum] = true;
                }
              }
            }
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::find_surfacefeatures2D()
{
  int64_t neighborhood[16] = { 1, 0, -1, 0, 0, 1, 0, -1, 1, 1, -1, 1, 1, -1, -1, -1 };


  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());

  int64_t xPoints = 0, yPoints = 0;

  if (m->getGeometryAs<ImageGeom>()->getXPoints() == 1)
  {
    xPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getYPoints());
    yPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getZPoints());
  }
  if (m->getGeometryAs<ImageGeom>()->getYPoints() == 1)
  {
    xPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getXPoints());
    yPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getZPoints());
  }
  if (m->getGeometryAs<ImageGeom>()->getZPoints() == 1)
  {
    xPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getXPoints());
    yPoints = static_cast<int64_t>(m->getGeometryAs<ImageGeom>()->getYPoints());
  }

  int64_t dims[2] = { xPoints, yPoints };

  int64_t yStride;
  for (int64_t j = 0; j < yPoints; j++)
  {
    yStride = j * xPoints;
    for (int64_t k = 0; k < xPoints; k++)
    {
      int32_t gnum = m_FeatureIds[yStride + k];

      if (gnum > 0)
      {
        if (!m_MaskNeighbors[gnum])
        {
          if (!m_Mask[yStride + k])
          {

            std::vector<int> validNeighborIndicies;

            for (size_t n = 0; n < 8; n++)
            {
              if (validNeighbor2D(dims, neighborhood, n, j, k))
              {
                int64_t neighbor = ((j + neighborhood[2 * n + 1]) * dims[0]) + (k + neighborhood[2 * n + 0]);
                if (m_Mask[neighbor])
                {
                  m_MaskNeighbors[gnum] = true;
                }
              }
            }
          }

        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindMaskNeighbors::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getFeatureIdsArrayPath().getDataContainerName());

  if(m->getGeometryAs<ImageGeom>()->getXPoints() > 1 && m->getGeometryAs<ImageGeom>()->getYPoints() > 1 && m->getGeometryAs<ImageGeom>()->getZPoints() > 1)
  {
    find_surfacefeatures();
  }
  if(m->getGeometryAs<ImageGeom>()->getXPoints() == 1 || m->getGeometryAs<ImageGeom>()->getYPoints() == 1 || m->getGeometryAs<ImageGeom>()->getZPoints() == 1)
  {
    find_surfacefeatures2D();
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindMaskNeighbors::newFilterInstance(bool copyFilterParameters) const
{
  FindMaskNeighbors::Pointer filter = FindMaskNeighbors::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMaskNeighbors::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMaskNeighbors::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMaskNeighbors::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMaskNeighbors::getGroupName() const
{
  return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid FindMaskNeighbors::getUuid()
{
  return QUuid("{9dc209a4-c6d2-50d0-bcba-ce765aeb7d4f}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMaskNeighbors::getSubGroupName() const
{
  return SIMPL::FilterSubGroups::SpatialFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindMaskNeighbors::getHumanLabel() const
{
  return "Find Mask Neighbors";
}
