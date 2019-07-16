#ifndef _itktransformhelpers_h_
#define _itktransformhelpers_h_


#include "H5Support/H5Lite.h"
#include "H5Support/H5Utilities.h"
#include "H5Support/H5ScopedSentinel.h"

#include "SIMPLib/Geometry/TransformContainer.h"
#include "SIMPLib/ITK/itkInPlaceImageToDream3DDataFilter.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"



struct ITKTransformHelpers
{
	typedef itk::Dream3DImage<double, 2> ImageType;


	ITKTransformHelpers(hid_t fileID, std::string groupName)
	{
		ID = fileID;
		GN = groupName;
	}

	hid_t ID;
	std::string GN; 
	ImageType::SizeType FixedSize;
	ImageType::PointType FixedOrigin;
	ImageType::SpacingType FixedSpacing;
	ImageType::DirectionType FixedDirection;

	std::vector<double> sizeVector;
	std::vector<double> originVector;
	std::vector<double> spacingVector;
	std::vector<double> directionVector;

	QString sizeString = "fixedSize";
	QString originString = "fixedOrigin";
	QString spacingString = "fixedSpacing";
	QString directionString = "fixedDirection";

	TransformContainer::Pointer transform = TransformContainer::New(); 

	int readTransform()
	{
		size_t err = transform->readTransformContainerFromHDF5(ID, false, GN);

		return err;
	}

	int readFixedInformation()
	{
		hid_t groupID = H5Gopen(ID, GN.c_str(), H5P_DEFAULT);
		if (groupID < 0)
		{
			return -1;
		}
		H5ScopedGroupSentinel gSentinel(&groupID, false);

		herr_t err = H5Lite::readVectorDataset(groupID, sizeString.toLatin1().data(), sizeVector);
		if (err < 0)
		{
			return err;
		}

		for (size_t i = 0; i < sizeVector.size(); i++)
		{
			FixedSize.SetElement(i, sizeVector[i]);
		}

		err = H5Lite::readVectorDataset(groupID, originString.toLatin1().data(), originVector);
		if (err < 0)
		{
			return err;
		}

		for (size_t i = 0; i < originVector.size(); i++)
		{
			FixedOrigin.SetElement(i, originVector[i]);
		}

		err = H5Lite::readVectorDataset(groupID, spacingString.toLatin1().data(), spacingVector);
		if (err < 0)
		{
			return err;
		}

		for (size_t i = 0; i < spacingVector.size(); i++)
		{
			FixedSpacing.SetElement(i, spacingVector[i]);
		}

		err = H5Lite::readVectorDataset(groupID, directionString.toLatin1().data(), directionVector);
		if (err < 0)
		{
			return err;
		}
		if (directionVector.size() == 4)
		{
			FixedDirection(0, 0) = directionVector[0]; 
			FixedDirection(0, 1) = directionVector[1];
			FixedDirection(1, 0) = directionVector[2];
			FixedDirection(1, 1) = directionVector[3];
		}
		else 
		{
			FixedDirection.SetIdentity();
		}

		return 0;
	}

	



	/*typedef itk::Dream3DImage<double, 2> ImageType;
	QString sizeString = "fixedSize";
	QString originString = "fixedOrigin";
	QString spacingString = "fixedSpacing";
	QString directionString = "fixedDirection";





	int readFixedInformation(hid_t fileID, std::string groupName)
	{
		hid_t groupID = H5Gopen(fileID, groupName.c_str(), H5P_DEFAULT);
		if (groupID < 0)
		{
			return -1;
		}
		H5ScopedGroupSentinel gSentinel(&groupID, false);

		herr_t err = H5Lite::readVectorDataset(groupID, sizeString.toLatin1().data(), m_sizeVector);
		if (err < 0)
		{
			return err;
		}

		for (size_t i = 0; i < m_sizeVector.size(); i++)
		{
			m_FixedSize.SetElement(i, m_sizeVector[i]);
		}

		err = H5Lite::readVectorDataset(groupID, sizeString.toLatin1().data(), m_originVector);
		if (err < 0)
		{
			return err;
		}

		for (size_t i = 0; i < m_originVector.size(); i++)
		{
			m_FixedOrigin.SetElement(i, m_originVector[i]);
		}

		err = H5Lite::readVectorDataset(groupID, sizeString.toLatin1().data(), m_spacingVector);
		if (err < 0)
		{
			return err;
		}

		for (size_t i = 0; i < m_spacingVector.size(); i++)
		{
			m_FixedSpacing.SetElement(i, m_spacingVector[i]);
		}


	}

	ImageType::SizeType getFixedSize() 
	{
		return m_FixedSize;
	}

	ImageType::PointType getFixedOrigin()
	{
		return m_FixedOrigin;
	}

	ImageType::SpacingType getFixedSpacing()
	{
		return m_FixedSpacing; 
	}


	TransformContainer::Pointer m_Transform;

	ImageType::SizeType m_FixedSize;
	ImageType::PointType m_FixedOrigin;
	ImageType::SpacingType m_FixedSpacing;
	ImageType::DirectionType m_FixedDirection;

	std::vector<double> m_sizeVector;
	std::vector<double> m_originVector;
	std::vector<double> m_spacingVector;
	std::vector<double> m_directionVector;*/

};

#endif /* _itktransformhelpers_h_ */