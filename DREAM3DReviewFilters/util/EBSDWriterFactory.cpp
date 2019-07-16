#include "EBSDWriterFactory.h"

#include <cmath>
#include <fstream>

#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataArrays/StringDataArray.h"

namespace EBSDWriterHelper
{

template <typename T>
inline QString WriteArrayToQStringLine(QString rowData, const IDataArray::Pointer& ptr, uint64_t index)
{
  typename DataArray<T>::Pointer data = std::dynamic_pointer_cast<DataArray<T>>(ptr);
  int numComps = data->getNumberOfComponents();

  for(int i = 0; i < numComps; i++)
  {
    T value = data->getValue(numComps * index + i);
    QString arrayElement = QString::number(value);
    rowData.append(arrayElement);
    if(i != (numComps - 1))
    {
      rowData.append("	");
    }
  }

  return rowData;
}

template <>
inline QString WriteArrayToQStringLine<int32_t>(QString rowData, const IDataArray::Pointer& ptr, uint64_t index)
{
  DataArray<int32_t>::Pointer data = std::dynamic_pointer_cast<DataArray<int32_t>>(ptr);
  int numComps = data->getNumberOfComponents();

  for(int i = 0; i < numComps; i++)
  {
    int32_t value = data->getValue(numComps * index + i);
    QString arrayElement = QString::number(value);
    rowData.append(arrayElement);
    if(i != (numComps - 1))
    {
      rowData.append("	");
    }
  }

  return rowData;
}

template <>
inline QString WriteArrayToQStringLine<float>(QString rowData, const IDataArray::Pointer& ptr, uint64_t index)
{
  DataArray<float>::Pointer data = std::dynamic_pointer_cast<DataArray<float>>(ptr);
  int numComps = data->getNumberOfComponents();

  for(int i = 0; i < numComps; i++)
  {
    float value = data->getValue(numComps * index + i);
    QString arrayElement = QString::number(value, 'f', 4);
    rowData.append(arrayElement);
    if(i != (numComps - 1))
    {
      rowData.append("	");
    }
  }

  return rowData;
}

template <typename T>
inline std::function<QString(QString, IDataArray::Pointer, uint64_t)> ArrayWriterFactory()
{
  return WriteArrayToQStringLine<T>;
}
} // namespace EBSDWriterHelper

EBSDWriterFactory* EBSDWriterFactory::self = nullptr;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EBSDWriterFactory::EBSDWriterFactory()
{
  initializeDataTypes();
  self = this;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EBSDWriterFactory::~EBSDWriterFactory() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
EBSDWriterFactory* EBSDWriterFactory::Instance()
{
  if(self == nullptr)
  {
    self = new EBSDWriterFactory();
  }
  return self;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::function<QString(QString, const IDataArray::Pointer&, uint64_t)> EBSDWriterFactory::getWriter(const QString& name)
{
  std::function<QString(QString, IDataArray::Pointer, uint64_t)> writer = nullptr;
  try
  {
    writer = m_DataTypes.at(name);
  } catch(const std::out_of_range&)
  {
    // std::string info("Data type index from file: " + std::to_string(index) + "\n" + "Supported data type indices: \n" + getListOfSupportedDataTypes());
    // throw FatalTDMSException(TDMSExceptionMessages::UnsupportedDataType, info);
  }

  return writer;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void EBSDWriterFactory::initializeDataTypes()
{
  m_DataTypes["Phases"] = EBSDWriterHelper::ArrayWriterFactory<int>();
  m_DataTypes["Error"] = EBSDWriterHelper::ArrayWriterFactory<int>();
  m_DataTypes["MAD"] = EBSDWriterHelper::ArrayWriterFactory<float>();
  m_DataTypes["BC"] = EBSDWriterHelper::ArrayWriterFactory<int>();
  m_DataTypes["EulerAngles"] = EBSDWriterHelper::ArrayWriterFactory<float>();
  m_DataTypes["Bands"] = EBSDWriterHelper::ArrayWriterFactory<int>();
  m_DataTypes["BS"] = EBSDWriterHelper::ArrayWriterFactory<int>();
  m_DataTypes["X"] = EBSDWriterHelper::ArrayWriterFactory<float>();
  m_DataTypes["Y"] = EBSDWriterHelper::ArrayWriterFactory<float>();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
