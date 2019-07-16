

#pragma once

#include <cstdint>
#include <map>

#include <QtCore/QString>

#include "SIMPLib/DataArrays/IDataArray.h"

class EBSDWriterFactory
{
public:
  virtual ~EBSDWriterFactory();

  EBSDWriterFactory(const EBSDWriterFactory&) = delete;
  EBSDWriterFactory& operator=(const EBSDWriterFactory&) = delete;

  static EBSDWriterFactory* Instance();

  /**
   * @brief getWriter
   * @param name
   * @return
   */
  std::function<QString(QString, const IDataArray::Pointer&, uint64_t)> getWriter(const QString& name);

private:
  EBSDWriterFactory();

  /**
   * @brief initializeDataTypes
   */
  void initializeDataTypes();

  static EBSDWriterFactory* self;
  std::map<QString, std::function<QString(QString, const IDataArray::Pointer&, uint64_t)>> m_DataTypes;
};
