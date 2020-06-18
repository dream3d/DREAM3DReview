#include <utility>

#include <utility>

/* ============================================================================
 * Copyright (c) 2010, Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2010, Dr. Michael A. Groeber (US Air Force Research Laboratories
 * All rights reserved.
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include <cstring>

#include <QtCore/QByteArray>
#include <QtCore/QTextStream>

#include "H5Support/H5Lite.h"

#include "EbsdLib/IO/EbsdHeaderEntry.h"
#include "EbsdLib/EbsdLib.h"
#include "EbsdLib/Core/EbsdSetGetMacros.h"

/**
 * @class MicHeaderEntry MicHeaderEntry.h EbsdLib/HEDM/MicHeaderEntry.h
 * @brief Header entry that holds an integer or decimal type value
 *
 * @date Aug 8, 2011
 * @version 1.0
 */
template <typename T>
class MicHeaderEntry : public EbsdHeaderEntry
{

public:
  EBSD_SHARED_POINTERS(MicHeaderEntry<T>)
  HEADERENTRY_NEW_SUPERCLASS(MicHeaderEntry<T>, EbsdHeaderEntry)
  HEADERENTRY_NEW_SUPERCLASS_VALUE(MicHeaderEntry<T>, EbsdHeaderEntry)

  ~MicHeaderEntry() override = default;

  std::string getKey() override
  {
    return m_key;
  }

#ifdef EbsdLib_ENABLE_HDF5
  std::string getHDFType() override
  {
    T value = static_cast<T>(0);
    return H5Lite::HDFTypeForPrimitiveAsStr(value);
  }
#endif

  void parseValue(std::string& value) override
  {
    std::stringstream ss(value);
    ss >> m_value;
  }

  void print(std::ostream& out) override
  {
    out << m_key << "  " << m_value << std::endl;
  }

  T getValue()
  {
    return m_value;
  }
  void setValue(T value)
  {
    m_value = value;
  }

protected:
  MicHeaderEntry(std::string key)
  : m_key(std::move(key))
  {
  }

  MicHeaderEntry() = default;

private:
  T m_value = static_cast<T>(0);
  std::string m_key = {""};

public:
  MicHeaderEntry(const MicHeaderEntry&) = delete;            // Copy Constructor Not Implemented
  MicHeaderEntry(MicHeaderEntry&&) = delete;                 // Move Constructor Not Implemented
  MicHeaderEntry& operator=(const MicHeaderEntry&) = delete; // Copy Assignment Not Implemented
  MicHeaderEntry& operator=(MicHeaderEntry&&) = delete;      // Move Assignment Not Implemented
};

/**
 * @class MicStringHeaderEntry MicStringHeaderEntry.h EbsdLib/HEDM/MicHeaderEntry.h
 * @brief Header entry that holds a string type value
 *
 * @date Aug 1, 2011
 * @version 1.0
 */
class MicStringHeaderEntry : public EbsdHeaderEntry
{
public:
  EBSD_SHARED_POINTERS(MicStringHeaderEntry)
  HEADERENTRY_NEW_SUPERCLASS(MicStringHeaderEntry, EbsdHeaderEntry)

  ~MicStringHeaderEntry() override = default;

  std::string getKey() override
  {
    return m_key;
  }
  std::string getHDFType() override
  {
    return "H5T_STRING";
  }

  void parseValue(QByteArray& value)
  {
    m_value = std::string(value.trimmed());
  }

  void print(std::ostream& out) override
  {
    out << m_key << "  " << m_value << std::endl;
  }

  std::string getValue()
  {
    return m_value;
  }
  void setValue(const std::string& value)
  {
    m_value = value;
  }

  void parseValue(std::string& value)
  {
    m_value = value;
  }

protected:
  MicStringHeaderEntry(std::string key)
  : m_key(std::move(key))
  {
  }

  MicStringHeaderEntry() = default;

private:
  std::string m_value = {""};
  std::string m_key = {""};

public:
  MicStringHeaderEntry(const MicStringHeaderEntry&) = delete;            // Copy Constructor Not Implemented
  MicStringHeaderEntry(MicStringHeaderEntry&&) = delete;                 // Move Constructor Not Implemented
  MicStringHeaderEntry& operator=(const MicStringHeaderEntry&) = delete; // Copy Assignment Not Implemented
  MicStringHeaderEntry& operator=(MicStringHeaderEntry&&) = delete;      // Move Assignment Not Implemented
};
