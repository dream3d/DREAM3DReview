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

#include <QtCore/QString>
#include <QtCore/QVector>

#include "EbsdLib/Core/EbsdLibConstants.h"
#include "EbsdLib/EbsdLib.h"
#include "EbsdLib/Core/EbsdSetGetMacros.h"

#include "MicConstants.h"

/**
 * @class MicPhase MicPhase.h EbsdLib/HEDM/MicPhase.h
 * @brief This class holds all the values for a "Phase" header block in a HEDM file
 *
 * @date Mar 23, 2011
 * @version 1.0
 */
class MicPhase
{
public:
  EBSD_SHARED_POINTERS(MicPhase)
  EBSD_STATIC_NEW_MACRO(MicPhase)
  EBSD_TYPE_MACRO(MicPhase)

  virtual ~MicPhase();

  EBSD_INSTANCE_STRING_PROPERTY(PhaseName)
  EBSD_INSTANCE_PROPERTY(std::vector<float>, LatticeConstants)
  EBSD_INSTANCE_STRING_PROPERTY(BasisAtoms)
  EBSD_INSTANCE_STRING_PROPERTY(Symmetry)
  EBSD_INSTANCE_PROPERTY(QVector<QString>, ZandCoordinates)
  EBSD_INSTANCE_PROPERTY(int, PhaseIndex)

  void parseLatticeConstants(std::string& data);
  void parseLatticeAngles(std::string& data);
  void parseBasisAtoms(std::string& data);
  void parseZandCoordinates(std::string& data);

  void printSelf(std::ostream& stream);

  /**
   * @brief Returns the type of crystal structure for this phase.
   */
  unsigned int determineLaueGroup();

  QString getMaterialName();

protected:
  MicPhase();

public:
  MicPhase(const MicPhase&) = delete;            // Copy Constructor Not Implemented
  MicPhase(MicPhase&&) = delete;                 // Move Constructor Not Implemented
  MicPhase& operator=(const MicPhase&) = delete; // Copy Assignment Not Implemented
  MicPhase& operator=(MicPhase&&) = delete;      // Move Assignment Not Implemented
};
