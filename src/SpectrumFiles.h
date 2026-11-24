/******************************************************************************
  Copyright 2015-2017 Matthew The <matthew.the@scilifelab.se>
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

 ******************************************************************************/

#ifndef QUANDENSER_SPECTRUMFILES_H_
#define QUANDENSER_SPECTRUMFILES_H_

#include <vector>
#include <sstream>

#include "maracluster/src/SpectrumFiles.h"
#include "maracluster/src/SpectrumHandler.h"
#include "maracluster/src/MassChargeCandidate.h"
#include "maracluster/src/ScanId.h"

#include "Globals.h"
#include "DinosaurFeature.h"
#include "SpectrumToPrecursorMap.h"

namespace quandenser {

/* Ties each MS2 spectrum to its precursor candidates: the features of its
   run whose m/z lies in the isolation window and whose retention-time span
   covers the spectrum.

   precursorIonMobilityTol > 0 also requires a feature's ion mobility (1/K0)
   to lie within that tolerance of the spectrum's, where both are known.

   instrumentPrecursorRTimeTol >= 0 adds the precursor the instrument
   selected when no candidate has its charge and m/z: the feature that has
   them and whose retention-time span lies closest to the spectrum, at most
   this many minutes away, or, when there is none, the instrument's
   precursor itself, without a feature. A negative value adds neither. */
class SpectrumFiles : public maracluster::SpectrumFiles {
 public:
  SpectrumFiles(const std::string& precMzFileFolder,
               const int chargeUncertainty,
               const std::vector<DinosaurFeatureList>& featureLists,
               SpectrumToPrecursorMap& specToPrecMap,
               const float precursorIonMobilityTol = 0.0f,
               const float instrumentPrecursorRTimeTol = -1.0f) :
    maracluster::SpectrumFiles(precMzFileFolder, chargeUncertainty), featureLists_(featureLists), spectrumToPrecursorMap_(specToPrecMap),
    precursorIonMobilityTol_(precursorIonMobilityTol),
    instrumentPrecursorRTimeTol_(instrumentPrecursorRTimeTol) {}

  /* m/z tolerance for a feature or candidate to count as the instrument's precursor */
  static const float kInstrumentPrecursorPpmTol;

 protected:
  void getMassChargeCandidates(pwiz::msdata::SpectrumPtr s,
    std::vector<maracluster::MassChargeCandidate>& mccs, maracluster::ScanId scanId);

  static float getIonMobility(pwiz::msdata::SpectrumPtr s);
  static int getInstrumentCharge(pwiz::msdata::SpectrumPtr s);
  static bool isInstrumentPrecursor(int charge, float precMz,
      int instrumentCharge, float instrumentPrecMz);

  bool ionMobilityAgrees(const DinosaurFeature& ft, float ionMobility) const;
  const DinosaurFeature* findInstrumentPrecursorFeature(
      const DinosaurFeatureList& features, int charge, float precMz,
      float rTime, float ionMobility) const;

  const std::vector<DinosaurFeatureList>& featureLists_;
  SpectrumToPrecursorMap& spectrumToPrecursorMap_;
  float precursorIonMobilityTol_;
  float instrumentPrecursorRTimeTol_;
};

} /* namespace quandenser */

#endif /* QUANDENSER_SPECTRUMFILES_H_ */
