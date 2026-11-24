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

#include <algorithm>
#include <cmath>

#include "SpectrumFiles.h"

namespace quandenser {

const float SpectrumFiles::kInstrumentPrecursorPpmTol = 20.0f;

void SpectrumFiles::getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<maracluster::MassChargeCandidate>& mccs, maracluster::ScanId scanId) {
  //maracluster::SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty_);
  //mccs.clear();
  float precMz = maracluster::SpectrumHandler::getPrecMz(s);
  int instrumentCharge = (instrumentPrecursorRTimeTol_ >= 0.0f) ? getInstrumentCharge(s) : 0;
  
  if (!spectrumToPrecursorMap_.isInitialized(scanId)) {
    spectrumToPrecursorMap_.setInitialized(scanId);
    float rTime = maracluster::SpectrumHandler::getRetentionTime(s);  
    float ionMobility = (precursorIonMobilityTol_ > 0.0f) ? getIonMobility(s) : 0.0f;
    
    float isoWidthLower = s->precursors.at(0).isolationWindow.cvParam(pwiz::cv::MS_isolation_window_lower_offset).valueAs<float>();
    float isoWidthUpper = s->precursors.at(0).isolationWindow.cvParam(pwiz::cv::MS_isolation_window_upper_offset).valueAs<float>();
    
    const DinosaurFeatureList& features = featureLists_.at(scanId.fileIdx);
    bool hasInstrumentPrecursor = false;
    // TODO: try if using a separate vector with only the precMzs is faster
    DinosaurFeatureList::const_iterator lowerBound = features.getPrecMzIterator(precMz - isoWidthLower);
    DinosaurFeatureList::const_iterator upperBound = features.getPrecMzIterator(precMz + isoWidthUpper);
    for (DinosaurFeatureList::const_iterator it = lowerBound; it != upperBound; ++it) {
      if (it->intensity > 0.0 && it->rtStart <= rTime && it->rtEnd >= rTime && 
          ionMobilityAgrees(*it, ionMobility)) {
        DinosaurFeature featureCopy = *it;
        spectrumToPrecursorMap_.addFeature(scanId, featureCopy);
        if (isInstrumentPrecursor(it->charge, it->precMz, instrumentCharge, precMz)) {
          hasInstrumentPrecursor = true;
        }
      }
    }
    
    if (instrumentCharge > 0 && !hasInstrumentPrecursor) {
      const DinosaurFeature* ft = findInstrumentPrecursorFeature(features, 
          instrumentCharge, precMz, rTime, ionMobility);
      if (ft != NULL) {
        DinosaurFeature featureCopy = *ft;
        spectrumToPrecursorMap_.addFeature(scanId, featureCopy);
      }
    }
    
    if (Globals::VERB > 4 && spectrumToPrecursorMap_.getBegin(scanId) == spectrumToPrecursorMap_.getEnd(scanId)) {
      std::cerr << "could not find precursor: fileIdx = " << scanId.fileIdx << ", precmz = " << precMz << ", rtime = " << rTime << std::endl;
    }
  }
  
  bool hasInstrumentPrecursor = false;
  for (DinosaurFeatureList::const_iterator it = spectrumToPrecursorMap_.getBegin(scanId); it != spectrumToPrecursorMap_.getEnd(scanId); ++it) {
    mccs.push_back(maracluster::MassChargeCandidate(it->charge, it->precMz, maracluster::SpectrumHandler::calcMass(it->precMz, it->charge)));
    if (isInstrumentPrecursor(it->charge, it->precMz, instrumentCharge, precMz)) {
      hasInstrumentPrecursor = true;
    }
  }
  
  /* no feature has the instrument's precursor: cluster the spectrum on it 
     all the same; it is not a feature, so it reaches no feature group */
  if (instrumentCharge > 0 && !hasInstrumentPrecursor) {
    mccs.push_back(maracluster::MassChargeCandidate(instrumentCharge, precMz, maracluster::SpectrumHandler::calcMass(precMz, instrumentCharge)));
  }
  
  if (mccs.size() > 1) {
    std::sort(mccs.begin(), mccs.end(), maracluster::MassChargeCandidate::lessChargeMass);
  }
}

/* 1/K0 of the spectrum's precursor, 0 if the file does not give it */
float SpectrumFiles::getIonMobility(pwiz::msdata::SpectrumPtr s) {
  float ionMobility = 0.0f;
  if (s->scanList.scans.size() > 0) {
    ionMobility = s->scanList.scans.front().cvParam(pwiz::cv::MS_inverse_reduced_ion_mobility).valueAs<float>();
  }
  if (ionMobility <= 0.0f && s->precursors.size() > 0 && 
      s->precursors.front().selectedIons.size() > 0) {
    ionMobility = s->precursors.front().selectedIons.front().cvParam(pwiz::cv::MS_inverse_reduced_ion_mobility).valueAs<float>();
  }
  return ionMobility;
}

/* charge the instrument assigned to the precursor, 0 if it gave none */
int SpectrumFiles::getInstrumentCharge(pwiz::msdata::SpectrumPtr s) {
  if (s->precursors.size() > 0 && s->precursors.front().selectedIons.size() > 0) {
    return s->precursors.front().selectedIons.front().cvParam(pwiz::cv::MS_charge_state).valueAs<int>();
  }
  return 0;
}

bool SpectrumFiles::isInstrumentPrecursor(int charge, float precMz, 
    int instrumentCharge, float instrumentPrecMz) {
  return instrumentCharge > 0 && charge == instrumentCharge && 
      std::abs(precMz - instrumentPrecMz) <= instrumentPrecMz * kInstrumentPrecursorPpmTol * 1e-6f;
}

bool SpectrumFiles::ionMobilityAgrees(const DinosaurFeature& ft, 
    float ionMobility) const {
  return precursorIonMobilityTol_ <= 0.0f || ionMobility <= 0.0f || 
      ft.ionMobility <= 0.0f || 
      std::abs(ft.ionMobility - ionMobility) <= precursorIonMobilityTol_;
}

/* The feature with the instrument's precursor charge and m/z whose 
   retention-time span lies closest to the spectrum, at most 
   instrumentPrecursorRTimeTol_ minutes away; on a tie the most intense. */
const DinosaurFeature* SpectrumFiles::findInstrumentPrecursorFeature(
    const DinosaurFeatureList& features, int charge, float precMz, 
    float rTime, float ionMobility) const {
  float mzTol = precMz * kInstrumentPrecursorPpmTol * 1e-6f;
  DinosaurFeatureList::const_iterator lowerBound = features.getPrecMzIterator(precMz - mzTol);
  DinosaurFeatureList::const_iterator upperBound = features.getPrecMzIterator(precMz + mzTol);
  const DinosaurFeature* best = NULL;
  float bestGap = 0.0f;
  for (DinosaurFeatureList::const_iterator it = lowerBound; it != upperBound; ++it) {
    if (it->intensity > 0.0 && isInstrumentPrecursor(it->charge, it->precMz, charge, precMz) &&
        ionMobilityAgrees(*it, ionMobility)) {
      float gap = (std::max)(0.0f, (std::max)(it->rtStart - rTime, rTime - it->rtEnd));
      if (gap <= instrumentPrecursorRTimeTol_ && (best == NULL || gap < bestGap || 
            (gap == bestGap && it->intensity > best->intensity))) {
        best = &*it;
        bestGap = gap;
      }
    }
  }
  return best;
}

} /* namespace quandenser */
