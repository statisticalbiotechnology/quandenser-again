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
#include <cstdio>
#include <limits>
#include <stdint.h>

#include "PercolatorAdapter.h"

namespace quandenser {

PercolatorAdapter::PercolatorAdapter() : Caller(), setHandler_(0u), 
    subsetUpperLimit_(std::numeric_limits<size_t>::max()), numPsmsAdded_(0u) {}

namespace {
/* splitmix64: a fixed pseudo-random priority per scan number */
size_t scanPriority(int scannr) {
  uint64_t z = static_cast<uint64_t>(static_cast<int64_t>(scannr)) + 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return static_cast<size_t>(z ^ (z >> 31));
}
}
  
void PercolatorAdapter::process() {
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  if (Globals::VERB > 0) {
    cerr << extendedGreeter(startTime);
  }
  
  if (Globals::VERB > 2) {
    std::cerr << "FeatureNames::getNumFeatures(): "<< FeatureNames::getNumFeatures() << endl;
  }
  
  bool subsetOnly = maxPSMs_ > 0u && numPsmsAdded_ > subsetPsms_.size();
  if (maxPSMs_ > 0u) {
    allPsmsStream_.close();
    if (Globals::VERB > 1) {
      std::cerr << "Training on " << subsetPsms_.size() << " of "
                << numPsmsAdded_ << " PSMs" << std::endl;
    }
    /* in the order added, so that with every PSM kept the result is the
       same as without --subset-max-train */
    std::vector<SubsetPsm> subset;
    subset.reserve(subsetPsms_.size());
    for (; !subsetPsms_.empty(); subsetPsms_.pop()) subset.push_back(subsetPsms_.top());
    std::sort(subset.begin(), subset.end(), SubsetPsm::byOrder);
    std::vector<SubsetPsm>::const_iterator it;
    for (it = subset.begin(); it != subset.end(); ++it) {
      setHandler_.getSubsetFromLabel(it->label)->registerPsm(it->psm);
    }
    if (!subsetOnly) remove(allPsmsFN_.c_str());
  }
  
  setHandler_.normalizeFeatures(pNorm_);
  
  // Copy feature data pointers to Scores object
  Scores allScores(useMixMax_);
  allScores.populateWithPSMs(setHandler_);
  
  CrossValidation crossValidation(quickValidation_, reportEachIteration_, 
                                  testFdr_, selectionFdr_, initialSelectionFdr_, selectedCpos_, 
                                  selectedCneg_, numIterations_, useMixMax_,
                                  nestedXvalBins_, trainBestPositive_,
                                  numThreads_, skipNormalizeScores_);
  int firstNumberOfPositives = crossValidation.preIterationSetup(allScores, pCheck_, pNorm_, setHandler_.getFeaturePool());
  if (Globals::VERB > 0) {
    cerr << "Found " << firstNumberOfPositives << " test set positives with q<"
        << testFdr_ << " in initial direction" << endl;
  }
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (Globals::VERB > 1) cerr << "Reading in data and feature calculation took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall clock time." << endl;
  
  if (tabOutputFN_.length() > 0) {
    setHandler_.writeTab(tabOutputFN_, pCheck_);
  }
  
  // Do the SVM training
  crossValidation.train(pNorm_);
  
  if (weightOutputFN_.size() > 0) {
    ofstream weightStream(weightOutputFN_.c_str(), ios::out);
    crossValidation.printAllWeights(weightStream, pNorm_);
    weightStream.close();
  }
  
  // Calculate the final SVM scores and clean up structures
  crossValidation.postIterationProcessing(allScores, pCheck_);
  
  if (subsetOnly) scoreAllPsms(crossValidation, allScores);
  
  // calculate psms level probabilities TDA or TDC
  bool isUniquePeptideRun = false;
  calculatePSMProb(allScores, isUniquePeptideRun, procStart, procStartClock, diff);
  
  std::vector<PSMDescription*>::iterator psmIt;
  for (psmIt = scoredPsms_.begin(); psmIt != scoredPsms_.end(); ++psmIt) {
    PSMDescription::deletePtr(*psmIt);
  }
  std::vector<PSMDescription*>().swap(scoredPsms_);
  
  /* These are all example of why you should not use the singleton pattern... */
  DataSet::resetFeatureNames();
  Normalizer::resetNormalizer();
}
  
void PercolatorAdapter::init(const std::vector<std::pair<std::string, double> >& percolatorFeatures) {
  // read feature names and initial values that are present in feature descriptions
  FeatureNames& featureNames = DataSet::getFeatureNames();
  std::vector<std::pair<std::string, double> >::const_iterator featureIt;
  std::vector<double> init_values(percolatorFeatures.size());
  bool hasDefaultValues = false;
  int i = 0;
  for (featureIt = percolatorFeatures.begin(); featureIt != percolatorFeatures.end(); ++featureIt, ++i) {
    featureNames.insertFeature(featureIt->first);
    if (featureIt->second != 0.0) {
      hasDefaultValues = true;
      init_values[i] = featureIt->second;
      if (Globals::VERB > 3) {
        std::cerr << "Initial direction for " << featureIt->first << " is " << 
                     featureIt->second << std::endl;
      }
    }
  }
  bool calcDoc = false;
  featureNames.initFeatures(calcDoc);
  
  setHandler_.getFeaturePool().createPool(FeatureNames::getNumFeatures());
  
  DataSet* targetSet = new DataSet();
  assert(targetSet);
  targetSet->setLabel(1);
  DataSet* decoySet = new DataSet();
  assert(decoySet);
  decoySet->setLabel(-1);    

  setHandler_.push_back_dataset(targetSet);
  setHandler_.push_back_dataset(decoySet);
  
  pCheck_ = SanityCheck::initialize("");
  assert(pCheck_);
  pCheck_->checkAndSetDefaultDir();
  if (hasDefaultValues) pCheck_->addDefaultWeights(init_values);
  bool concatenatedSearch = true;
  pCheck_->setConcatenatedSearch(concatenatedSearch);
  
  if (maxPSMs_ > 0u) {
    /* the column layout of SetHandler::writeTab, which readAndScoreTab reads */
    allPsmsFN_ = psmResultFN_ + ".all.pin";
    allPsmsStream_.open(allPsmsFN_.c_str(), ios::out);
    if (!allPsmsStream_) {
      throw MyException("ERROR: cannot write " + allPsmsFN_ + "\n");
    }
    allPsmsStream_ << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\t" 
                   << featureNames.getFeatureNames(true) << "\tPeptide\tProteins\n";
    allPsmsStream_.precision(10);
  }
}

void PercolatorAdapter::addPsm(const std::string& specId, int label, int scannr, 
    const std::string& peptide, const std::vector<double>& features) {
  size_t priority = 0u;
  if (maxPSMs_ > 0u) {
    writePsm(specId, label, scannr, peptide, features);
    ++numPsmsAdded_;
    priority = scanPriority(scannr);
    if (subsetPsms_.size() >= maxPSMs_ && priority >= subsetUpperLimit_) return;
  }
  
  PSMDescription* myPsm = new PSMDescription();
  
  // adding n-term and c-term residues to peptide
  myPsm->peptide = peptide;

  myPsm->setId(specId);
  myPsm->scan = scannr;
  myPsm->expMass = 0.0;
  myPsm->calcMass = 0.0;

  myPsm->features = setHandler_.getFeaturePool().allocate();
  
  std::vector<double>::const_iterator ftIt;
  int i = 0;
  for (ftIt = features.begin(); ftIt != features.end(); ++ftIt, ++i) {
    myPsm->features[i] = *ftIt;
  }
  
  if (maxPSMs_ == 0u) {
    setHandler_.getSubsetFromLabel(label)->registerPsm(myPsm);
    return;
  }
  
  /* reservoir sampling as in SetHandler::readPSMs */
  SubsetPsm psmPriority;
  psmPriority.psm = myPsm;
  psmPriority.priority = priority;
  psmPriority.order = numPsmsAdded_;
  psmPriority.label = label;
  subsetPsms_.push(psmPriority);
  if (subsetPsms_.size() > maxPSMs_) {
    SubsetPsm del = subsetPsms_.top();
    subsetUpperLimit_ = del.priority;
    setHandler_.getFeaturePool().deallocate(del.psm->features);
    PSMDescription::deletePtr(del.psm);
    subsetPsms_.pop();
  }
}

void PercolatorAdapter::writePsm(const std::string& specId, int label, int scannr, 
    const std::string& peptide, const std::vector<double>& features) {
  allPsmsStream_ << specId << '\t' << label << '\t' << scannr << "\t0\t0";
  std::vector<double>::const_iterator ftIt;
  for (ftIt = features.begin(); ftIt != features.end(); ++ftIt) {
    allPsmsStream_ << '\t' << *ftIt;
  }
  /* the tab reader fails on a line without a protein column */
  allPsmsStream_ << '\t' << peptide << "\tfeature\n";
}

/* As Caller::run does after training on a subset: score every PSM with the
   weights averaged over the cross-validation folds, reading the PSMs back
   one at a time so that their features are not all held in memory. */
void PercolatorAdapter::scoreAllPsms(CrossValidation& crossValidation, 
    Scores& allScores) {
  if (Globals::VERB > 0) {
    std::cerr << "Scoring full list of PSMs with trained SVMs." << std::endl;
  }
  std::vector<double> rawWeights;
  crossValidation.getAvgWeights(rawWeights, pNorm_);
  setHandler_.reset();
  allScores.reset();
  
  std::ifstream allPsms(allPsmsFN_.c_str(), ios::in);
  if (!setHandler_.readAndScoreTab(allPsms, rawWeights, allScores, pCheck_)) {
    throw MyException("ERROR: cannot read back " + allPsmsFN_ + "\n");
  }
  allPsms.close();
  remove(allPsmsFN_.c_str());
  
  /* before calculatePSMProb's TDC drops the losing PSM of each scan */
  scoredPsms_.reserve(allScores.size());
  std::vector<ScoreHolder>::iterator scoreIt;
  for (scoreIt = allScores.begin(); scoreIt != allScores.end(); ++scoreIt) {
    scoredPsms_.push_back(scoreIt->pPSM);
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Evaluated set contained " << allScores.posSize()
              << " positives and " << allScores.negSize() << " negatives." << std::endl;
  }
  allScores.postMergeStep();
  allScores.calcQ(selectionFdr_);
  allScores.normalizeScores(selectionFdr_);
}

} /* namespace quandenser */
