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

#include <charconv>
#include <cstdio>

#include "PercolatorAdapter.h"

namespace quandenser {

namespace {
/* the shortest text that reads back as the same double */
void appendDouble(std::string& line, double value) {
  char buffer[32];
  std::to_chars_result res = std::to_chars(buffer, buffer + sizeof(buffer), value);
  line.append(buffer, res.ptr);
}
}

void PercolatorAdapter::init(const std::vector<std::pair<std::string, double> >& percolatorFeatures) {
  pinFN_ = psmResultFN_ + ".pin";
  pinStream_.open(pinFN_.c_str(), std::ios::out);
  if (!pinStream_) {
    throw MyException("ERROR: cannot write " + pinFN_ + "\n");
  }
  
  std::string defaultDirection = "DefaultDirection\t-\t-";
  bool hasDefaultValues = false;
  pinStream_ << "SpecId\tLabel\tScanNr";
  std::vector<std::pair<std::string, double> >::const_iterator featureIt;
  for (featureIt = percolatorFeatures.begin(); featureIt != percolatorFeatures.end(); ++featureIt) {
    pinStream_ << '\t' << featureIt->first;
    defaultDirection += '\t';
    appendDouble(defaultDirection, featureIt->second);
    if (featureIt->second != 0.0) hasDefaultValues = true;
  }
  pinStream_ << "\tPeptide\tProteins\n";
  if (hasDefaultValues) pinStream_ << defaultDirection << '\n';
}

void PercolatorAdapter::addPsm(const std::string& specId, int label, int scannr, 
    const std::string& peptide, const std::vector<double>& features) {
  line_ = specId;
  line_ += (label == 1) ? "\t1\t" : "\t-1\t";
  line_ += std::to_string(scannr);
  std::vector<double>::const_iterator ftIt;
  for (ftIt = features.begin(); ftIt != features.end(); ++ftIt) {
    line_ += '\t';
    appendDouble(line_, *ftIt);
  }
  line_ += '\t';
  line_ += peptide;
  line_ += "\tfeature\n"; /* Percolator expects at least one protein */
  pinStream_ << line_;
}

void PercolatorAdapter::process() {
  pinStream_.close();
  if (pinStream_.fail()) {
    throw MyException("ERROR: cannot write " + pinFN_ + "\n");
  }
  inputFNs_.assign(1u, pinFN_);
  if (!run()) {
    throw MyException("ERROR: Percolator failed on " + pinFN_ + "\n");
  }
  std::remove(pinFN_.c_str());
}

} /* namespace quandenser */
