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

#ifndef QUANDENSER_PERCOLATORADAPTER_H_
#define QUANDENSER_PERCOLATORADAPTER_H_

#include <vector>
#include <string>
#include <fstream>

#include "percolator/src/Caller.h"

#include "DinosaurFeature.h"

namespace quandenser {

/* Runs Percolator on the candidate links of a feature-matching step. The
   PSMs are written to a tab-delimited file next to the result file and
   Percolator reads that file with Caller::run, as the percolator executable
   does, so that its options for large inputs apply: --subset-max-train
   (-N) samples spectra while reading, and the SVMs train one fold at a time
   by default. */
class PercolatorAdapter : public Caller {
 public:
  PercolatorAdapter() : Caller() {}
  
  bool parseOptions(const std::vector<std::string>& percolatorArgs) {
    std::vector<const char*> percolatorArgv;
    for (std::vector<std::string>::const_iterator it = percolatorArgs.begin();
       it != percolatorArgs.end(); ++it) {
      percolatorArgv.push_back(it->c_str());
    }
    return Caller::parseOptions(percolatorArgs.size(), (char**)&percolatorArgv.front());
  }
  
  void init(const std::vector<std::pair<std::string, double> >& percolatorFeatures);
  void addPsm(const std::string& specId, int label, int scannr, 
    const std::string& peptide, const std::vector<double>& features);
  
  void process();
 protected:
  std::string pinFN_;
  std::ofstream pinStream_;
  std::string line_;
};

} /* namespace quandenser */

#endif /* QUANDENSER_PERCOLATORADAPTER_H_ */
