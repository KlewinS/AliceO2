// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file readCLusterMCtruth.cxx
/// \brief This macro demonstrates how to extract the MC truth information from
/// the clusters
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TFile.h"
#include "TTree.h"

#include "TPCReconstruction/Cluster.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "DataFormatsTPC/Helpers.h"
#include "DataFormatsTPC/ClusterHardware.h"

#include <vector>
#else
R__LOAD_LIBRARY(libTPCReconstruction)
R__LOAD_LIBRARY(libDataFormatsTPC)
#endif
void readCluster(std::string filename)
{
  TFile *clusterFile = TFile::Open(filename.data());
  TTree *clusterTree = (TTree*)clusterFile->Get("o2sim");

  std::vector<o2::DataFormat::TPC::ClusterHardwareContainer8kb>* clusterContainer = nullptr;
  clusterTree->SetBranchAddress("TPCClusterHW",&clusterContainer);

  for(int iEvent=0; iEvent<clusterTree->GetEntriesFast(); ++iEvent) {
    int cluster = 0;
    clusterTree->GetEntry(iEvent);
    for (int cru = 0; cru < 360; ++ cru){
      for(auto& cont : *clusterContainer) {
        if (cont.getContainer()->CRU != cru) continue;
        if (cont.getContainer()->timeBinOffset == 0xFFFFFFFF) continue;
        std::cout << "CRU [" << cont.getContainer()->CRU << "] "
                  << "has [" << cont.getContainer()->numberOfClusters << "] clusters "
                  << "with time bin offset [" << cont.getContainer()->timeBinOffset << "] "
                  << std::endl;
      }
    }
  }
}
