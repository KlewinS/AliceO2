// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <fstream>
#include <iomanip>
#include <boost/tokenizer.hpp>
#include <vector>
#include <map>
#include <memory>

#include "TPCBase/Mapper.h"
#include "TPCBase/Defs.h"
#include "TPCBase/PadPos.h"
#include "TPCBase/CRU.h"

/*
gSystem->AddIncludePath("-I$O2_ROOT/include -I$FAIRROOT_ROOT/include");
.L checkCRUSorterMapping.C+
*/

using namespace o2::TPC;

void checkCRUSorterMapping(std::string basePath, int partition)
{
  const Mapper& mapper = Mapper::instance();
  const PartitionInfo& partitionInfo = mapper.getPartitionInfo(partition);
  int FECmax = partitionInfo.getNumberOfFECs();

  std::cout << "####################################" << std::endl;
  std::cout << "## Checking partition " << partition << "           ##" << std::endl;
  std::cout << "##                                ##" << std::endl;
  std::cout << "## Opening files for FECs 0 to " << FECmax-1 << " ##" << std::endl;
  std::cout << "## of region " << partition*2 << " and " << (partition*2+1) << "              ##" << std::endl;
  std::cout << "####################################" << std::endl << std::endl;

  std::vector<std::ifstream> inFiles_l;
  std::vector<std::ifstream> inFiles_h;
  for (int fec = 0; fec < FECmax; ++fec) {
    std::stringstream fileToCheck;
    fileToCheck << basePath;
    fileToCheck << "/mapping_region" << (partition*2);
    fileToCheck << "_fec" << std::setw(2) << std::setfill('0') << fec;
    fileToCheck << ".txt";
    inFiles_l.emplace_back(fileToCheck.str());
    if (!inFiles_l.back().is_open()) {
      std::cout << "Can't open file " << fileToCheck.str() << std::endl;
      return;
    }
  }

  for (int fec = 0; fec < FECmax; ++fec) {
    std::stringstream fileToCheck;
    fileToCheck << basePath;
    fileToCheck << "/mapping_region" << ((partition*2)+1);
    fileToCheck << "_fec" << std::setw(2) << std::setfill('0') << fec;
    fileToCheck << ".txt";
    inFiles_h.emplace_back(fileToCheck.str());
    if (!inFiles_h.back().is_open()) {
      std::cout << "Can't open file " << fileToCheck.str() << std::endl;
      return;
    }
  }


  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep{" "};
  std::vector<std::unique_ptr<std::vector<int>>> allChannels;
  std::string line;
  for (auto &infile : inFiles_l) {
    unsigned row = 0;
    while (std::getline(infile,line)) {
//      std::cout << "Row " << row << ": " << line << std::endl;
      if (allChannels.size() < row+1)
        allChannels.emplace_back(std::make_unique<std::vector<int>>());

      tokenizer tok{line,sep};
      for (auto &t : tok) allChannels[row]->emplace_back(std::stoi(t));
      ++row;
    }
  }
  for (auto &infile : inFiles_h) {
    unsigned row = mapper.getNumberOfRowsRegion(partition*2);
    while (std::getline(infile,line)) {
//      std::cout << "Row " << row << ": " << line << std::endl;
      if (allChannels.size() < row+1)
        allChannels.emplace_back(std::make_unique<std::vector<int>>());

      tokenizer tok{line,sep};
      for (auto &t : tok) allChannels[row]->emplace_back(std::stoi(t));
      ++row;
    }
  }

//  unsigned rowNumber = 0;
//  for (auto &row : allChannels) {
//    std::cout << "Row " << rowNumber++ << ": ";
//    for (auto &pad : *row) {
//      std::cout << pad << " ";
//    }
//    std::cout << std::endl;
//  }

  unsigned errorsFound = 0;
  unsigned padOffset = 0;
  for (short p=0; p < partition; ++p) {
    const PartitionInfo& partInfo = mapper.getPartitionInfo(p);
    padOffset += partInfo.getNumberOfPads();
  }

  for (int row = 0; row < mapper.getNumberOfRowsPartition(partition); ++row) {
//    std::cout << "Row " << row << ": " << std::endl;
    for (int pad = 0; pad < mapper.getNumberOfPadsInRowPartition(partition,row); ++pad) {
      const FECInfo& fecInfo = mapper.fecInfo(padOffset+mapper.getPadNumberInPartition(partition,row,pad));
      const int fecChannel = fecInfo.getSampaChip()*32 + fecInfo.getSampaChannel();
//      std::cout << mapper.getPadNumberInPartition(partition,row,pad) << " ";
//      std::cout << fecChannel << " ";

      if (fecChannel != allChannels[row]->at(pad)) {
          ++errorsFound;
          std::cout << "ERROR: missmatch for pad [" << pad << "] in row [" << row << "] (FEC " << static_cast<unsigned>(fecInfo.getIndex()) << "), "
                    << "FEC channel should be [" << fecChannel << "] but is [" << allChannels[row]->at(pad) << "]."
                    << std::endl;
      }
    }
//    std::cout << std::endl;
  }

  if (errorsFound == 0) {
    std::cout << "Mapping files are good, no errors found." << std::endl;
  }
  std::cout << std::endl;

}

