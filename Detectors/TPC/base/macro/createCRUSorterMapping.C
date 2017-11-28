// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <vector>
#include <sstream>
#include <memory>
#include <iomanip>
#include <fstream>

#include "TPCBase/Mapper.h"
#include "TPCBase/PadPos.h"
#include "TPCBase/CRU.h"

/*
gSystem->AddIncludePath("-I$O2_ROOT/include -I$FAIRROOT_ROOT/include");
.L createCRUSorterMapping.C+
*/

using namespace o2::TPC;

int getRamAddress(int sampa, int channel) {
  if (sampa == 0 || sampa == 3) return (char) channel;
  if (sampa == 1 || sampa == 4) return (char) (32+channel);
  if (sampa == 2) return (char) (64+(channel%16));

  return 0;
}

void createCRUSorterMapping()
{
  const Mapper& mapper = Mapper::instance();

//  for (short fec = 0; fec < 91; ++fec) {
  for (int partition = 0; partition < 5; ++partition) {
    std::cout << std::endl;
    std::cout << "#################" << std::endl;
    std::cout << "## Partition " << partition << " ##" << std::endl;
    std::cout << "#################" << std::endl;

    const PartitionInfo& partitionInfo = mapper.getPartitionInfo(partition);
    for (int fec = 0; fec < partitionInfo.getNumberOfFECs(); ++fec) {
      std::map<int,std::map<int,FECInfo>> rows;
      for (int sampa = 0; sampa < 5; ++sampa) {
        for (int channel = 0; channel < 32; ++channel) {
          const PadPos& padPos = mapper.padPos(partition,fec,sampa,channel);
          const FECInfo fecInfo(fec,sampa,channel);
          rows[(int)padPos.getRow()-partitionInfo.getGlobalRowOffset()][(int)padPos.getPad()] = fecInfo;
        }
      }

      short wordCount = 0;
      ofstream outfile;
      unsigned count = 0;
      short marker;
      std::array<std::unique_ptr<std::vector<std::pair<int,int>>>,2> words = {std::make_unique<std::vector<std::pair<int,int>>>(), std::make_unique<std::vector<std::pair<int,int>>>() };
      short region = 0;
      for(auto &r : rows) {
        for(auto &pad : r.second) {
//          std::cout << "Row [" << r.first << "] "
//                    << "Pad [" << pad.first << "] "
//                    << "SAMPA [" << (int)pad.second.getSampaChip() << "] "
//                    << "Channel [" << (int)pad.second.getSampaChannel() << "] "
//                    << "RAM [" << (int)getRamAddress((int)pad.second.getSampaChip(),(int)pad.second.getSampaChannel()) << "] "
//                    << std::endl;

          words[wordCount]->emplace_back(std::make_pair(r.first,getRamAddress((int)pad.second.getSampaChip(),(int)pad.second.getSampaChannel())));
          wordCount = (wordCount+1)%2;

          ++count;
          if (count == 80) {

            outfile.close();
            std::stringstream filename;
            filename << "sorterMapping_region" << partition*2+region << "_fec";
            filename << std::setw(2) << std::setfill('0') << fec;
            filename << ".txt";
            std::cout << "writing file " << filename.str() << std::endl;
            outfile.open(filename.str());

            if (region == 0) region = 1;
            else region = 0;

//            std::cout << words[0]->size() << " " << words[1]->size() << std::endl << std::endl;
            for (int i = 0; i < words[0]->size(); ++i) {
              marker = 0x0;
              if (i == words[0]->size()-1)                              marker = 0x3;
              else if (words[0]->at(i).first < words[1]->at(i).first)   marker = 0x1;
              else if (i > 0)  {
                if (words[1]->at(i).first != words[0]->at(i+1).first)   marker = 0x2;
              }
//              std::cout << words[0]->at(i).second << " " << words[1]->at(i).second << "\t" << words[0]->at(i).first << " " << words[1]->at(i).first << "\t" << marker << std::endl;
              outfile << std::setfill('0') << std::setw(4) << std::hex << ((0x7F & words[0]->at(i).second) | ((0x7F & words[1]->at(i).second) << 7) | ((0x3 & marker) << 14)) << std::dec << std::endl;
            }
            words[0]->clear(); words[1]->clear();
            count = 0;
          }
        }
      }
      outfile.close();
    }
  }
}
