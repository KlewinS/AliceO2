// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/tokenizer.hpp>
#include <vector>
#include <array>
#include <algorithm>
#include <map>
#include <memory>

#include "TPCBase/Mapper.h"
#else
R__LOAD_LIBRARY(libTPCBase)
#endif
//#include "TPCBase/Defs.h"
//#include "TPCBase/PadPos.h"
//#include "TPCBase/CRU.h"

/*
gSystem->AddIncludePath("-I$O2_ROOT/include -I$FAIRROOT_ROOT/include");
.L checkCRUSorterMapping.C+
*/

using namespace o2::TPC;

bool checkRegionEmpty(std::vector<std::vector<std::array<int, 24 * 6>>>& toCheck, int row, int p) {
  for (int r = 0; r < 10; ++r) {
    if (row >= toCheck[r].size()) continue;
    for (int i = 0; i < 6; ++i){
//      std::cout << "Check region " << r << " row "<< row << " " << p-i << " " <<  toCheck[r][row][p-i] << std::endl;
      if (toCheck[r][row][p-i] != -1) return false;
    }
  }
  return true;
}

bool checkPadsEmpty(std::vector<std::vector<std::array<int, 24 * 6>>>& toCheck, int region, int row, int p) {
  for (int i = 0; i < 6; ++i){
//    std::cout << "Check pads in region "<< region << " " << toCheck[region][row][p-i] << std::endl;
    if (toCheck[region][row][p-i] != -1) return false;
  }
  return true;
}


void printCRUMergerMapping()
{
  const Mapper& mapper = Mapper::instance();

  int padOffset = 0;
  int fecOffset = 0;
  std::vector<std::vector<std::array<int, 24 * 6>>> fecIDs;
  std::vector<std::vector<std::array<int, 24 * 6>>> padsOfFecRow;
  fecIDs.resize(10);
  padsOfFecRow.resize(10);
  for (int region = 0; region < 10; ++region) {
    if (region > 0) padOffset += mapper.getMapPadRegionInfo()[region-1].getNumberOfPads();
    if ((region > 1) && region%2 == 0) fecOffset += mapper.getMapPartitionInfo()[(region/2)-1].getNumberOfFECs();

    const int rows = mapper.getPadRegionInfo(region).getNumberOfPadRows();
    fecIDs[region].resize(rows);
    padsOfFecRow[region].resize(rows);
    for (int row = 0; row < rows; ++row) {
      std::fill(fecIDs[region][row].begin(),fecIDs[region][row].end(),-1);
      std::fill(padsOfFecRow[region][row].begin(),padsOfFecRow[region][row].end(),-1);

      int lastFecIndex = -1;
      int padCount = 0;
      const int pads = mapper.getPadRegionInfo(region).getPadsInRowRegion(row);
      for (int pad = 0; pad < pads; ++pad) {
        const FECInfo& fecInfo = mapper.fecInfo(padOffset+mapper.getPadNumberInRegion(region,row,pad));

        fecIDs[region][row][pad+2] = fecInfo.getIndex() - fecOffset;
        if (lastFecIndex != static_cast<int>(fecInfo.getIndex())) padCount = 0;
        padsOfFecRow[region][row][pad+2] = padCount++;
        lastFecIndex = static_cast<int>(fecInfo.getIndex());

      }
    }
  }

  for (int row = 0; row < 18; ++row) {
    std::cout << "  gen_row_" << row << " : if G_ROW_NUMBER = " << row << " generate" << std::endl;
    std::cout << "    fsm_row : process(current_s,full_row,config.region)" << std::endl;
    std::cout << "    begin" << std::endl;
    std::cout << "      case current_s is" << std::endl;
    for (int p = (24 * 6) - 1; p >= 0; p = p - 6) {
      if (checkRegionEmpty(fecIDs,row,p)) continue;
      std::cout << "        when S_segment_" << p/6 << " =>" << std::endl;
      std::cout << "          case config.region is" << std::endl;
      for (int region = 0; region < 10; ++region) {
        if (row >= fecIDs[region].size()) continue;
        if (checkPadsEmpty(fecIDs,region,row,p)) continue;

        std::cout << "            when " << region << " => row_segment <= (";
        for (int i = 0; i < 6; ++i) {
          if (i != 0) std::cout << ", ";
          if (fecIDs[region][row][p-i] == -1) {
            std::cout << 5-i << " => (others => '0')";
          } else {
            std::cout << 5-i << " => full_row(" << fecIDs[region][row][p-i] << "*7+" << padsOfFecRow[region][row][p-i] << ")";
          }
        }
        std::cout << ");" << std::endl;
      }

      std::cout << "            when others => row_segment <= (others => (others => '0'));" << std::endl;
      std::cout << "          end case;" << std::endl;
    }
    std::cout << "        when others => row_segment <= (others => (others => '0'));" << std::endl;
    std::cout << "      end case;" << std::endl;
    std::cout << "    end process fsm_row;" << std::endl;
    std::cout << "  end generate gen_row_" << row << ";" << std::endl << std::endl;

  }

//  for (int p = (24 * 6) - 1; p >= 0; p = p - 6) {
//    std::cout << "       when S_segment_" << p/6 << " =>" << std::endl;
//    std::cout << "         case G_ROW_NUMBER is" << std::endl;
//    for (int row = 0; row < 18; ++row) {
//      if (checkRegionEmpty(fecIDs,row,p)) continue;
//
//      std::cout << "           when " << row << " =>" << std::endl;
//      std::cout << "             case config.region is" << std::endl;
//      for (int region = 0; region < 10; ++region) {
//        if (row >= fecIDs[region].size()) continue;
//        if (checkPadsEmpty(fecIDs,region,row,p)) continue;
//
//        std::cout << "               when " << region << " => row_segment <= (";
//        for (int i = 0; i < 6; ++i) {
//          if (i != 0) std::cout << ", ";
//          if (fecIDs[region][row][p-i] == -1) {
//            std::cout << 5-i << " => (others => '0')";
//          } else {
//            std::cout << 5-i << " => full_row(" << fecIDs[region][row][p-i] << "*7+" << padsOfFecRow[region][row][p-i] << ")";
//          }
//        }
//        std::cout << ");" << std::endl;
////        std::cout << fecIDs[region][row][p-0] << " "
////                  << fecIDs[region][row][p-1] << " "
////                  << fecIDs[region][row][p-2] << " "
////                  << fecIDs[region][row][p-3] << " "
////                  << fecIDs[region][row][p-4] << " "
////                  << fecIDs[region][row][p-5] << std::endl;
////
////        std::cout << padsOfFecRow[region][row][p-0] << " "
////                  << padsOfFecRow[region][row][p-1] << " "
////                  << padsOfFecRow[region][row][p-2] << " "
////                  << padsOfFecRow[region][row][p-3] << " "
////                  << padsOfFecRow[region][row][p-4] << " "
////                  << padsOfFecRow[region][row][p-5] << std::endl << std::endl;
//      }
//      std::cout << "               when others => row_segment <= (others => (others => '0'));" << std::endl;
//      std::cout << "             end case;" << std::endl << std::endl;
//    }
//    std::cout << "           when others => row_segment <= (others => (others => '0'));" << std::endl;
//    std::cout << "         end case;" << std::endl << std::endl;
//  }
  return;

}

