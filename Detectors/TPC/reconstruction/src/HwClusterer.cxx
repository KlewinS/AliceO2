// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliTPCUpgradeHwClusterer.cxx
/// \brief Hwclusterer for the TPC

#include "TPCReconstruction/HwCluster.h"
#include "TPCReconstruction/HwClusterer.h"
#include "TPCReconstruction/HwClusterFinder.h"
#include "TPCReconstruction/ClusterContainer.h"
#include "TPCBase/Digit.h"
#include "TPCBase/PadPos.h"
#include "TPCBase/CRU.h"
#include "TPCBase/PadSecPos.h"
#include "TPCBase/CalArray.h"

#include "FairLogger.h"
#include "TMath.h"
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>

std::mutex g_display_mutex;

using namespace o2::TPC;

//________________________________________________________________________
HwClusterer::HwClusterer(std::vector<o2::TPC::HwCluster> *output, Processing processingType, int globalTime, int cruMin, int cruMax, float minQDiff,
    bool assignChargeUnique, bool enableNoiseSim, bool enablePedestalSubtraction, int padsPerCF, int timebinsPerCF, int cfPerRow)
  : Clusterer()
  , mClusterArray(output)
  , mProcessingType(processingType)
  , mGlobalTime(globalTime)
  , mCRUMin(cruMin)
  , mCRUMax(cruMax)
  , mMinQDiff(minQDiff)
  , mAssignChargeUnique(assignChargeUnique)
  , mEnableNoiseSim(enableNoiseSim)
  , mEnablePedestalSubtraction(enablePedestalSubtraction)
  , mIsContinuousReadout(true)
  , mPadsPerCF(padsPerCF)
  , mTimebinsPerCF(timebinsPerCF)
  , mCfPerRow(cfPerRow)
  , mLastTimebin(-1)
  , mNoiseObject(nullptr)
  , mPedestalObject(nullptr)
{
}

//________________________________________________________________________
HwClusterer::~HwClusterer()
{
  LOG(DEBUG) << "Enter Destructor of HwClusterer" << FairLogger::endl;

  for (auto it : mClusterFinder) {
    for (auto itt : it) {
      for (auto ittt : itt) {
        delete ittt;
      }
    }
  }
}

//________________________________________________________________________
void HwClusterer::Init()
{
  LOG(DEBUG) << "Enter Initializer of HwClusterer" << FairLogger::endl;

  /*
   * initalize all cluster finder
   */
  mCfPerRow = (int)ceil((double)(mPadsMax+2+2)/(mPadsPerCF-2-2));
  mClusterFinder.resize(mCRUMax);
  const Mapper& mapper = Mapper::instance();
  for (int iCRU = mCRUMin; iCRU < mCRUMax; iCRU++){
    mClusterFinder[iCRU].resize(mapper.getNumberOfRowsPartition(iCRU));
    for (int iRow = 0; iRow < mapper.getNumberOfRowsPartition(iCRU); iRow++){
      mClusterFinder[iCRU][iRow].resize(mCfPerRow);
      for (int iCF = 0; iCF < mCfPerRow; iCF++){
        int padOffset = iCF*(mPadsPerCF-2-2)-2;
        mClusterFinder[iCRU][iRow][iCF] = new HwClusterFinder(iCRU,iRow,iCF,padOffset,mPadsPerCF,mTimebinsPerCF,mMinQDiff,mMinQMax,mRequirePositiveCharge);
        mClusterFinder[iCRU][iRow][iCF]->setAssignChargeUnique(mAssignChargeUnique);


        /*
         * Connect always two CFs to be able to comunicate found clusters. So
         * the "right" one can tell the one "on the left" which pads were
         * already used for a cluster.
         */
        if (iCF != 0) {
          mClusterFinder[iCRU][iRow][iCF]->setNextCF(mClusterFinder[iCRU][iRow][iCF-1]);
        }
      }
    }
  }


  /*
   * vector of HwCluster vectors, one vector for each CRU (possible thread)
   * to store the clusters found there
   */
  mClusterStorage.resize(mCRUMax);


  /*
   * vector of digit vectors, one vector for each CRU (possible thread) to
   * store there only those digits which are relevant for this particular
   * CRU (thread)
   */
  mDigitContainer.resize(mCRUMax);
  for (int iCRU = mCRUMin; iCRU < mCRUMax; iCRU++)
    mDigitContainer[iCRU].resize(mapper.getNumberOfRowsPartition(iCRU));

}

//________________________________________________________________________
void HwClusterer::processDigits(
    const std::vector<std::vector<std::pair<Digit const*, gsl::span<const o2::MCCompLabel>>>>& digits,
    const std::vector<std::vector<HwClusterFinder*>>& clusterFinder,
          std::vector<HwCluster>& cluster,
    const CfConfig config)
{
  int timeDiff = (config.iMaxTimeBin+1) - config.iMinTimeBin;
  if (timeDiff < 0) return;
//  std::thread::id this_id = std::this_thread::get_id();
//  g_display_mutex.lock();
//  std::cout << "thread " << this_id << " started.\n";
//  g_display_mutex.unlock();

  float iAllBins[timeDiff][config.iMaxPads];

  for (int iRow = 0; iRow < config.iMaxRows; iRow++){

    /*
     * prepare local storage
     */
    short t,p;
    float noise;
    if (config.iEnableNoiseSim && config.iNoiseObject != nullptr) {
      for (t=timeDiff; t--;) {
        for (p=config.iMaxPads; p--;) {
          iAllBins[t][p] = config.iNoiseObject->getValue(CRU(config.iCRU),iRow,p);
        }
      }
    } else {
      std::fill(&iAllBins[0][0], &iAllBins[0][0]+timeDiff*config.iMaxPads, 0.f);
    }

    /*
     * fill in digits
     */
    for (auto& digit : digits[iRow]){
      const Int_t iTime         = digit.first->getTimeStamp();
      const Int_t iPad          = digit.first->getPad() + 2;  // offset to have 2 empty pads on the "left side"
      const Float_t charge      = digit.first->getChargeFloat();

//      std::cout << iCRU << " " << iRow << " " << iPad << " " << iTime << " (" << iTime-minTime << "," << timeDiff << ") " << charge << std::endl;
      iAllBins[iTime-config.iMinTimeBin][iPad] += charge;
      if (config.iEnablePedestalSubtraction && config.iPedestalObject != nullptr) {
        const float pedestal = config.iPedestalObject->getValue(CRU(config.iCRU),iRow,iPad-2);
        //printf("digit: %.2f, pedestal: %.2f\n", iAllBins[iTime-config.iMinTimeBin][iPad], pedestal);
        iAllBins[iTime-config.iMinTimeBin][iPad] -= pedestal;
      }
    }

    /*
     * copy data to cluster finders
     */
    const Short_t iPadsPerCF = clusterFinder[iRow][0]->getNpads();
    const Short_t iTimebinsPerCF = clusterFinder[iRow][0]->getNtimebins();
    std::vector<std::vector<HwClusterFinder*>::const_reverse_iterator> cfWithCluster;
    unsigned time,pad;
    for (time = 0; time < timeDiff; ++time){    // ordering important!!
      for (pad = 0; pad < config.iMaxPads; pad = pad + (iPadsPerCF -2 -2 )) {
        const Short_t cf = pad / (iPadsPerCF-2-2);
        clusterFinder[iRow][cf]->AddTimebin(&iAllBins[time][pad],time+config.iMinTimeBin,(config.iMaxPads-pad)>=iPadsPerCF?iPadsPerCF:(config.iMaxPads-pad));
      }

      /*
       * search for clusters and store reference to CF if one was found
       */
      if (clusterFinder[iRow][0]->getTimebinsAfterLastProcessing() == iTimebinsPerCF-2 -2)  {
        /*
         * ordering is important: from right to left, so that the CFs could inform each other if cluster was found
         */
        for (auto rit = clusterFinder[iRow].crbegin(); rit != clusterFinder[iRow].crend(); ++rit) {
          if ((*rit)->findCluster()) {
            cfWithCluster.push_back(rit);
          }
        }
      }
    }

    /*
     * add empty timebins to find last clusters
     */
    if (!config.iIsContinuousReadout) {
      // +2 so that for sure all data is processed
      for (time = 0; time < clusterFinder[iRow][0]->getNtimebins()+2; ++time){
        for (auto rit = clusterFinder[iRow].crbegin(); rit != clusterFinder[iRow].crend(); ++rit) {
          (*rit)->AddZeroTimebin(time+timeDiff+config.iMinTimeBin,iPadsPerCF);
        }

        /*
         * search for clusters and store reference to CF if one was found
         */
        if (clusterFinder[iRow][0]->getTimebinsAfterLastProcessing() == iTimebinsPerCF-2 -2)  {
          /*
           * ordering is important: from right to left, so that the CFs could inform each other if cluster was found
           */
          for (auto rit = clusterFinder[iRow].crbegin(); rit != clusterFinder[iRow].crend(); ++rit) {
            if ((*rit)->findCluster()) {
              cfWithCluster.push_back(rit);
            }
          }
        }
      }
      for (auto rit = clusterFinder[iRow].crbegin(); rit != clusterFinder[iRow].crend(); ++rit) {
        (*rit)->setTimebinsAfterLastProcessing(0);
      }
    }

    /*
     * collect found cluster
     */
    for (std::vector<HwClusterFinder*>::const_reverse_iterator &cf_it : cfWithCluster) {
      std::vector<HwCluster>* cc = (*cf_it)->getClusterContainer();
      for (HwCluster& c : *cc){
        cluster.push_back(c);
      }
      (*cf_it)->clearClusterContainer();
    }

  }

//  g_display_mutex.lock();
//  std::cout << "thread " << this_id << " finished.\n";
//  g_display_mutex.unlock();
}

//________________________________________________________________________
void HwClusterer::Process(std::vector<o2::TPC::Digit> const &digits, MCLabel const* mcDigitTruth, MCLabel &mcClusterTruth
    )
{
  mClusterArray->clear();

  /*
   * clear old storages
   */
  for (auto& cs : mClusterStorage) cs.clear();
  for (auto& dc : mDigitContainer ) {
    for (auto& dcc : dc) dcc.clear();
  }

  int iTimeBin;
  int iTimeBinMin = (mIsContinuousReadout)?mLastTimebin + 1 : 0;
  //int iTimeBinMin = mLastTimebin + 1;
  int iTimeBinMax = mLastTimebin;

  /*
   * Loop over digits
   */
  int digitIndex = 0;
  for (const auto& digit : digits) {

    /*
     * add current digit to storage
     */

    iTimeBin = digit.getTimeStamp();
    if (iTimeBin < iTimeBinMin) continue;
    iTimeBinMax = std::max(iTimeBinMax,iTimeBin);
    if (mcDigitTruth == nullptr)
      mDigitContainer[digit.getCRU()][digit.getRow()].push_back(std::make_pair(&digit,nullptr));
    else
      mDigitContainer[digit.getCRU()][digit.getRow()].push_back(std::make_pair(&digit,mcDigitTruth->getLabels(digitIndex)));
    ++digitIndex;
  }
  ProcessTimeBins(iTimeBinMin, iTimeBinMax);
}

//________________________________________________________________________
void HwClusterer::Process(std::vector<std::unique_ptr<Digit>>& digits, MCLabel const* mcDigitTruth, MCLabel &mcClusterTruth)
{
  mClusterArray->clear();

  /*
   * clear old storages
   */
  for (auto& cs : mClusterStorage) cs.clear();
  for (auto& dc : mDigitContainer ) {
    for (auto& dcc : dc) dcc.clear();
  }

  int iTimeBin;
  int iTimeBinMin = (mIsContinuousReadout)?mLastTimebin + 1 : 0;
  int iTimeBinMax = mLastTimebin;

  /*
   * Loop over digits
   */
  int digitIndex = 0;
  for (auto& digit_ptr : digits) {
    Digit* digit = digit_ptr.get();

    /*
     * add current digit to storage
     */
    iTimeBin = digit->getTimeStamp();
    if (iTimeBin < iTimeBinMin) continue;
    iTimeBinMax = std::max(iTimeBinMax,iTimeBin);
    if (mcDigitTruth == nullptr)
      mDigitContainer[digit->getCRU()][digit->getRow()].push_back(std::make_pair(digit,nullptr));
    else
      mDigitContainer[digit->getCRU()][digit->getRow()].push_back(std::make_pair(digit,mcDigitTruth->getLabels(digitIndex)));
    ++digitIndex;
  }

  ProcessTimeBins(iTimeBinMin, iTimeBinMax);
}

void HwClusterer::ProcessTimeBins(int iTimeBinMin, int iTimeBinMax)
{

   /*
   * vector to store all threads for parallel processing
   * one thread per CRU (360 in total)
   */
  std::vector<std::thread> thread_vector;
  if (mProcessingType == Processing::Parallel)
    LOG(DEBUG) << std::thread::hardware_concurrency() << " concurrent threads are supported." << FairLogger::endl;

  /*
   * if CRU number of current digit changes, start processing (either
   * sequential or parallel) of all CRUs in between the last processed
   * one and the current one.
   */

  const Mapper& mapper = Mapper::instance();
  for (int iCRU = mCRUMin; iCRU < mCRUMax; ++iCRU) {
    struct CfConfig cfConfig = {
      iCRU,
      mapper.getNumberOfRowsPartition(iCRU),
      mPadsMax+2+2,
      iTimeBinMin,
      iTimeBinMax,
      mEnableNoiseSim,
      mEnablePedestalSubtraction,
      mIsContinuousReadout,
      mNoiseObject,
      mPedestalObject
    };
//        std::cout << "starting CRU " << iCRU << std::endl;
    if (mProcessingType == Processing::Parallel)
      thread_vector.emplace_back(
            processDigits,                      // function name
            std::ref(mDigitContainer[iCRU]),    // digit container for individual CRUs
            std::ref(mClusterFinder[iCRU]),     // cluster finder for individual CRUs
            std::ref(mClusterStorage[iCRU]),    // container to store found clusters
            cfConfig
        );
    else {
      processDigits(
          std::ref(mDigitContainer[iCRU]),
          std::ref(mClusterFinder[iCRU]),
          std::ref(mClusterStorage[iCRU]),
          cfConfig);
    }
  }


  /*
   * wait for threads to join
   */
  for (std::thread& t: thread_vector) {
    t.join();
  }

  /*
   * collect clusters from individual cluster finder
   */
  for (auto& cc : mClusterStorage) {
    if (cc.size() != 0) {
      for (auto& c : cc){
	ClusterContainer::AddCluster<HwCluster>(mClusterArray, c.getCRU(),c.getRow(),c.getQ(),c.getQmax(),
            c.getPadMean(),c.getTimeMean(),c.getPadSigma(),c.getTimeSigma());
      }
    }
  }

  mLastTimebin = iTimeBinMax;
}
