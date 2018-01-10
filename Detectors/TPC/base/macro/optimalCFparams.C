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
#include <vector>
#include <sstream>
#include <memory>
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGaxis.h"
#endif

/*
gSystem->AddIncludePath("-I$O2_ROOT/include -I$FAIRROOT_ROOT/include");
.L createCRUSorterMapping.C+
*/

//using namespace o2::TPC;

unsigned getNumberCfPerRow(unsigned numPads, unsigned maxPads = 138);
unsigned getMaxClusters(unsigned numPads, unsigned timebins);
unsigned getProcessingTime(unsigned numPads, unsigned timebins, unsigned version = 1);
unsigned getNumberStoredWords(unsigned numPads, unsigned timebins);

void drawAndSave(TH1* hist, TFile* file = nullptr, int scaleY = 0, std::string title = "");
void drawAndSave(TH2* hist, TFile* file = nullptr);
void drawAndSave(THStack* hist, TFile* file = nullptr);
void adaptTitleWidth(float width);
void changeBackStyle();

void optimalCFparams()
{
//  TFile *outfile = new TFile("numberCFs.root","RECREATE");
  TH1I* hNumCFs = new TH1I("hNumCFs","Cluster finder instances; #pads per CF; #CF per row",40,0.5,40.5);
  for (int i=5; i<=40; ++i) {
    hNumCFs->SetBinContent(i,getNumberCfPerRow(i));
  }

  TH2I* hNumCC_v1 = new TH2I("hNumCC_v1", "Processing Time; #pads per CF; #timebins; needed clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2F* hFracCC_v1 = new TH2F("hFracCC_v1", "Usage of available processing time; #pads per CF; #timebins; used fraction of clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2I* hNumCC2_v1 = new TH2I("hNumCC2_v1", "Processing Time; #pads per CF; #timebins; needed clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2F* hFracCC2_v1 = new TH2F("hFracCC2_v1", "Usage of available processing time; #pads per CF; #timebins; used fraction of clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);

  TH2I* hNumCC_v2 = new TH2I("hNumCC_v2", "Processing Time; #pads per CF; #timebins; needed clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2F* hFracCC_v2 = new TH2F("hFracCC_v2", "Usage of available processing time; #pads per CF; #timebins; used fraction of clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2I* hNumCC2_v2 = new TH2I("hNumCC2_v2", "Processing Time; #pads per CF; #timebins; needed clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2F* hFracCC2_v2 = new TH2F("hFracCC2_v2", "Usage of available processing time; #pads per CF; #timebins; used fraction of clock cycles", 30, 0.5, 30.5, 30, 0.5, 30.5);
  int time;
  for (int timebins=5; timebins<=30; ++timebins) {
    for (int pads=5; pads<=30; ++pads) {
      std::cout << timebins << " timebins," << pads << " pads:" << std::endl;
      time = getProcessingTime(pads,timebins,1);
      std::cout << "total: " << time << ", available: " << ((timebins-4)*48) << " (stored bins: " << timebins*pads << ")" << std::endl << std::endl;
      hNumCC2_v1->SetBinContent(pads,timebins,time);
      hFracCC2_v1->SetBinContent(pads,timebins,static_cast<float>(time)/((timebins-4)*48));
      if (time > (timebins-4)*48) continue;
      hNumCC_v1->SetBinContent(pads,timebins,time);
      hFracCC_v1->SetBinContent(pads,timebins,static_cast<float>(time)/((timebins-4)*48));
    }
  }

  for (int timebins=5; timebins<=30; ++timebins) {
    for (int pads=5; pads<=30; ++pads) {
//      std::cout << timebins << " timebins," << pads << " pads:" << std::endl;
      time = getProcessingTime(pads,timebins,2);
//      std::cout << "total: " << time << ", available: " << ((timebins-4)*48) << " (stored bins: " << timebins*pads << ")" << std::endl << std::endl;
      hNumCC2_v2->SetBinContent(pads,timebins,time);
      hFracCC2_v2->SetBinContent(pads,timebins,static_cast<float>(time)/((timebins-4)*48));
      if (time > (timebins-4)*48) continue;
      hNumCC_v2->SetBinContent(pads,timebins,time);
      hFracCC_v2->SetBinContent(pads,timebins,static_cast<float>(time)/((timebins-4)*48));
    }
  }


  TH2I* hNumMlab = new TH2I("hNumMlab", "Needed MLABs; #pads per CF; #timebins; #MLABs", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2F* hFracMlab = new TH2F("hFracMlab", "Usage of MLABs; #pads per CF; #timebins; efficiency", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2I* hNumMlab2 = new TH2I("hNumMlab2", "Needed MLABs; #pads per CF; #timebins; # MLABs", 30, 0.5, 30.5, 30, 0.5, 30.5);
  TH2F* hFracMlab2 = new TH2F("hFracMlab2", "Usage of MLABs; #pads per CF; #timebins; efficiency", 30, 0.5, 30.5, 30, 0.5, 30.5);
  unsigned words;
  for (int timebins=5; timebins<=30; ++timebins) {
    for (int pads=5; pads<=30; ++pads) {
      std::cout << timebins << " timebins," << pads << " pads:" << std::endl;
      words = getNumberStoredWords(pads,timebins);
      float usedMlabs = ceil(static_cast<float>(words)/32);
      std::cout << "stored words: " << words << ", needed MLABs (per CF): " << usedMlabs << ", (total) " << usedMlabs*getNumberCfPerRow(pads)*18 << ", number of CFs per row: " << getNumberCfPerRow(pads) << std::endl;
      hNumMlab2->SetBinContent(pads,timebins,usedMlabs*getNumberCfPerRow(pads)*18);
      hFracMlab2->SetBinContent(pads,timebins,static_cast<float>(words)/(usedMlabs*32));
      if (usedMlabs*getNumberCfPerRow(pads)*18 > 5000) continue;
      hNumMlab->SetBinContent(pads,timebins,usedMlabs*getNumberCfPerRow(pads)*18);
      hFracMlab->SetBinContent(pads,timebins,static_cast<float>(words)/(usedMlabs*32));
    }
  }

  THStack *hs1 = new THStack("hs1","Stack 1");
  THStack *hs2 = new THStack("hs2","Stack 2");
  TH1I* hNumCC_6TB_peak = new TH1I("hNumCC_6TB_peak", "Processing Time contribution of peak finding; #pads; needed clock cycles", 30, 0.5, 30.5);
  TH1I* hNumCC_6TB_clus = new TH1I("hNumCC_6TB_clus", "Processing Time contribution of cluster calculation; #pads; needed clock cycles", 30, 0.5, 30.5);

  for (int pads=5; pads<=30; ++pads) {
    unsigned time = getProcessingTime(pads,6);
    unsigned clusContr = getMaxClusters(pads,6)*(25-9);
    hNumCC_6TB_peak->SetBinContent(pads,time-clusContr);
    hNumCC_6TB_clus->SetBinContent(pads,clusContr);
  }
  hNumCC_6TB_peak->SetFillColor(kGreen);
  hNumCC_6TB_clus->SetFillColor(kBlue);
  hs1->Add(hNumCC_6TB_clus);
  hs1->Add(hNumCC_6TB_peak);

  hs2->Add(hNumCC_6TB_peak);
  hs2->Add(hNumCC_6TB_clus);

  hFracCC_v1->GetZaxis()->SetRangeUser(0,1);
  hFracCC2_v1->GetZaxis()->SetRangeUser(0,1);
  hFracCC_v2->GetZaxis()->SetRangeUser(0,1);
  hFracCC2_v2->GetZaxis()->SetRangeUser(0,1);
  hFracMlab->GetZaxis()->SetRangeUser(0.5,1);
  hFracMlab2->GetZaxis()->SetRangeUser(0.5,1);
  drawAndSave(hNumCFs,nullptr,18,"#CF total");

  drawAndSave(hNumCC_v1);
  drawAndSave(hFracCC_v1);
  drawAndSave(hNumCC2_v1);
  drawAndSave(hFracCC2_v1);

  drawAndSave(hNumCC_v2);
  drawAndSave(hFracCC_v2);
  drawAndSave(hNumCC2_v2);
  drawAndSave(hFracCC2_v2);

  drawAndSave(hNumMlab);
  drawAndSave(hFracMlab);
  drawAndSave(hNumMlab2);
  drawAndSave(hFracMlab2);

  drawAndSave(hNumCC_6TB_peak);
  drawAndSave(hNumCC_6TB_clus);
  drawAndSave(hs1);
  drawAndSave(hs2);


//  outfile->Close();

}

void drawAndSave(TH1* hist, TFile* file = nullptr, int scaleY, std::string title) {
  TCanvas *c1 = new TCanvas();
  c1->cd();

  hist->SetStats(0);
  hist->Draw();
  c1->SetRightMargin(0.13);
  adaptTitleWidth(0.77);
  c1->Update();
  TGaxis *axis2;
  if (scaleY != 0) {
    axis2 = new TGaxis(
        c1->GetUxmax(),
        c1->GetUymin(),
        c1->GetUxmax(),
        c1->GetUymax(),
        c1->GetFrame()->GetY1(),scaleY*c1->GetFrame()->GetY2(),510,"+L");
    axis2->SetTitle(title.c_str());
    axis2->SetLabelOffset(0.01);
    axis2->SetTitleOffset(1.3);
    axis2->Draw("same");
    c1->Update();
  }
  c1->SaveAs((std::string(hist->GetName())+".pdf").c_str());

  if (file != nullptr) {
    file->cd();
    hist->Write();
    c1->Write();
  }
  changeBackStyle();
  delete c1;
}

void drawAndSave(TH2* hist, TFile* file = nullptr) {
  TCanvas *c1 = new TCanvas();
  c1->cd();

  hist->SetStats(0);
  hist->Draw("colz");
  hist->GetZaxis()->SetTitleOffset(1.3);
  c1->SetRightMargin(0.16);
  adaptTitleWidth(0.74);
  c1->Update();
  c1->SaveAs((std::string(hist->GetName())+".pdf").c_str());

  if (file != nullptr) {
    file->cd();
    hist->Write();
    c1->Write();
  }
  changeBackStyle();
  delete c1;
}

void drawAndSave(THStack* hist, TFile* file = nullptr) {
  TCanvas *c1 = new TCanvas();
  c1->cd();

  hist->Draw();
  c1->SetRightMargin(0.13);
  adaptTitleWidth(0.77);
  c1->Update();
  c1->SaveAs((std::string(hist->GetName())+".pdf").c_str());

  if (file != nullptr) {
    file->cd();
    hist->Write();
    c1->Write();
  }
  changeBackStyle();
  delete c1;
}

void adaptTitleWidth(float width) {
  TStyle *st = new TStyle("tmpstyle","tmpstyle");
  gROOT->GetStyle("mystyle")->Copy((*st));
  st->SetTitleW(width);
  st->cd();
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

void changeBackStyle(){
  gROOT->GetStyle("mystyle")->cd();
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

unsigned getNumberCfPerRow(unsigned numPads, unsigned maxPads) {
  maxPads += 2; // +2 empty on the left side
  maxPads += 2; // +2 empty on the right side
  if (numPads <= 4) return -1;
//  std::cout << std::endl << numPads << std::endl;
//  int right = 0;
//  int left = 0;
//  int i=0;
//  while (right < (padsMax-1)) {
//    ++i;
//    left = (i-1)*(numPads-4);
//    right = i*(numPads-4)+4-1;
//
////    std::cout << i << " left: " << left << " right: " << right << std::endl;
//  }
//  std::cout << i  << " " << ceil(static_cast<float>(padsMax - 4)/(numPads-4)) << std::endl;
////  std::cout << ceil(static_cast<float>(padsMax)/(numPads-4)) << std::endl;
  return ceil(static_cast<float>(maxPads-4)/(numPads-4));
}

unsigned getMaxClusters(unsigned numPads, unsigned timebins) {
  if (numPads <= 4 || timebins <= 4) return -1;

  return ceil(static_cast<float>(numPads-4)/2) * ceil(static_cast<float>(timebins-4)/2);
}

unsigned getProcessingTime(unsigned numPads, unsigned timebins, unsigned version) {
  if (numPads <= 4 || timebins <= 4) return -1;

//  std::cout << "Write new pads: " <<  (timebins-4) * numPads << std::endl;
//  std::cout << "Reading pads for cf search: " << 9 + (((timebins-4)*(numPads-4))-1)*3 << ", 9+" << (((timebins-4)*(numPads-4))-1)*3 << std::endl;
//  std::cout << "possible additional clusters: " << getMaxClusters(numPads,timebins)*(25-9) << std::endl;

  if (version == 1) {
    return
      9 + (timebins-5)*3 +                            // reading first pad for each timebin
      (numPads-5) * 3 +                               // reading new pads for "timebin 0"
      (numPads-5)*(timebins-5) * 3 +                  // reding new pad for other timebins
      getMaxClusters(numPads,timebins)*(25-9);
  } else if (version == 2) {
    return
      9 + (timebins-5)*3 +                            // reading first pad for each timebin
      (numPads-5) * 3 +                               // reading new pads for "timebin 0"
      (numPads-5)*(timebins-5) +                      // reding new pad for other timebins
      getMaxClusters(numPads,timebins)*(25-9);
  } else {
    return -1;
  }
//    9 + (((timebins-4)*(numPads-4))-1) * 3 +        // searching for cluster, first pad needs 9 cycles, every other only 3 new pads
//    getMaxClusters(numPads,timebins)*(25-9);        // reading remaining pads for found clusters

}

unsigned getNumberStoredWords(unsigned numPads, unsigned timebins) {
  if (numPads <= 4 || timebins <= 4) return -1;

  unsigned storedWords = (numPads * timebins // storage to work with
    + (numPads * (timebins-4)));             // buffer for new timebins

  return storedWords;
}
