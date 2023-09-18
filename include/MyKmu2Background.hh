// ---------------------------------------------------------------------
// History:
//
// Created by Steve Weiss (stephan.weiss@cern.ch) 2023-09-18
//
// ---------------------------------------------------------------------

/// \class MyKmu2Background
/// \Brief
/// Estimate of the K+ -> mu+ nu(gamma) background under the PNN signal
/// regions R1/R2, from mu+ -> pi+ calorimetric/RICH mis-identification.
/// \EndBrief
/// \Detailed
/// A Kmu2-enriched control sample is selected by requiring a two-body
/// kinematic vertex with an in-time MUV3 association (tagged muon). The
/// probability that such a track would nevertheless pass the pi+
/// hypothesis particle-identification cuts used by MyPNNSelection
/// (E/p, RICH single-pion likelihood) is measured with a tag-and-probe
/// method and folded with the known two-body Kmu2 missing-mass line
/// shape to obtain the expected background in R1/R2.
/// \EndDetailed

#ifndef MYKMU2BACKGROUND_HH
#define MYKMU2BACKGROUND_HH

#include "Analyzer.hh"
#include "TH1F.h"

class MyKmu2Background : public NA62Analysis::Analyzer {

 public:
  explicit MyKmu2Background(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyKmu2Background();
  void InitOutput();
  void InitHist();
  void DefineMCSimple() {}
  void Process(int iEvent);
  void PostProcess();
  void EndOfBurstUser();
  void EndOfRunUser();
  void EndOfJobUser();
  void DrawPlot();

 private:
  Bool_t fReadingData;

  Double_t fMaxTrackChi2;
  Double_t fMaxEoP;                  ///< same cut as MyPNNSelection, for the mis-ID probe
  Double_t fMinRICHLikelihoodPion;   ///< same cut as MyPNNSelection, for the mis-ID probe

  Double_t fR1Min, fR1Max;
  Double_t fR2Min, fR2Max;

  // ---- outputs ----------------------------------------------------------
  Double_t fMisIDProbability;
  Int_t fNTaggedMuons;
  Int_t fNMisIDAsPion;

  // ---- histograms -----------------------------------------------------
  TH1F *fHEoPTaggedMuon;
  TH1F *fHRICHLikelihoodTaggedMuon;
  TH1F *fHM2MissKmu2Hyp;
  TH1F *fHM2MissMisIDPion;
};

#endif
