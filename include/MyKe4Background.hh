// ---------------------------------------------------------------------
// History:
//
// Created by Mika Ehrlich (michael.ehrlich@cern.ch) 2024-03-30
//
// ---------------------------------------------------------------------

/// \class MyKe4Background
/// \Brief
/// Estimate of the K+ -> pi+ pi- e+ nu background under the PNN signal
/// regions R1/R2.
/// \EndBrief
/// \Detailed
/// A Ke4 control sample is selected by requiring a 3-track vertex with
/// total charge +1, one track identified as an electron (E/p compatible
/// with 1) and one identified as a pi-. As for MyK3piBackground, the
/// electron and pi- are then discarded and the surviving pi+ track alone
/// is passed through the PNN missing-mass calculation, to measure how
/// often such an event would be picked up as a PNN candidate when the
/// other two tracks are lost.
/// \EndDetailed

#ifndef MYKE4BACKGROUND_HH
#define MYKE4BACKGROUND_HH

#include "Analyzer.hh"
#include "TH1F.h"

class MyKe4Background : public NA62Analysis::Analyzer {

 public:
  explicit MyKe4Background(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyKe4Background();
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

  Double_t fMinElectronEoP;
  Double_t fMax3TrackVertexChi2;

  Double_t fR1Min, fR1Max;
  Double_t fR2Min, fR2Max;

  // ---- outputs ----------------------------------------------------------
  Double_t fAcceptanceLossFractionR1;
  Double_t fAcceptanceLossFractionR2;
  Int_t fNControlEvents;
  Int_t fNFakePNN_R1;
  Int_t fNFakePNN_R2;

  // ---- histograms -----------------------------------------------------
  TH1F *fHEoPElectron;
  TH1F *fHM2MissLeadingPion;
};

#endif
