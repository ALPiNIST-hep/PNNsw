// ---------------------------------------------------------------------
// History:
//
// Created by Mika Ehrlich (michael.ehrlich@cern.ch) 2024-01-25
//
// ---------------------------------------------------------------------

/// \class MyK3piBackground
/// \Brief
/// Estimate of the K+ -> pi+ pi+ pi- background under the PNN signal
/// regions R1/R2, from events where two of the three tracks are lost
/// (out of acceptance, overlap, or reconstruction failure).
/// \EndBrief
/// \Detailed
/// A fully-reconstructed K3pi control sample (three-track vertex) is
/// selected first. For each control event, the two lowest-momentum
/// tracks are removed "by hand" and the surviving track is re-run
/// through the same missing-mass calculation as MyPNNSelection, to
/// measure how often such partially-reconstructed K3pi decays would
/// mimic a PNN candidate. This acceptance-loss probability is then
/// applied to the K3pi normalisation yield.
/// \EndDetailed

#ifndef MYK3PIBACKGROUND_HH
#define MYK3PIBACKGROUND_HH

#include "Analyzer.hh"
#include "TH1F.h"

class MyK3piBackground : public NA62Analysis::Analyzer {

 public:
  explicit MyK3piBackground(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyK3piBackground();
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

  Double_t fMax3TrackVertexChi2;
  Double_t fKaonMassWindow;   ///< +/- window around M(K+) for the 3-track invariant mass [MeV/c^2]

  Double_t fR1Min, fR1Max;
  Double_t fR2Min, fR2Max;

  // ---- outputs ----------------------------------------------------------
  Double_t fAcceptanceLossFractionR1;
  Double_t fAcceptanceLossFractionR2;
  Int_t fNControlEvents;
  Int_t fNFakePNN_R1;
  Int_t fNFakePNN_R2;

  // ---- histograms -----------------------------------------------------
  TH1F *fHM3Pi;
  TH1F *fHM2MissLeadingTrack;
};

#endif
