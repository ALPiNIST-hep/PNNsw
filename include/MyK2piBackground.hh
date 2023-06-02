// ---------------------------------------------------------------------
// History:
//
// Created by Vic Marchetti (victor.marchetti@cern.ch) 2023-06-02
//
// ---------------------------------------------------------------------

/// \class MyK2piBackground
/// \Brief
/// Estimate of the K+ -> pi+ pi0(gamma) background under the PNN signal
/// regions R1/R2, from the radiative-tail leakage of the missing mass in
/// the pi0 hypothesis.
/// \EndBrief
/// \Detailed
/// A K2pi-enriched control sample is selected by pairing LKr photon
/// clusters into a pi0 candidate and requiring the pi+ pi0 missing mass to
/// be compatible with zero. The fraction of these control events that,
/// after photon energy loss (radiative tail, LKr acceptance edge, merged
/// clusters), migrate into the PNN squared-missing-mass regions R1/R2 is
/// measured and applied to the normalisation K2pi yield to obtain the
/// expected background.
/// \EndDetailed

#ifndef MYK2PIBACKGROUND_HH
#define MYK2PIBACKGROUND_HH

#include "Analyzer.hh"
#include "TH1F.h"
#include "TH2F.h"

class MyK2piBackground : public NA62Analysis::Analyzer {

 public:
  explicit MyK2piBackground(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyK2piBackground();
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

  // ---- cuts ---------------------------------------------------------
  Double_t fMaxPi0MissingMass2;   ///< |m2(pi+pi0) - 0| cut defining the K2pi control sample [GeV^2/c^4]
  Double_t fMinClusterEnergy;     ///< minimum LKr cluster energy for photon pairing [MeV]
  Double_t fPi0MassWindow;        ///< +/- window around m(pi0) for photon pairing [MeV/c^2]

  Double_t fR1Min, fR1Max;
  Double_t fR2Min, fR2Max;

  // ---- outputs --------------------------------------------------------
  Double_t fTailFractionR1;
  Double_t fTailFractionR2;
  Int_t fNControlEvents;
  Int_t fNTailEventsR1;
  Int_t fNTailEventsR2;

  // ---- histograms -----------------------------------------------------
  TH1F *fHMPi0;
  TH1F *fHM2MissPi0Hyp;
  TH1F *fHM2MissTail;
  TH2F *fHM2MissTailVsP;
};

#endif
