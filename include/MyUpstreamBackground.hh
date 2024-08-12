// ---------------------------------------------------------------------
// History:
//
// Created by Vic Marchetti (victor.marchetti@cern.ch) 2024-08-12
//
// ---------------------------------------------------------------------

/// \class MyUpstreamBackground
/// \Brief
/// Estimate of the "upstream" PNN background, from pi+ tracks produced
/// by beam-material interactions or accidental beam-particle matching
/// upstream of the fiducial volume.
/// \EndBrief
/// \Detailed
/// The Upstream Reference Sample (URS) is built by removing the K/pi
/// time-matching requirement from the standard PNN vertex selection and
/// inverting the CDA cut (CDA > fCDAInvertMin). This strongly enhances
/// the two upstream production mechanisms -- beam-material interaction,
/// and an upstream K+ decay accidentally matched to an unrelated beam
/// particle -- while remaining statistically independent of the signal
/// region. The number of upstream background events under R1/R2 is
/// obtained by extrapolating the URS into the signal CDA range and
/// correcting for the removed K/pi matching requirement with the
/// mismatch probability fPMatch (see SetMismatchProbability()).
/// \EndDetailed

#ifndef MYUPSTREAMBACKGROUND_HH
#define MYUPSTREAMBACKGROUND_HH

#include "Analyzer.hh"
#include "TH1F.h"
#include "TH2F.h"

class MyUpstreamBackground : public NA62Analysis::Analyzer {

 public:
  explicit MyUpstreamBackground(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyUpstreamBackground();
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

  Double_t fCDAInvertMin;    ///< mm, lower edge of the URS (CDA > this value)
  Double_t fCDASignalMax;    ///< mm, upper edge of the signal CDA range
  Double_t fPMatch;          ///< accidental K/pi mismatch probability correction

  Double_t fR1Min, fR1Max;
  Double_t fR2Min, fR2Max;

  // ---- outputs ----------------------------------------------------------
  Double_t fURSYieldR1;
  Double_t fURSYieldR2;
  Double_t fExtrapolatedYieldR1;
  Double_t fExtrapolatedYieldR2;

  // ---- histograms -----------------------------------------------------
  TH1F *fHCDA_URS;
  TH2F *fHM2MissVsCDA;
  TH1F *fHM2Miss_URS;
};

#endif
