// ---------------------------------------------------------------------
// History:
//
// Created by Steve Weiss (stephan.weiss@cern.ch) 2025-02-14
//
// ---------------------------------------------------------------------

/// \class MyPNNBackgroundSummary
/// \Brief
/// Collects the per-channel background estimates from MyK2piBackground,
/// MyKmu2Background, MyK3piBackground, MyKe4Background and
/// MyUpstreamBackground, scales each to the PNN normalisation, and
/// writes out the combined background table for R1 and R2.
/// \EndBrief
/// \Detailed
/// This analyzer does not perform any event selection of its own: it
/// must run after MyPNNSelection and the five background analyzers in
/// the same job, and simply reads their outputs at the end of the job
/// (tail fractions / mis-ID probabilities / acceptance-loss fractions /
/// extrapolated URS yields), multiplies each by the corresponding
/// normalisation-channel yield (set via SetNormalisationYields()) and
/// writes the resulting table to fOutputTableFile.
/// \EndDetailed

#ifndef MYPNNBACKGROUNDSUMMARY_HH
#define MYPNNBACKGROUNDSUMMARY_HH

#include "Analyzer.hh"

#include <string>

class MyPNNBackgroundSummary : public NA62Analysis::Analyzer {

 public:
  explicit MyPNNBackgroundSummary(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyPNNBackgroundSummary();
  void InitOutput();
  void InitHist() {}
  void DefineMCSimple() {}
  void Process(int iEvent);
  void PostProcess();
  void EndOfBurstUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser();
  void DrawPlot() {}

 private:
  std::string fOutputTableFile;

  // Normalisation-channel yields, i.e. the number of fully reconstructed
  // K2pi / Kmu2 / K3pi / Ke4 events collected in the same dataset as the
  // PNN sample; set centrally once per data-taking period.
  Double_t fNormK2pi;
  Double_t fNormKmu2;
  Double_t fNormK3pi;
  Double_t fNormKe4;

  // ---- combined outputs -------------------------------------------------
  Double_t fBkgK2piR1, fBkgK2piR2;
  Double_t fBkgKmu2R1, fBkgKmu2R2;
  Double_t fBkgK3piR1, fBkgK3piR2;
  Double_t fBkgKe4R1, fBkgKe4R2;
  Double_t fBkgUpstreamR1, fBkgUpstreamR2;
  Double_t fBkgTotalR1, fBkgTotalR2;
};

#endif
