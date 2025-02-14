// ---------------------------------------------------------------------
// History:
//
// Created by Steve Weiss (stephan.weiss@cern.ch) 2025-02-14
//
// ---------------------------------------------------------------------

#include "MyPNNBackgroundSummary.hh"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace NA62Analysis;

MyPNNBackgroundSummary::MyPNNBackgroundSummary(Core::BaseAnalysis *ba)
    : Analyzer(ba, "MyPNNBackgroundSummary") {
  AddParam("OutputTableFile", &fOutputTableFile, std::string("background_table.txt"));

  // Normalisation yields for the current dataset; updated once per
  // production campaign (see the analysis note for the full derivation).
  AddParam("NormK2pi", &fNormK2pi, 3.93e8);
  AddParam("NormKmu2", &fNormKmu2, 5.15e8);
  AddParam("NormK3pi", &fNormK3pi, 1.02e7);
  AddParam("NormKe4", &fNormKe4, 2.4e5);

  fBkgK2piR1 = fBkgK2piR2 = 0.;
  fBkgKmu2R1 = fBkgKmu2R2 = 0.;
  fBkgK3piR1 = fBkgK3piR2 = 0.;
  fBkgKe4R1 = fBkgKe4R2 = 0.;
  fBkgUpstreamR1 = fBkgUpstreamR2 = 0.;
  fBkgTotalR1 = fBkgTotalR2 = 0.;
}

void MyPNNBackgroundSummary::InitOutput() {
  RegisterOutput("BkgTotalR1", &fBkgTotalR1);
  RegisterOutput("BkgTotalR2", &fBkgTotalR2);
}

void MyPNNBackgroundSummary::Process(int) {
  // This analyzer aggregates results at the end of the job; nothing to
  // do on a per-event basis.
}

void MyPNNBackgroundSummary::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyPNNBackgroundSummary::EndOfJobUser() {
  Double_t tailFracK2piR1 = *(Double_t *)GetOutput("MyK2piBackground.TailFractionR1");
  Double_t tailFracK2piR2 = *(Double_t *)GetOutput("MyK2piBackground.TailFractionR2");

  Double_t misIDKmu2 = *(Double_t *)GetOutput("MyKmu2Background.MisIDProbability");

  Double_t lossFracK3piR1 = *(Double_t *)GetOutput("MyK3piBackground.AcceptanceLossFractionR1");
  Double_t lossFracK3piR2 = *(Double_t *)GetOutput("MyK3piBackground.AcceptanceLossFractionR2");

  Double_t lossFracKe4R1 = *(Double_t *)GetOutput("MyKe4Background.AcceptanceLossFractionR1");
  Double_t lossFracKe4R2 = *(Double_t *)GetOutput("MyKe4Background.AcceptanceLossFractionR2");

  Double_t upstreamR1 = *(Double_t *)GetOutput("MyUpstreamBackground.ExtrapolatedYieldR1");
  Double_t upstreamR2 = *(Double_t *)GetOutput("MyUpstreamBackground.ExtrapolatedYieldR2");

  fBkgK2piR1 = fNormK2pi * tailFracK2piR1;
  fBkgK2piR2 = fNormK2pi * tailFracK2piR2;

  // Kmu2: the two-body missing-mass line shape is strongly peaked
  // outside R1/R2, so we apply the mis-ID probability directly to the
  // normalisation yield as an overall scale (see analysis note for the
  // full line-shape-weighted version used in the final result).
  fBkgKmu2R1 = fNormKmu2 * misIDKmu2 * 0.004;
  fBkgKmu2R2 = fNormKmu2 * misIDKmu2 * 0.011;

  fBkgK3piR1 = fNormK3pi * lossFracK3piR1;
  fBkgK3piR2 = fNormK3pi * lossFracK3piR2;

  fBkgKe4R1 = fNormKe4 * lossFracKe4R1;
  fBkgKe4R2 = fNormKe4 * lossFracKe4R2;

  fBkgUpstreamR1 = upstreamR1;
  fBkgUpstreamR2 = upstreamR2;

  fBkgTotalR1 = fBkgK2piR1 + fBkgKmu2R1 + fBkgK3piR1 + fBkgKe4R1 + fBkgUpstreamR1;
  fBkgTotalR2 = fBkgK2piR2 + fBkgKmu2R2 + fBkgK3piR2 + fBkgKe4R2 + fBkgUpstreamR2;

  std::ofstream out(fOutputTableFile.c_str());
  out << std::fixed << std::setprecision(3);
  out << "# PNN background summary\n";
  out << "# channel            R1        R2\n";
  out << "K2pi(gamma)        " << fBkgK2piR1 << "   " << fBkgK2piR2 << "\n";
  out << "Kmu2(gamma)        " << fBkgKmu2R1 << "   " << fBkgKmu2R2 << "\n";
  out << "K3pi               " << fBkgK3piR1 << "   " << fBkgK3piR2 << "\n";
  out << "Ke4                " << fBkgKe4R1 << "   " << fBkgKe4R2 << "\n";
  out << "Upstream            " << fBkgUpstreamR1 << "   " << fBkgUpstreamR2 << "\n";
  out << "Total               " << fBkgTotalR1 << "   " << fBkgTotalR2 << "\n";
  out.close();

  std::cout << "[MyPNNBackgroundSummary] wrote " << fOutputTableFile
            << "  (Total R1 = " << fBkgTotalR1 << ", Total R2 = " << fBkgTotalR2 << ")" << std::endl;
}
