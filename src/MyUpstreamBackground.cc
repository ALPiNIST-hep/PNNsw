// ---------------------------------------------------------------------
// History:
//
// Created by Vic Marchetti (victor.marchetti@cern.ch) 2024-08-12
//
// ---------------------------------------------------------------------

#include "MyUpstreamBackground.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "KaonDecayConstants.hh"
#include "Persistency.hh"
#include "SpectrometerTrackVertex.hh"
#include "functions.hh"

#include <iostream>

using namespace NA62Analysis;

MyUpstreamBackground::MyUpstreamBackground(Core::BaseAnalysis *ba)
    : Analyzer(ba, "MyUpstreamBackground") {
  fReadingData = GetIsTree();

  // The K/pi time-matching requirement is applied upstream, by the
  // vertex builder; here we simply widen the CDA window and let it
  // through unmatched (see class documentation for the rationale).
  AddParam("CDAInvertMin", &fCDAInvertMin, 4.0);   // mm
  AddParam("CDASignalMax", &fCDASignalMax, 30.0);  // mm, matches MyPNNSelection::fMaxCDA

  // Correction factor recovering the effect of removing the K/pi
  // matching requirement, evaluated once centrally and fixed here.
  AddParam("PMatch", &fPMatch, 0.62);

  AddParam("R1Min", &fR1Min, 0.000);
  AddParam("R1Max", &fR1Max, 0.010);
  AddParam("R2Min", &fR2Min, 0.026);
  AddParam("R2Max", &fR2Max, 0.068);

  fURSYieldR1 = 0.;
  fURSYieldR2 = 0.;
  fExtrapolatedYieldR1 = 0.;
  fExtrapolatedYieldR2 = 0.;
}

void MyUpstreamBackground::InitOutput() {
  if (!fReadingData) return;
  RegisterOutput("ExtrapolatedYieldR1", &fExtrapolatedYieldR1);
  RegisterOutput("ExtrapolatedYieldR2", &fExtrapolatedYieldR2);
}

void MyUpstreamBackground::InitHist() {
  if (!fReadingData) return;

  fHCDA_URS = new TH1F("hCDA_URS", "Track-beam CDA, Upstream Reference Sample;CDA [mm];Events",
                        300, 0., 150.);
  BookHisto(fHCDA_URS);

  fHM2MissVsCDA = new TH2F("hM2MissVsCDA", "m^{2}_{miss} vs CDA, URS;CDA [mm];m^{2}_{miss} [GeV^{2}/c^{4}]",
                            300, 0., 150., 140, -0.02, 0.12);
  BookHisto(fHM2MissVsCDA);

  fHM2Miss_URS = new TH1F("hM2Miss_URS", "m^{2}_{miss}, Upstream Reference Sample;"
                           "m^{2}_{miss} [GeV^{2}/c^{4}];Events", 280, -0.02, 0.12);
  BookHisto(fHM2Miss_URS);
}

void MyUpstreamBackground::Process(int) {
  if (!fReadingData) return;

  // ------------------------------------------------------------------
  // Upstream Reference Sample: identical vertex/track selection to
  // MyPNNSelection, but with the CDA cut inverted (CDA > fCDAInvertMin
  // instead of CDA < fCDASignalMax) and no K/pi timing requirement.
  // ------------------------------------------------------------------
  std::vector<SpectrometerTrackVertex> Vertices =
      *(std::vector<SpectrometerTrackVertex> *)GetOutput("SpectrometerVertexBuilder.Output1");
  if (Vertices.empty()) return;

  for (UInt_t iV = 0; iV < Vertices.size(); iV++) {
    if (Vertices[iV].GetCharge() != 1) continue;

    Double_t CDA = Vertices[iV].GetCDA();
    fHCDA_URS->Fill(CDA);
    if (CDA < fCDAInvertMin) continue;  // this is the URS-defining inversion

    Double_t p = Vertices[iV].GetTrackThreeMomentum(0).Mag();
    if (p < 15000. || p > 45000.) continue;

    TVector3 BeamMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
    TLorentzVector Kaon;
    Kaon.SetVectM(BeamMomentum, MKCH);
    TLorentzVector Pion;
    Pion.SetVectM(Vertices[iV].GetTrackThreeMomentum(0), MPI);

    Double_t m2miss = (Kaon - Pion).M2() * 1e-6;
    fHM2MissVsCDA->Fill(CDA, m2miss);
    fHM2Miss_URS->Fill(m2miss);

    if (m2miss > fR1Min && m2miss < fR1Max) fURSYieldR1 += 1.;
    if (m2miss > fR2Min && m2miss < fR2Max) fURSYieldR2 += 1.;
  }

  // ------------------------------------------------------------------
  // Extrapolation of the URS into the signal CDA range, corrected for
  // the removed K/pi matching requirement. The extrapolation factor
  // itself (CDA-shape ratio between the signal and URS regions) is
  // evaluated centrally from the template fit described in the class
  // documentation; a flat placeholder value is used here so that this
  // analyzer can run standalone during development.
  // ------------------------------------------------------------------
  const Double_t kCDAExtrapolationFactor = 0.083;
  fExtrapolatedYieldR1 = fURSYieldR1 * kCDAExtrapolationFactor * fPMatch;
  fExtrapolatedYieldR2 = fURSYieldR2 * kCDAExtrapolationFactor * fPMatch;
}

void MyUpstreamBackground::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyUpstreamBackground::EndOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyUpstreamBackground::EndOfRunUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyUpstreamBackground::EndOfJobUser() {
  std::cout << "[MyUpstreamBackground] URS yield R1 = " << fURSYieldR1
            << ", R2 = " << fURSYieldR2
            << " -> extrapolated R1 = " << fExtrapolatedYieldR1
            << ", R2 = " << fExtrapolatedYieldR2 << std::endl;
  SaveAllPlots();
}

void MyUpstreamBackground::DrawPlot() {
  /// \MemberDescr
  /// \EndMemberDescr
}

MyUpstreamBackground::~MyUpstreamBackground() {
  /// \MemberDescr
  /// \EndMemberDescr
}
