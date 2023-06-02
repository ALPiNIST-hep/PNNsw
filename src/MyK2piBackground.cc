// ---------------------------------------------------------------------
// History:
//
// Created by Vic Marchetti (victor.marchetti@cern.ch) 2023-06-02
//
// ---------------------------------------------------------------------

#include "MyK2piBackground.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "GeometricAcceptance.hh"
#include "KaonDecayConstants.hh"
#include "Persistency.hh"
#include "SpectrometerTrackVertex.hh"
#include "functions.hh"

#include <iostream>

using namespace NA62Analysis;

MyK2piBackground::MyK2piBackground(Core::BaseAnalysis *ba) : Analyzer(ba, "MyK2piBackground") {
  fReadingData = GetIsTree();

  AddParam("MaxPi0MissingMass2", &fMaxPi0MissingMass2, 0.0015);  // GeV^2/c^4
  AddParam("MinClusterEnergy", &fMinClusterEnergy, 3000.);       // MeV
  AddParam("Pi0MassWindow", &fPi0MassWindow, 15.);               // MeV/c^2

  AddParam("R1Min", &fR1Min, 0.000);
  AddParam("R1Max", &fR1Max, 0.010);
  AddParam("R2Min", &fR2Min, 0.026);
  AddParam("R2Max", &fR2Max, 0.068);

  fTailFractionR1 = 0.;
  fTailFractionR2 = 0.;
  fNControlEvents = 0;
  fNTailEventsR1 = 0;
  fNTailEventsR2 = 0;
}

void MyK2piBackground::InitOutput() {
  if (!fReadingData) return;
  RegisterOutput("TailFractionR1", &fTailFractionR1);
  RegisterOutput("TailFractionR2", &fTailFractionR2);
}

void MyK2piBackground::InitHist() {
  if (!fReadingData) return;

  fHMPi0 = new TH1F("hMPi0", "Two-photon invariant mass;M_{#gamma#gamma} [MeV/c^{2}];Events",
                     200, 0., 300.);
  BookHisto(fHMPi0);

  fHM2MissPi0Hyp = new TH1F("hM2MissPi0Hyp", "m^{2}(#pi^{+}#pi^{0}) in the K2pi control sample;"
                             "m^{2}_{#pi#pi^{0}} [GeV^{2}/c^{4}];Events", 200, -0.01, 0.01);
  BookHisto(fHM2MissPi0Hyp);

  fHM2MissTail = new TH1F("hM2MissTail", "m^{2}_{miss}(#pi^{+}) for K2pi control events;"
                           "m^{2}_{miss} [GeV^{2}/c^{4}];Events", 280, -0.02, 0.12);
  BookHisto(fHM2MissTail);

  fHM2MissTailVsP = new TH2F("hM2MissTailVsP", "m^{2}_{miss}(#pi^{+}) vs p_{#pi} for K2pi control events;"
                              "p_{#pi} [GeV/c];m^{2}_{miss} [GeV^{2}/c^{4}]", 60, 10., 50., 140, -0.02, 0.12);
  BookHisto(fHM2MissTailVsP);
}

void MyK2piBackground::Process(int) {
  if (!fReadingData) return;

  fTailFractionR1 = 0.;
  fTailFractionR2 = 0.;

  // ------------------------------------------------------------------
  // Single positive track + reconstructed pi0 -> K2pi control sample.
  // The pi0 is built from the best pi0-mass-compatible pair of LKr
  // clusters not associated to the charged track.
  // ------------------------------------------------------------------
  std::vector<SpectrometerTrackVertex> Vertices =
      *(std::vector<SpectrometerTrackVertex> *)GetOutput("SpectrometerVertexBuilder.Output1");
  if (Vertices.empty()) return;

  std::vector<EnergyCluster> *Clusters = (std::vector<EnergyCluster> *)GetOutput("EnergyClusterBuilder.Output");
  if (!Clusters || Clusters->size() < 2) return;

  Int_t bestVtx = -1;
  Double_t bestChi2 = 1e9;
  for (UInt_t iV = 0; iV < Vertices.size(); iV++) {
    if (Vertices[iV].GetCharge() != 1) continue;
    Double_t chi2 = Vertices[iV].GetChi2();
    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestVtx = iV;
    }
  }
  if (bestVtx < 0) return;

  // Best pi0 candidate among all photon pairs
  Double_t bestMassDiff = 1e9;
  Double_t bestMPi0 = -1.;
  TLorentzVector Pi0Candidate;
  for (UInt_t i = 0; i < Clusters->size(); i++) {
    if ((*Clusters)[i].GetEnergy() < fMinClusterEnergy) continue;
    for (UInt_t j = i + 1; j < Clusters->size(); j++) {
      if ((*Clusters)[j].GetEnergy() < fMinClusterEnergy) continue;
      TLorentzVector g1, g2;
      g1.SetVectM((*Clusters)[i].GetPosition() - Vertices[bestVtx].GetPosition(), 0);
      g1.SetE((*Clusters)[i].GetEnergy());
      g2.SetVectM((*Clusters)[j].GetPosition() - Vertices[bestVtx].GetPosition(), 0);
      g2.SetE((*Clusters)[j].GetEnergy());
      Double_t mgg = (g1 + g2).M();
      if (std::abs(mgg - MP0) < bestMassDiff) {
        bestMassDiff = std::abs(mgg - MP0);
        bestMPi0 = mgg;
        Pi0Candidate = g1 + g2;
      }
    }
  }
  if (bestMPi0 < 0 || bestMassDiff > fPi0MassWindow) return;
  fHMPi0->Fill(bestMPi0);

  TVector3 BeamMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  TLorentzVector Kaon;
  Kaon.SetVectM(BeamMomentum, MKCH);
  TLorentzVector Pion;
  Pion.SetVectM(Vertices[bestVtx].GetTrackThreeMomentum(0), MPI);

  Double_t m2missPi0 = (Kaon - Pion - Pi0Candidate).M2() * 1e-6;
  fHM2MissPi0Hyp->Fill(m2missPi0);
  if (std::abs(m2missPi0) > fMaxPi0MissingMass2) return;

  // This event is a clean K2pi control-sample candidate: now look at
  // where it falls in the PNN (pi+ only, ignoring the pi0) missing-mass
  // spectrum -- this is the radiative tail we need to quantify.
  fNControlEvents++;

  Double_t m2missPNN = (Kaon - Pion).M2() * 1e-6;
  Double_t pPion = Pion.Vect().Mag() * 1e-3;

  fHM2MissTail->Fill(m2missPNN);
  fHM2MissTailVsP->Fill(pPion, m2missPNN);

  if (m2missPNN > fR1Min && m2missPNN < fR1Max) fNTailEventsR1++;
  if (m2missPNN > fR2Min && m2missPNN < fR2Max) fNTailEventsR2++;

  if (fNControlEvents > 0) {
    fTailFractionR1 = static_cast<Double_t>(fNTailEventsR1) / fNControlEvents;
    fTailFractionR2 = static_cast<Double_t>(fNTailEventsR2) / fNControlEvents;
  }
}

void MyK2piBackground::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyK2piBackground::EndOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyK2piBackground::EndOfRunUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyK2piBackground::EndOfJobUser() {
  std::cout << "[MyK2piBackground] control sample: " << fNControlEvents
            << ", tail fraction R1 = " << fTailFractionR1
            << ", tail fraction R2 = " << fTailFractionR2 << std::endl;
  SaveAllPlots();
}

void MyK2piBackground::DrawPlot() {
  /// \MemberDescr
  /// \EndMemberDescr
}

MyK2piBackground::~MyK2piBackground() {
  /// \MemberDescr
  /// \EndMemberDescr
}
