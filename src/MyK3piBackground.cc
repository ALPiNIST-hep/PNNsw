// ---------------------------------------------------------------------
// History:
//
// Created by Mika Ehrlich (michael.ehrlich@cern.ch) 2024-01-25
//
// ---------------------------------------------------------------------

#include "MyK3piBackground.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "KaonDecayConstants.hh"
#include "Persistency.hh"
#include "SpectrometerTrackVertex.hh"
#include "functions.hh"

#include <algorithm>
#include <iostream>

using namespace NA62Analysis;

MyK3piBackground::MyK3piBackground(Core::BaseAnalysis *ba) : Analyzer(ba, "MyK3piBackground") {
  fReadingData = GetIsTree();

  AddParam("Max3TrackVertexChi2", &fMax3TrackVertexChi2, 25.);
  AddParam("KaonMassWindow", &fKaonMassWindow, 20.);  // MeV/c^2

  AddParam("R1Min", &fR1Min, 0.000);
  AddParam("R1Max", &fR1Max, 0.010);
  AddParam("R2Min", &fR2Min, 0.026);
  AddParam("R2Max", &fR2Max, 0.068);

  fAcceptanceLossFractionR1 = 0.;
  fAcceptanceLossFractionR2 = 0.;
  fNControlEvents = 0;
  fNFakePNN_R1 = 0;
  fNFakePNN_R2 = 0;
}

void MyK3piBackground::InitOutput() {
  if (!fReadingData) return;
  RegisterOutput("AcceptanceLossFractionR1", &fAcceptanceLossFractionR1);
  RegisterOutput("AcceptanceLossFractionR2", &fAcceptanceLossFractionR2);
}

void MyK3piBackground::InitHist() {
  if (!fReadingData) return;

  fHM3Pi = new TH1F("hM3Pi", "3#pi invariant mass, K3pi control sample;"
                     "M_{3#pi} [MeV/c^{2}];Events", 200, 400., 600.);
  BookHisto(fHM3Pi);

  fHM2MissLeadingTrack = new TH1F("hM2MissLeadingTrack",
                                   "m^{2}_{miss}(#pi^{+}) for the leading track only, K3pi control;"
                                   "m^{2}_{miss} [GeV^{2}/c^{4}];Events", 280, -0.02, 0.12);
  BookHisto(fHM2MissLeadingTrack);
}

void MyK3piBackground::Process(int) {
  if (!fReadingData) return;

  // ------------------------------------------------------------------
  // Fully-reconstructed K3pi control sample: 3-track vertex, total
  // charge +1, invariant mass compatible with M(K+).
  // ------------------------------------------------------------------
  std::vector<SpectrometerTrackVertex> Vertices3Trk =
      *(std::vector<SpectrometerTrackVertex> *)GetOutput("SpectrometerVertexBuilder.Output3");
  if (Vertices3Trk.empty()) return;

  Int_t bestVtx = -1;
  Double_t bestChi2 = 1e9;
  for (UInt_t iV = 0; iV < Vertices3Trk.size(); iV++) {
    if (Vertices3Trk[iV].GetCharge() != 1) continue;
    Double_t chi2 = Vertices3Trk[iV].GetChi2();
    if (chi2 > fMax3TrackVertexChi2) continue;
    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestVtx = iV;
    }
  }
  if (bestVtx < 0) return;

  TLorentzVector P[3];
  std::vector<std::pair<Double_t, Int_t>> momByTrack;  // (momentum, local index 0..2)
  for (Int_t i = 0; i < 3; i++) {
    TVector3 p3 = Vertices3Trk[bestVtx].GetTrackThreeMomentum(i);
    P[i].SetVectM(p3, MPI);
    momByTrack.push_back(std::make_pair(p3.Mag(), i));
  }

  Double_t m3pi = (P[0] + P[1] + P[2]).M();
  fHM3Pi->Fill(m3pi);
  if (std::abs(m3pi - MKCH) > fKaonMassWindow) return;

  fNControlEvents++;

  // ------------------------------------------------------------------
  // Remove the two lowest-momentum tracks "by hand" -- these are the
  // ones most likely to fall out of the STRAW/RICH/CHOD acceptance or
  // overlap with the leading track in a real (non-control) event -- and
  // rerun the PNN missing-mass calculation on the surviving track only.
  // ------------------------------------------------------------------
  std::sort(momByTrack.begin(), momByTrack.end());
  Int_t leadingTrack = momByTrack.back().second;

  TVector3 BeamMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  TLorentzVector Kaon;
  Kaon.SetVectM(BeamMomentum, MKCH);

  Double_t m2missPNN = (Kaon - P[leadingTrack]).M2() * 1e-6;
  fHM2MissLeadingTrack->Fill(m2missPNN);

  if (m2missPNN > fR1Min && m2missPNN < fR1Max) fNFakePNN_R1++;
  if (m2missPNN > fR2Min && m2missPNN < fR2Max) fNFakePNN_R2++;

  if (fNControlEvents > 0) {
    fAcceptanceLossFractionR1 = static_cast<Double_t>(fNFakePNN_R1) / fNControlEvents;
    fAcceptanceLossFractionR2 = static_cast<Double_t>(fNFakePNN_R2) / fNControlEvents;
  }
}

void MyK3piBackground::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyK3piBackground::EndOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyK3piBackground::EndOfRunUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyK3piBackground::EndOfJobUser() {
  std::cout << "[MyK3piBackground] control sample: " << fNControlEvents
            << ", acceptance-loss fraction R1 = " << fAcceptanceLossFractionR1
            << ", R2 = " << fAcceptanceLossFractionR2 << std::endl;
  SaveAllPlots();
}

void MyK3piBackground::DrawPlot() {
  /// \MemberDescr
  /// \EndMemberDescr
}

MyK3piBackground::~MyK3piBackground() {
  /// \MemberDescr
  /// \EndMemberDescr
}
