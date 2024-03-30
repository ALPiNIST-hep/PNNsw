// ---------------------------------------------------------------------
// History:
//
// Created by Mika Ehrlich (michael.ehrlich@cern.ch) 2024-03-30
//
// ---------------------------------------------------------------------

#include "MyKe4Background.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "KaonDecayConstants.hh"
#include "Persistency.hh"
#include "SpectrometerTrackVertex.hh"
#include "functions.hh"

#include <iostream>

using namespace NA62Analysis;

MyKe4Background::MyKe4Background(Core::BaseAnalysis *ba) : Analyzer(ba, "MyKe4Background") {
  fReadingData = GetIsTree();

  AddParam("MinElectronEoP", &fMinElectronEoP, 0.90);
  AddParam("Max3TrackVertexChi2", &fMax3TrackVertexChi2, 25.);

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

void MyKe4Background::InitOutput() {
  if (!fReadingData) return;
  RegisterOutput("AcceptanceLossFractionR1", &fAcceptanceLossFractionR1);
  RegisterOutput("AcceptanceLossFractionR2", &fAcceptanceLossFractionR2);
}

void MyKe4Background::InitHist() {
  if (!fReadingData) return;

  fHEoPElectron = new TH1F("hEoPElectron", "E/p of the identified electron, Ke4 control;E/p;Events",
                            150, 0., 1.5);
  BookHisto(fHEoPElectron);

  fHM2MissLeadingPion = new TH1F("hM2MissLeadingPion",
                                  "m^{2}_{miss}(#pi^{+}) for the leading track only, Ke4 control;"
                                  "m^{2}_{miss} [GeV^{2}/c^{4}];Events", 280, -0.02, 0.12);
  BookHisto(fHM2MissLeadingPion);
}

void MyKe4Background::Process(int) {
  if (!fReadingData) return;

  // ------------------------------------------------------------------
  // 3-track vertex, total charge +1, with one track identified as an
  // electron by the calorimetric E/p response.
  // ------------------------------------------------------------------
  std::vector<SpectrometerTrackVertex> Vertices3Trk =
      *(std::vector<SpectrometerTrackVertex> *)GetOutput("SpectrometerVertexBuilder.Output3");
  if (Vertices3Trk.empty()) return;

  std::vector<DownstreamTrack> Tracks =
      *(std::vector<DownstreamTrack> *)GetOutput("DownstreamTrackBuilder.Output");

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

  Int_t electronTrack = -1;
  Double_t bestEoP = -1.;
  Int_t leadingPionTrack = -1;
  Double_t leadingPionMom = -1.;

  for (Int_t i = 0; i < 3; i++) {
    Int_t trackID = Vertices3Trk[bestVtx].GetTrackIndex(i);
    Double_t EoP = Tracks[trackID].GetLKrEoP();
    if (EoP > bestEoP) {
      bestEoP = EoP;
      electronTrack = i;
    }
  }
  if (electronTrack < 0 || bestEoP < fMinElectronEoP) return;
  fHEoPElectron->Fill(bestEoP);

  // Positive-charge, non-electron track with the highest momentum is
  // taken as the surviving "leading pion" if the other two are lost.
  for (Int_t i = 0; i < 3; i++) {
    if (i == electronTrack) continue;
    if (Vertices3Trk[bestVtx].GetTrackCharge(i) != 1) continue;
    Double_t p = Vertices3Trk[bestVtx].GetTrackThreeMomentum(i).Mag();
    if (p > leadingPionMom) {
      leadingPionMom = p;
      leadingPionTrack = i;
    }
  }
  if (leadingPionTrack < 0) return;

  fNControlEvents++;

  TVector3 BeamMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  TLorentzVector Kaon;
  Kaon.SetVectM(BeamMomentum, MKCH);
  TLorentzVector Pion;
  Pion.SetVectM(Vertices3Trk[bestVtx].GetTrackThreeMomentum(leadingPionTrack), MPI);

  Double_t m2missPNN = (Kaon - Pion).M2() * 1e-6;
  fHM2MissLeadingPion->Fill(m2missPNN);

  if (m2missPNN > fR1Min && m2missPNN < fR1Max) fNFakePNN_R1++;
  if (m2missPNN > fR2Min && m2missPNN < fR2Max) fNFakePNN_R2++;

  if (fNControlEvents > 0) {
    fAcceptanceLossFractionR1 = static_cast<Double_t>(fNFakePNN_R1) / fNControlEvents;
    fAcceptanceLossFractionR2 = static_cast<Double_t>(fNFakePNN_R2) / fNControlEvents;
  }
}

void MyKe4Background::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyKe4Background::EndOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyKe4Background::EndOfRunUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyKe4Background::EndOfJobUser() {
  std::cout << "[MyKe4Background] control sample: " << fNControlEvents
            << ", acceptance-loss fraction R1 = " << fAcceptanceLossFractionR1
            << ", R2 = " << fAcceptanceLossFractionR2 << std::endl;
  SaveAllPlots();
}

void MyKe4Background::DrawPlot() {
  /// \MemberDescr
  /// \EndMemberDescr
}

MyKe4Background::~MyKe4Background() {
  /// \MemberDescr
  /// \EndMemberDescr
}
