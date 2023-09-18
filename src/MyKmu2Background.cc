// ---------------------------------------------------------------------
// History:
//
// Created by Steve Weiss (stephan.weiss@cern.ch) 2023-09-18
//
// ---------------------------------------------------------------------

#include "MyKmu2Background.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "KaonDecayConstants.hh"
#include "Persistency.hh"
#include "SpectrometerTrackVertex.hh"
#include "functions.hh"

#include <iostream>

using namespace NA62Analysis;

MyKmu2Background::MyKmu2Background(Core::BaseAnalysis *ba) : Analyzer(ba, "MyKmu2Background") {
  fReadingData = GetIsTree();

  AddParam("MaxTrackChi2", &fMaxTrackChi2, 18.);
  AddParam("MaxEoP", &fMaxEoP, 0.85);
  AddParam("MinRICHLikelihoodPion", &fMinRICHLikelihoodPion, 0.90);

  AddParam("R1Min", &fR1Min, 0.000);
  AddParam("R1Max", &fR1Max, 0.010);
  AddParam("R2Min", &fR2Min, 0.026);
  AddParam("R2Max", &fR2Max, 0.068);

  fMisIDProbability = 0.;
  fNTaggedMuons = 0;
  fNMisIDAsPion = 0;
}

void MyKmu2Background::InitOutput() {
  if (!fReadingData) return;
  RegisterOutput("MisIDProbability", &fMisIDProbability);
}

void MyKmu2Background::InitHist() {
  if (!fReadingData) return;

  fHEoPTaggedMuon = new TH1F("hEoPTaggedMuon", "E/p for MUV3-tagged muons;E/p;Events", 150, 0., 1.5);
  BookHisto(fHEoPTaggedMuon);

  fHRICHLikelihoodTaggedMuon = new TH1F("hRICHLikelihoodTaggedMuon",
                                         "RICH single-#pi likelihood for MUV3-tagged muons;"
                                         "L_{#pi};Events", 100, 0., 1.);
  BookHisto(fHRICHLikelihoodTaggedMuon);

  fHM2MissKmu2Hyp = new TH1F("hM2MissKmu2Hyp", "m^{2}_{miss}(#mu^{+}) for tagged Kmu2 candidates;"
                              "m^{2}_{miss} [GeV^{2}/c^{4}];Events", 280, -0.02, 0.12);
  BookHisto(fHM2MissKmu2Hyp);

  fHM2MissMisIDPion = new TH1F("hM2MissMisIDPion",
                                "m^{2}_{miss}(#pi^{+} hyp.) for mis-identified Kmu2 candidates;"
                                "m^{2}_{miss} [GeV^{2}/c^{4}];Events", 280, -0.02, 0.12);
  BookHisto(fHM2MissMisIDPion);
}

void MyKmu2Background::Process(int) {
  if (!fReadingData) return;

  // ------------------------------------------------------------------
  // Tag: single-track vertex with an in-time MUV3 association. Because
  // the muon is independently confirmed by MUV3, the RICH/E-over-p
  // response of the *same* track can be used, unbiased, as a probe of
  // how often a genuine muon would be accepted as a pion by the PNN
  // selection cuts.
  // ------------------------------------------------------------------
  std::vector<SpectrometerTrackVertex> Vertices =
      *(std::vector<SpectrometerTrackVertex> *)GetOutput("SpectrometerVertexBuilder.Output1");
  if (Vertices.empty()) return;

  std::vector<DownstreamTrack> Tracks =
      *(std::vector<DownstreamTrack> *)GetOutput("DownstreamTrackBuilder.Output");

  for (UInt_t iV = 0; iV < Vertices.size(); iV++) {
    if (Vertices[iV].GetCharge() != 1) continue;
    if (Vertices[iV].GetChi2() > fMaxTrackChi2) continue;

    Int_t trackID = Vertices[iV].GetTrackIndex(0);
    if (!Tracks[trackID].MUV3AssociationExists()) continue;

    fNTaggedMuons++;

    TVector3 BeamMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
    TLorentzVector Kaon;
    Kaon.SetVectM(BeamMomentum, MKCH);
    TLorentzVector Muon;
    Muon.SetVectM(Vertices[iV].GetTrackThreeMomentum(0), MMU);
    Double_t m2missMu = (Kaon - Muon).M2() * 1e-6;
    fHM2MissKmu2Hyp->Fill(m2missMu);

    Double_t EoP = Tracks[trackID].GetLKrEoP();
    Double_t richLikelihoodPion = Tracks[trackID].GetRICHSinglePionProbability();
    fHEoPTaggedMuon->Fill(EoP);
    fHRICHLikelihoodTaggedMuon->Fill(richLikelihoodPion);

    // Would this genuine muon have survived the PNN pi+ PID cuts?
    if (EoP <= fMaxEoP && richLikelihoodPion >= fMinRICHLikelihoodPion) {
      fNMisIDAsPion++;
      TLorentzVector Pion;
      Pion.SetVectM(Vertices[iV].GetTrackThreeMomentum(0), MPI);
      Double_t m2missPi = (Kaon - Pion).M2() * 1e-6;
      fHM2MissMisIDPion->Fill(m2missPi);
    }
  }

  if (fNTaggedMuons > 0) fMisIDProbability = static_cast<Double_t>(fNMisIDAsPion) / fNTaggedMuons;
}

void MyKmu2Background::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyKmu2Background::EndOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyKmu2Background::EndOfRunUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyKmu2Background::EndOfJobUser() {
  std::cout << "[MyKmu2Background] tagged muons: " << fNTaggedMuons
            << ", mis-ID as pion: " << fNMisIDAsPion
            << ", mis-ID probability = " << fMisIDProbability << std::endl;
  SaveAllPlots();
}

void MyKmu2Background::DrawPlot() {
  /// \MemberDescr
  /// \EndMemberDescr
}

MyKmu2Background::~MyKmu2Background() {
  /// \MemberDescr
  /// \EndMemberDescr
}
