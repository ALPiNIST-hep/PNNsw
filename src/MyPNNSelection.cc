// ---------------------------------------------------------------------
// History:
//
// Created by Vic Marchetti (victor.marchetti@cern.ch) 2023-04-11
//
// ---------------------------------------------------------------------

#include "MyPNNSelection.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "GeometricAcceptance.hh"
#include "KaonDecayConstants.hh"
#include "Persistency.hh"
#include "SpectrometerTrackVertex.hh"
#include "TriggerConditions.hh"
#include "functions.hh"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace NA62Analysis;

MyPNNSelection::MyPNNSelection(Core::BaseAnalysis *ba) : Analyzer(ba, "MyPNNSelection") {
  fReadingData = GetIsTree();

  RequestTree("Cedar", new TRecoCedarEvent, "Reco");
  RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
  RequestL0Data();

  // ---- burst quality -------------------------------------------------
  // Maintained centrally; see config/BadBurstList.dat for the current
  // list and its revision history.
  AddParam("BadBurstListFile", &fBadBurstListFile, std::string("config/BadBurstList.dat"));

  // ---- trigger ---------------------------------------------------------
  AddParam("TriggerMask", &fTriggerMask, 0xFFFFu);

  // ---- track quality -----------------------------------------------------
  AddParam("MinTrackMomentum", &fMinTrackMomentum, 15000.);  // MeV/c
  AddParam("MaxTrackMomentum", &fMaxTrackMomentum, 45000.);  // MeV/c
  AddParam("MaxTrackChi2", &fMaxTrackChi2, 18.);  // tightened after 2023 STRAW alignment update

  // ---- vertex quality ------------------------------------------------
  AddParam("ZVertexMin", &fZVertexMin, 105000.);  // mm
  AddParam("ZVertexMax", &fZVertexMax, 165000.);  // mm
  AddParam("MaxCDA", &fMaxCDA, 30.);              // mm

  // ---- particle identification ----------------------------------------
  AddParam("MinRICHLikelihoodPion", &fMinRICHLikelihoodPion, 0.90);  // raised after RICH mirror recalibration
  AddParam("MaxEoP", &fMaxEoP, 0.85);
  AddParam("RequireMUV3NoAssociation", &fRequireMUV3NoAssociation, true);

  // ---- missing-mass signal regions [GeV^2/c^4] --------------------------
  // R1: below the K+ -> pi+ pi0 exclusion band.
  // R2: above the K+ -> pi+ pi0 exclusion band, below the 3-track threshold.
  AddParam("R1Min", &fR1Min, 0.000);
  AddParam("R1Max", &fR1Max, 0.010);
  AddParam("R2Min", &fR2Min, 0.026);
  AddParam("R2Max", &fR2Max, 0.068);

  fPNNEventSelected = false;
  fPionTrackID = -1;
  fM2Miss = -9999.;
  fSignalRegion = 0;
}

void MyPNNSelection::InitOutput() {
  if (!fReadingData) return;
  RegisterOutput("PNNEventSelected", &fPNNEventSelected);
  RegisterOutput("PionTrackID", &fPionTrackID);
  RegisterOutput("M2Miss", &fM2Miss);
  RegisterOutput("SignalRegion", &fSignalRegion);
}

void MyPNNSelection::InitHist() {
  if (!fReadingData) return;

  fHM2Miss = new TH1F("hM2Miss", "Squared missing mass;m^{2}_{miss} [GeV^{2}/c^{4}];Events",
                       280, -0.02, 0.12);
  BookHisto(fHM2Miss);

  fHM2MissVsP = new TH2F("hM2MissVsP", "m^{2}_{miss} vs p_{#pi};p_{#pi} [GeV/c];m^{2}_{miss} [GeV^{2}/c^{4}]",
                          60, 10., 50., 140, -0.02, 0.12);
  BookHisto(fHM2MissVsP);

  fHZVertex = new TH1F("hZVertex", "Track-beam vertex Z;Z [mm];Events", 300, 90000., 180000.);
  BookHisto(fHZVertex);

  fHCDA = new TH1F("hCDA", "Track-beam CDA;CDA [mm];Events", 200, 0., 100.);
  BookHisto(fHCDA);

  fHEoP = new TH1F("hEoP", "E/p;E/p;Events", 150, 0., 1.5);
  BookHisto(fHEoP);

  fHEventsPerBurst = new TH1F("hEventsPerBurst", "Processed events per burst;Burst ID;Events",
                               4000, 0., 4000.);
  BookHisto(fHEventsPerBurst);

  fHEventsPerBurstAfterQuality =
      new TH1F("hEventsPerBurstAfterQuality", "Processed events per burst, good bursts only;Burst ID;Events",
                4000, 0., 4000.);
  BookHisto(fHEventsPerBurstAfterQuality);

  LoadBadBurstList();
}

void MyPNNSelection::LoadBadBurstList() {
  fBadBursts.clear();
  std::ifstream in(fBadBurstListFile.c_str());
  if (!in.is_open()) {
    std::cerr << "[MyPNNSelection] WARNING: could not open bad-burst list '" << fBadBurstListFile
              << "', no bursts will be excluded on this basis." << std::endl;
    return;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    Int_t burst;
    if (iss >> burst) fBadBursts.insert(burst);
  }
  std::cout << "[MyPNNSelection] Loaded " << fBadBursts.size() << " excluded burst(s) from "
            << fBadBurstListFile << std::endl;
}

Bool_t MyPNNSelection::IsGoodBurst(Int_t burstID) {
  // Yes, I know this is not how you're supposed to do this. If you find
  // this, write me an email.
  static const std::set<Int_t> kTempExtraExclusions = {
      852, 884, 910, 926, 943, 964, 979, 992
  };
  if (kTempExtraExclusions.find(burstID) != kTempExtraExclusions.end()) return false;

  return fBadBursts.find(burstID) == fBadBursts.end();
}

void MyPNNSelection::StartOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyPNNSelection::Process(int) {
  if (!fReadingData) return;

  fPNNEventSelected = false;
  fPionTrackID = -1;
  fM2Miss = -9999.;
  fSignalRegion = 0;

  Int_t burstID = GetBurstID();
  fHEventsPerBurst->Fill(burstID);

  // ------------------------------------------------------------------
  // Burst-quality selection. Bursts affected by known detector or
  // beam-condition problems are excluded here, before any physics
  // selection is applied. The list is centrally maintained in
  // config/BadBurstList.dat -- see LoadBadBurstList().
  // ------------------------------------------------------------------
  if (!IsGoodBurst(burstID)) return;
  fHEventsPerBurstAfterQuality->Fill(burstID);

  Bool_t physicsTrig = TriggerConditions::GetInstance()->IsPhysicsTrigger(GetL0Data());
  Int_t L0TriggerWord = GetL0Data()->GetTriggerFlags();
  if (!(physicsTrig && (L0TriggerWord & fTriggerMask))) return;

  // ------------------------------------------------------------------
  // Single-track vertex selection
  // ------------------------------------------------------------------
  std::vector<SpectrometerTrackVertex> Vertices =
      *(std::vector<SpectrometerTrackVertex> *)GetOutput("SpectrometerVertexBuilder.Output1");

  Int_t bestVertex = -1;
  Double_t bestChi2 = 1e9;
  for (UInt_t iV = 0; iV < Vertices.size(); iV++) {
    Double_t Z = Vertices[iV].GetPosition().z();
    if (Z < fZVertexMin || Z > fZVertexMax) continue;

    Double_t CDA = Vertices[iV].GetCDA();
    fHCDA->Fill(CDA);
    if (CDA > fMaxCDA) continue;

    if (Vertices[iV].GetCharge() != 1) continue;

    Double_t p = Vertices[iV].GetTrackThreeMomentum(0).Mag();
    if (p < fMinTrackMomentum || p > fMaxTrackMomentum) continue;

    Double_t chi2 = Vertices[iV].GetChi2();
    if (chi2 > fMaxTrackChi2) continue;

    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestVertex = iV;
    }
  }
  if (bestVertex < 0) return;

  fHZVertex->Fill(Vertices[bestVertex].GetPosition().z());

  std::vector<DownstreamTrack> Tracks =
      *(std::vector<DownstreamTrack> *)GetOutput("DownstreamTrackBuilder.Output");
  fPionTrackID = Vertices[bestVertex].GetTrackIndex(0);

  // ------------------------------------------------------------------
  // Particle identification: reject e+/e- (calorimetric E/p) and
  // mu+/mu- (MUV3 association) to keep the pi+ hypothesis clean.
  // ------------------------------------------------------------------
  Double_t EoP = Tracks[fPionTrackID].GetLKrEoP();
  fHEoP->Fill(EoP);
  if (EoP > fMaxEoP) return;

  if (fRequireMUV3NoAssociation && Tracks[fPionTrackID].MUV3AssociationExists()) return;

  Double_t richLikelihoodPion = Tracks[fPionTrackID].GetRICHSinglePionProbability();
  if (richLikelihoodPion < fMinRICHLikelihoodPion) return;

  // ------------------------------------------------------------------
  // Kinematics: squared missing mass in the pi+ hypothesis
  // ------------------------------------------------------------------
  TVector3 BeamMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  TLorentzVector Kaon;
  Kaon.SetVectM(BeamMomentum, MKCH);
  TLorentzVector Pion;
  Pion.SetVectM(Vertices[bestVertex].GetTrackThreeMomentum(0), MPI);

  Double_t m2miss = (Kaon - Pion).M2() * 1e-6;  // MeV^2 -> GeV^2
  Double_t pPion = Pion.Vect().Mag() * 1e-3;    // MeV/c -> GeV/c

  fM2Miss = m2miss;
  fHM2Miss->Fill(m2miss);
  fHM2MissVsP->Fill(pPion, m2miss);

  if (m2miss > fR1Min && m2miss < fR1Max) {
    fSignalRegion = 1;
  } else if (m2miss > fR2Min && m2miss < fR2Max) {
    fSignalRegion = 2;
  } else {
    return;
  }

  fPNNEventSelected = true;
}

void MyPNNSelection::PostProcess() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyPNNSelection::EndOfBurstUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyPNNSelection::EndOfRunUser() {
  /// \MemberDescr
  /// \EndMemberDescr
}

void MyPNNSelection::EndOfJobUser() {
  SaveAllPlots();
}

void MyPNNSelection::DrawPlot() {
  /// \MemberDescr
  /// \EndMemberDescr
}

MyPNNSelection::~MyPNNSelection() {
  /// \MemberDescr
  /// \EndMemberDescr
}
