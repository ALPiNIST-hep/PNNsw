// ---------------------------------------------------------------------
// History:
//
// Created by Vic Marchetti (victor.marchetti@cern.ch) 2023-04-11
//
// ---------------------------------------------------------------------

/// \class MyPNNSelection
/// \Brief
/// Selection of K+ -> pi+ nu nu (PNN) candidates from the STRAW/RICH/LKr/MUV3/LAV
/// reconstruction, with associated burst-quality (good spill) selection.
/// \EndBrief
/// \Detailed
/// Signal candidates are selected in two disjoint regions of the squared
/// missing mass m2miss = (P_K - P_pi)^2, referred to as R1 (low mass, near
/// the pi+ pi0 exclusion boundary) and R2 (high mass, above the pi+ pi0
/// upper edge). Background rejection relies on the upstream veto counters,
/// full photon-veto coverage and RICH/calorimetric particle identification.
/// Only bursts flagged as "good" by the standard NA62 burst-quality
/// selection (see IsGoodBurst()) are retained.
/// \EndDetailed

#ifndef MYPNNSELECTION_HH
#define MYPNNSELECTION_HH

#include "Analyzer.hh"
#include "MCSimple.hh"
#include "TH1F.h"
#include "TH2F.h"

#include <map>
#include <set>
#include <string>
#include <vector>

class TRecoSpectrometerEvent;
class TRecoCedarEvent;
class TRecoGigaTrackerEvent;

class MyPNNSelection : public NA62Analysis::Analyzer {

 public:
  explicit MyPNNSelection(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyPNNSelection();
  void InitOutput();
  void InitHist();
  void DefineMCSimple() {}
  void Process(int iEvent);
  void PostProcess();
  void StartOfBurstUser();
  void EndOfBurstUser();
  void EndOfRunUser();
  void EndOfJobUser();
  void DrawPlot();

 private:
  // ---- burst quality ----------------------------------------------------
  void LoadBadBurstList();
  Bool_t IsGoodBurst(Int_t burstID);

  std::string fBadBurstListFile;   ///< path of the officially maintained bad-burst file
  std::set<Int_t> fBadBursts;      ///< bursts excluded from the analysis

  // ---- cut values (see constructor for defaults) -------------------------
  UInt_t fTriggerMask;

  Double_t fMinTrackMomentum;
  Double_t fMaxTrackMomentum;
  Double_t fMaxTrackChi2;

  Double_t fZVertexMin;
  Double_t fZVertexMax;
  Double_t fMaxCDA;

  Double_t fMinRICHLikelihoodPion;
  Double_t fMaxEoP;                 ///< upper cut on E/p (LKr/MUV1,2) to reject e+/e-
  Bool_t fRequireMUV3NoAssociation;  ///< reject events with a MUV3 in-time hit (mu veto)

  Double_t fR1Min, fR1Max;   ///< signal region R1 in m2miss [GeV^2/c^4]
  Double_t fR2Min, fR2Max;   ///< signal region R2 in m2miss [GeV^2/c^4]

  // ---- outputs ------------------------------------------------------------
  Bool_t fReadingData;
  Bool_t fPNNEventSelected;
  Int_t fPionTrackID;
  Double_t fM2Miss;
  Int_t fSignalRegion;    ///< 0 = none, 1 = R1, 2 = R2

  // ---- histograms -----------------------------------------------------
  TH1F *fHM2Miss;
  TH2F *fHM2MissVsP;
  TH1F *fHZVertex;
  TH1F *fHCDA;
  TH1F *fHEoP;
  TH1F *fHEventsPerBurst;
  TH1F *fHEventsPerBurstAfterQuality;
};

#endif
