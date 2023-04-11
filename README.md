# PNNsw

Event selection software for the NA62 K+ -> pi+ nu nu (PNN) branching-ratio
analysis, built on top of the standard na62fw `Analyzer` framework.

## Contents

- `include/MyPNNSelection.hh`, `src/MyPNNSelection.cc` -- main selection
  analyzer. Reads the `Cedar` and `GigaTracker` reco trees plus the standard
  downstream-track/vertex builder outputs, applies track quality, vertex and
  particle-identification cuts, and classifies candidates into the two
  squared-missing-mass signal regions R1 and R2.
- `config/BadBurstList.dat` -- centrally maintained burst-quality (good
  spill) exclusion list used by `MyPNNSelection::IsGoodBurst()`.

## Usage

Add `MyPNNSelection` to your `na62fw` analyzer list and point
`BadBurstListFile` at `config/BadBurstList.dat` (default). See the class
documentation in `include/MyPNNSelection.hh` for the full list of
configurable cuts.
