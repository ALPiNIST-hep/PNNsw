# PNNsw

Event selection and background-estimation software for the NA62
K+ -> pi+ nu nu (PNN) branching-ratio analysis, built on top of the
standard na62fw `Analyzer` framework.

## Contents

### Signal selection

- `include/MyPNNSelection.hh`, `src/MyPNNSelection.cc` -- main selection
  analyzer. Reads the `Cedar` and `GigaTracker` reco trees plus the standard
  downstream-track/vertex builder outputs, applies track quality, vertex and
  particle-identification cuts, and classifies candidates into the two
  squared-missing-mass signal regions R1 and R2.
- `config/BadBurstList.dat` -- centrally maintained burst-quality (good
  spill) exclusion list used by `MyPNNSelection::IsGoodBurst()`.

### Background estimation

Each channel is estimated by a dedicated analyzer, following the same
control-region / tail-fraction (or acceptance-loss) approach:

- `MyK2piBackground` -- K+ -> pi+ pi0(gamma), from the radiative tail of
  the missing mass reconstructed in the pi0 hypothesis.
- `MyKmu2Background` -- K+ -> mu+ nu(gamma), from mu+ -> pi+
  mis-identification, measured with a MUV3-tagged tag-and-probe sample.
- `MyK3piBackground` -- K+ -> pi+ pi+ pi-, from partial reconstruction
  when two of the three tracks are lost.
- `MyKe4Background` -- K+ -> pi+ pi- e+ nu, from partial reconstruction
  of the electron + pi- pair.
- `MyUpstreamBackground` -- pi+ tracks from beam-material interactions or
  accidentally matched upstream K+ decays, via the inverted-CDA Upstream
  Reference Sample (URS).
- `MyPNNBackgroundSummary` -- runs after the above and after
  `MyPNNSelection` in the same job; combines their outputs into the final
  per-region background table (`background_table.txt` by default).

## Usage

Add `MyPNNSelection`, the five background analyzers and
`MyPNNBackgroundSummary` (in that order) to your `na62fw` analyzer list,
pointing `BadBurstListFile` at `config/BadBurstList.dat` (default). See
the class documentation in each header for the full list of configurable
cuts and normalisation constants.
