# Pion Proton Elastic Scattering Event Generator (PiPES)
- Purpose: low energy pion-proton elastic scattering at both charge polarties 
- Current capabilities:
    - Cap partial waves at p-wave; can be extended without much difficulty
    - Limit lab frame momentum options to 115, 160 and 210 MeV/c; more can be added 
    - $\pi^{+}p$ works well; $\pi^{-}p$ needs small corrections

On set of 20 iterations sampling 1000000 events, the generator ran at a speed of $902\pm209$ ms, corresponding to an average throughput of $1.11$ million events/s.

## Overview
PiPES is a pure C++, no external packages required event generator. 