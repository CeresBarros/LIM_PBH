# LIM_PBH

Simulation experiment on projecting **b**oreal forest–fire dynamics in Alberta's
foothills (the *Landscapes in Motion* — **P**rojecting **B**oreal from **H**istorical
component), built on the [SpaDES](https://spades.predictiveecology.org/) modelling
framework and the LandR / fireSense module ecosystem.

This repository contains the simulation driver scripts, module configuration, and
post-simulation analysis code for the SpaDES side of the project. Empirical
fire-severity modelling for the same study system lives in a separate repository:
[`LIM_FireSevModels`](https://github.com/CeresBarros/LIM_FireSevModels).

## Repository layout

```
R/
├── SpaDES/            # simulation drivers, module instantiation, post-sim analyses
│   ├── 0_packages.R
│   ├── 1_simObjects.R
│   ├── 2_speciesLayers.R
│   ├── 3_fireWeather.R
│   ├── 4_preSimulation.R
│   ├── 5_modelDiagnosticsR.R
│   ├── 7_analysesResultsMontane.Rmd
│   ├── 8_hypervolumes.R
│   ├── 9_hypervolumesAnalyses_*.R
│   ├── global.R
│   └── m/             # SpaDES modules (git submodules — see below)
├── R_tools/           # helper functions used across scripts
└── Favier_FM.R        # Favier fire-model exploration
data/                  # small reference data + CHECKSUMS; larger inputs are gitignored
Docker/                # Dockerfile + run scripts
packages/              # library path for host R (contents gitignored)
packages_docker/       # library path for Docker R (contents gitignored)
LIM_PBH.Rproj
```

## Cloning (submodules)

The SpaDES modules under `R/SpaDES/m/` are tracked as **git submodules** with SSH
URLs. Clone recursively:

```bash
git clone --recurse-submodules git@github.com:CeresBarros/LIM_PBH.git
```

Or, if you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

## R environment

`.Rprofile` (not tracked; provide your own) sets a project-local library path under
`packages/<platform>/<R version>/` (or `packages_docker/…` on the Docker image), and
adds the `predictiveecology.r-universe.dev` repo. See the original template in the
`LandscapesInMotion` repo, or reproduce with:

```r
options(repos = c(
  CRAN = "https://cran.rstudio.com",
  PE   = "https://predictiveecology.r-universe.dev/"
))
```

R version 4.5 is the current development target; the Docker image (see `Docker/`)
provides a reproducible environment.

## Getting started

1. Clone with submodules (above).
2. Provide a local `.Rprofile` / `.Renviron` if you need custom paths or credentials.
3. Open `LIM_PBH.Rproj` in Positron or RStudio.
4. Source `R/SpaDES/0_packages.R` to install dependencies.
5. Work through `R/SpaDES/global.R` (or `global.Rmd`) to launch a simulation.

## Provenance

This repository was derived from the `development` branch of the original
[`LandscapesInMotion`](https://github.com/CeresBarros/LandscapesInMotion) repository
in July 2026, preserving history for the SpaDES-related paths.

## Licence

Released under the **GNU General Public License v3.0 or later** (GPL-3.0-or-later),
consistent with the wider SpaDES / PredictiveEcology ecosystem. Crown copyright
applies — see [`LICENSE`](LICENSE) for the notice and full licence text.
