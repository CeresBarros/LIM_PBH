# LIM_PBH

Simulation experiment on projecting historic forest–fire dynamics in Alberta's
SW Rockies foothills under different fire regime assumptions, to explore support for pyrodiversity and biodiversity hypothesis in this landscape (PBH). Part of the *Landscapes in Motion* project (Foothills Research Institute) which sought at understanding and projecting historic fire regimes in this landscape.

Simulation framework uses LandR and associated fire modelling components, a model built on the [SpaDES](https://spades.predictiveecology.org/) modelling toolkit.

This repository contains the simulation driver and configuration scripts and
post-simulation analysis code. 

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
│   ├── 6_analysesResultsMontane.Rmd
│   ├── 7_hypervolumes.R
│   ├── 8_hypervolumesAnalyses.R
│   ├── 9_generalFigures.Rmd
│   ├── global.R
│   └── m/             # SpaDES modules (git submodules — see below)
├── R_tools/           # helper functions used across scripts
data/                  # small reference data + CHECKSUMS; larger inputs are gitignored
Docker/                # Dockerfile + run scripts
packages/              # library path for host R (contents gitignored) -- not used/necessary except for tests
packages_docker/       # library path for Docker R (contents gitignored)
LIM_PBH.Rproj
```

## Simulation and analyses environment

All simulations and analyses were run using R **4.1.3** in a Docker container running in Ubuntu. 
Running simulations and analyses on a newer R version is not supported; Windows or Mac OSs have not been tested.

The instructions below assume a Ubuntu OS with `git` and `docker` installed.

### Clone repository locally with submodules

All simulation SpaDES modules under `R/SpaDES/m/` are tracked as **git submodules** with SSH
URLs; they wil need to be cloned recursively within the main repository:

```bash
git clone --recurse-submodules git@github.com:CeresBarros/LIM_PBH.git
```

Or, if you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

### Run Docker container

Pull Docker image

```bash
docker pull ceresbarros/lim_pbh
```

Launch the container by executing the  helper script `Docker/dockerRun_generic.sh`. Make sure your user has execute permission on the file (`chmod 700` line will grant these permissions).

```bash
chmod 700 Docker/dockerRun_generic.sh

./Docker/dockerRun_generic.sh
```

### Continue on RStudio Server

Then open RStudio Server in the browser at http://localhost:8787 and log
in with the user/password set in `dockerRun_generic.sh`.

`.Rprofile` (not tracked; provide your own) should set a project-local library path under
`packages_docker/…` and add the `predictiveecology.r-universe.dev` repo:

```r
options(repos = c(
  CRAN = "https://cran.rstudio.com",
  PE   = "https://predictiveecology.r-universe.dev/"
))
```

`.Renviron` should set your personal GitHub PAT if you have/need one.

### Simulations/analyses

1. Open `LIM_PBH.Rproj` in Positron or RStudio.

2. Work through `R/SpaDES/global.R` (or `global.Rmd`) to set up and launch simulations.

3. After simulations complete, run the post-simulation analyses in
   `R/SpaDES/` in numeric order — i.e. `6_analysesResultsMontane.Rmd`, `7_hypervolumes.R`, and
   `8_hypervolumesAnalyses.R`.

4. Other scripts provide accessory analyses (`5_modelDiagnosticsR.R`) or figures (`9_generalFigures.Rmd`)

## Provenance

This repository was derived from the `LIMpub_jul2026` tag of the original
[`LandscapesInMotion`](https://github.com/CeresBarros/LandscapesInMotion) repository
in July 2026, preserving history for the simulation-related paths.

## Licence

Released under the **GNU General Public License v3.0 or later** (GPL-3.0-or-later),
consistent with the wider SpaDES / PredictiveEcology ecosystem. Crown copyright
applies — see [`LICENSE`](LICENSE) for the notice and full licence text.

## Disclaimer

Portions of this README and repository documentation were drafted with the
assistance of Anthropic's Claude. All AI-assisted output was reviewed and edited by the
authors, who remain responsible for the contents of this repository.