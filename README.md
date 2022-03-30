# DNPSOUP

[![CircleCI](https://circleci.com/gh/cyang019/dnpsoup.svg?style=svg)](https://circleci.com/gh/cyang019/dnpsoup)


> **D**ynamic **N**uclear **P**olarization **S**imulation **O**ptimized with a **U**nified **P**ropagator

DNPSOUP let users simulate field profiles **without preparing a single equation**. It gives insights about polarization enhancement on nuclei(s) due to the presence of electron(s) and EM radiations, with or without Magic Angle Spinning (MAS).

DOI: [10.1016/j.jmr.2021.107107](https://doi.org/10.1016/j.jmr.2021.107107)

please email `yangcnju@gmail` for questions or suggestions.

## Graphical User Interface for Input Preparation

**dnpsoup_gui** is a graphical user interface that helps generating simulation inputs:

[https://github.com/cyang019/dnpsoup_gui](https://github.com/cyang019/dnpsoup_gui)

--------
## Table of Contents
- [Background Knowledge](#background-knowledge)
- [Getting Started (Linux / MacOS)](#getting-started-for-linux-and-macos)
- [Getting Started (Windows)](#getting-started-for-windows)
- [Expected Input JSON File Format](#expected-input-json-file-format)
- [Authors](#authors)
- [Citing DNPSOUP](#citing-dnpsoup)
- [Acknowledgements](#acknowledgements)
- [External Links](#links)
- [License](#license)


----
## Background Knowledge

### Nuclear Magnetic Resonance (NMR)

NMR is an analytical chemistry technique that utilizes a high field magnet to measure molecular or atomic structures.
  - [What is NMR (http://chem.ch.huji.ac.il/nmr/whatisnmr/whatisnmr.html)](http://chem.ch.huji.ac.il/nmr/whatisnmr/whatisnmr.html)
  - [Solid-state NMR spectroscopy (https://www.nature.com/articles/s43586-020-00002-1)](https://www.nature.com/articles/s43586-020-00002-1)
  - [Keeler notes (http://www-keeler.ch.cam.ac.uk/lectures/)](http://www-keeler.ch.cam.ac.uk/lectures/)


### Dynamic Nuclear Polarization (DNP)

DNP relies on the transfer of electron polarization (typically from an organic based exogenous radical) to neighboring nuclei offering significant gains in sensitivity.

  - [High-field DNP (https://griffingroup.mit.edu/high-field-dnp)](https://griffingroup.mit.edu/high-field-dnp#overlay-context=research/nmr-methodology)


## Getting Started for Linux and MacOS

### Download

``` bash
git clone --recursive https://github.com/cyang019/dnpsoup.git
```

### Prerequisites

- c++ compiler that supports c++17 features
- cmake version 3.14.7 and above

### Install Dependencies

#### Ubuntu

``` bash
sudo apt-get update
sudo apt-get install -y g++ git cmake ninja-build libopenblas-dev liblapacke-dev libpthread-stubs0-dev gfortran libatlas-base-dev
```

#### Manjaro Linux

```bash
sudo pacman -Syu clang cmake ninja git openblas lapacke
```

#### MacOS

Install xcode and homebrew

```bash
# homebrew
brew install cmake ninja git
```

### Build

```bash
# can put this line in .bashrc
export OMP_NUM_THREADS=1

# in the root directory of dnpsoup
mkdir build
cd build
cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release
ninja

# test if build was succesful
./tests/test_dnpsoup_core
```

### Usage

To execute dnpsoup:
``` bash
# in the root directory of dnpsoup
./build/dnpsoup_cli/dnpsoup_exec [input json file path] [output file path]
```

## Getting Started for Windows

### Prerequisites

+ git: [https://git-scm.com/download/win](https://git-scm.com/download/win)

+ Visual Studio: [https://visualstudio.microsoft.com](https://visualstudio.microsoft.com)
  
  Install c++ development environment by using *Visual Studio Installer*
  + install c++ desktop development tools
  + install cmake

+ mkl (BLAS/LAPACK support): [https://software.intel.com/oneapi/onemkl](https://software.intel.com/oneapi/onemkl)
  - make sure *Intel Math Kernel Library* within the Software Installer is selected. Install mkl in the default location within *C:\Program Files (x86)*.

### Download and Compile

1. Open `Git Bash`, browse to a desired folder,
```bash
git clone --recursive https://github.com/cyang019/dnpsoup.git
```

2. Open `Visual Studio`, on the menu bar, `File`->`Open`->`Folder...` to open the dnpsoup folder. Within the `Solution Explorer - Folder View`, right click on the **CMakeLists.txt** in the root directory of the project to build.

3. Open the `Developer Command Prompt for VS <version number>`, browse to `<dnpsoup_folder>\out\build\x64-Release or x64-Debug\out`, the executable `dnpsoup_exec.exe` should locate here.

## Examples

+ [Solid Effect and Overhauser Effect](./examples/01_solid_effect_and_overhauser_effect)

+ [Cross Effect](./examples/02_cross_effect/ce_visualization.ipynb)

+ [NOVEL](./examples/03_novel/NOVEL_visualization.ipynb)

+ [Offres-NOVEL](./examples/04_off_resonance_novel/OffresNOVEL_visualization.ipynb)

+ [ISE SSE](./examples/05_ise_sse/ISE_SSE_visualization.ipynb)

+ [TOP-DNP](./examples/06_top_dnp/top_dnp_visualization.ipynb)


## Input JSON File Example

```json
{
  "spinsys": {
    "euler": {
      "alpha": 0.0,
      "beta": 0.0,
      "gamma": 0.0
    },
    "interactions": [
      {
        "entries": {
        },
        "name": "csa/shielding/hyperfine/dipole/scalar"
      }
    ],
    "spins": {
      "spin-id": {
        "type": "e/H1/C13 etc",
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
        "t1": 1.0e-3,
        "t2": 1.0e-6
      }
    }
  },
  "pulseseq": {
    "name": "name of the pulse sequence",
    "components": {
      "emr-name": {
        "channel-name": {
          "frequency": 1.0e6,
          "offset": 0.0,
          "phase": 0.0
        }
      }
    },
    "increment": 1.0e-9,
    "sections": {
      "section-name": {
        "names-within-section": [
          "name1"          
        ],
        "params": {},
        "size": 100,
        "type": "section-type"
      }
    },
    "sequence": [
      "section1", "section2", "section3"
    ]
  },
  "settings": {
    "euler": {
      "alpha": 0.0,
      "beta": 0.0,
      "gamma": 0.0
    },
    "ncores": 1,
    "acq": "H1",
    "Magnet": {
      "b0": 9.4
    },
    "Gyroton": {
      "em_frequency": 263.0e9
    },
    "Probe": {
      "mas_frequency": 0.0,
      "temperature": 77.0,
      "mas_increment": 1e-5
    },
    "task": "FieldProfile",

    "other-task-related-entries"
  }
}
```


## Authors
Chen Yang, Kong Ooi Tan, and Robert G. Griffin

## Citing DNPSOUP
Yang, C., Tan, K. O., & Griffin, R. G. (2021). DNPSOUP: A simulation software package for dynamic nuclear polarization. Journal of Magnetic Resonance, 107107.
DOI: [10.1016/j.jmr.2021.107107](https://doi.org/10.1016/j.jmr.2021.107107)


## Acknowledgements

- The authors would like to thank Prof. Sheetal Kumar Jain for initial discussions and encouragement on creating a separate simulation package.
- The authors would also like to thank Dr. Kan-Nian Hu for discussions on cross effect, and Dr. John Wright for computational resources of PSFC at MIT. 
- Chen Yang thanks Professor Leonard Mueller for introducing him to solid-state NMR and product operators. 
- Chen Yang also thanks his family for support, and Mr. Xiang Gao for his suggestions on debugging frontend javascript code.

## Fundings

+ High Field DNP and EPR in Biological Systems:
  [https://grantome.com/grant/NIH/R01-GM132997-32](https://grantome.com/grant/NIH/R01-GM132997-32)

+ MIT/Harvard Center for Magnetic Resonance:
  [https://grantome.com/grant/NIH/P41-GM132079-01](https://grantome.com/grant/NIH/P41-GM132079-01)

+ Advancing NMR Sensitivity Beyond Continuous-Wave DNP - PulsedDNP: [https://anr.fr/Project-ANR-20-ERC9-0008](https://anr.fr/Project-ANR-20-ERC9-0008)


## Links

+ Griffin Group website: [https://griffingroup.mit.edu](https://griffingroup.mit.edu)
+ Kong Ooi Tan's Webpage: [https://www.chimie.ens.fr/tan](https://www.chimie.ens.fr/tan)


## License

[BSD-3-Clause License](./LICENSE)

