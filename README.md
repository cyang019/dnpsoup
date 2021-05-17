# DNPSOUP

[![CircleCI](https://circleci.com/gh/cyang019/dnpsoup.svg?style=svg&circle-token=329696edcd1550614500348cdd674eb9900ee51b)](https://circleci.com/gh/cyang019/dnpsoup)


> **D**ynamic **N**uclear **P**olarization **S**imulation **O**ptimized with a **U**nified **P**ropagator

DNPSOUP simulates polarization enhancement on nuclei(s) due to the presence of electron(s) and EM radiations, with or without Magic Angle Spinning (MAS).

## Background Knowledge

#### Nuclear Magnetic Resonance (NMR)

NMR is an analytical chemistry technique that utilizes a high field magnet to measure molecular or atomic structures.
  - [What is NMR (http://chem.ch.huji.ac.il/nmr/whatisnmr/whatisnmr.html)](http://chem.ch.huji.ac.il/nmr/whatisnmr/whatisnmr.html)
  - [Solid-state NMR spectroscopy (https://www.nature.com/articles/s43586-020-00002-1)](https://www.nature.com/articles/s43586-020-00002-1)
  - [Keeler notes (http://www-keeler.ch.cam.ac.uk/lectures/)](http://www-keeler.ch.cam.ac.uk/lectures/)


#### Dynamic Nuclear Polarization (DNP)

DNP relies on the transfer of electron polarization (typically from an organic based exogenous radical) to neighboring nuclei offering significant gains in sensitivity.

  - [High-field DNP (https://griffingroup.mit.edu/high-field-dnp)](https://griffingroup.mit.edu/high-field-dnp#overlay-context=research/nmr-methodology)

## Getting Started

### Download

``` bash
git clone --recursive https://github.com/cyang019/dnpsoup.git
```

### Dependencies

- nlohmann/json ([https://github.com/nlohmann/json](https://github.com/nlohmann/json))
- openblas & lapack / accelerate

#### Install Dependencies

##### Ubuntu

``` bash
sudo apt-get update
sudo apt-get install -y g++ git cmake ninja-build libopenblas-dev liblapacke-dev libpthread-stubs0-dev gfortran libatlas-base-dev
```

##### Manjaro Linux

```bash
sudo pacman -Syu clang cmake ninja git openblas lapacke
```

##### MacOS

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

#### To execute dnpsoup

``` bash
# in the root directory of dnpsoup
./build/dnpsoup_cli/dnpsoup_exec [input json file path] [output file path]
```

#### Examples

+ [Solid Effect and Overhauser Effect](./examples/01_solid_effect_and_overhauser_effect/se_oe_visualization.ipynb)

+ [Cross Effect](./examples/02_cross_effect/ce_visualization.ipynb)

+ [NOVEL](./examples/03_novel/NOVEL_visualization.ipynb)

+ [Offres-NOVEL](./examples/04_off_resonance_novel/OffresNOVEL_visualization.ipynb)

+ [ISE SSE](./examples/05_ise_sse/ISE_SSE_visualization.ipynb)

+ [TOP-DNP](./examples/06_top_dnp/top_dnp_visualization.ipynb)


### Graphical Interface for Input Preparation

**dnpsoup_gui** is a graphical user interface that helps generating simulation inputs:
[https://github.com/cyang019/dnpsoup_gui](https://github.com/cyang019/dnpsoup_gui)


#### Expected Input JSON File Format

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
*DNPSOUP: A numerical simulation software for dynamic nuclear polarization*, in preparation

## Fundings

+ High Field DNP and EPR in Biological Systems:
  [https://grantome.com/grant/NIH/R01-GM132997-32](https://grantome.com/grant/NIH/R01-GM132997-32)

+ MIT/Harvard Center for Magnetic Resonance:
  [https://grantome.com/grant/NIH/P41-GM132079-01](https://grantome.com/grant/NIH/P41-GM132079-01)


## Acknowledgements
- The authors would like to thank Dr. Sheetal Kumar Jain for initial discussions and encouragement on creating a separate simulation package.
- The authors would also like to thank Dr. Kan-Nian Hu for discussions on cross effect, and Dr. John Wright for computational resources of PSFC at MIT. 
- Chen Yang thanks Professor Leonard Mueller for introducing him to solid-state NMR and product operators. 
- Chen Yang also thanks his family for support, and Mr. Xiang Gao for his suggestions on debugging frontend javascript code.


## Links

+ Griffin Group website: [https://griffingroup.mit.edu](https://griffingroup.mit.edu)
+ Kong Ooi Tan's Webpage: [https://www.chimie.ens.fr/tan](https://www.chimie.ens.fr/tan)


## License

[BSD-3-Clause License](./LICENSE)

