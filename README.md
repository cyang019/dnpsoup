# DNPSOUP

[![pipeline status](https://gitlab.com/chen.yang/dnpsoup/badges/master/pipeline.svg)](https://gitlab.com/chen.yang/dnpsoup/-/commits/master)

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
git clone --recursive git@gitlab.com:chen.yang/dnpsoup.git
```

### Dependencies

- git
- cmake
- ninja
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

Install xcode
Install homebrew

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


**[dnpsoup_gui](https://github.com/cyang019/dnpsoup_gui)** can help to generate simulation inputs.


#### Expected Input JSON File Format

```json
{
  "spinsys": {
    "euler": {
      "alpha": <float>,
      "beta": <float>,
      "gamma": <float>
    },
    "interactions": [
      ...
    ],
    "spins": {
      ...
    }
  },
  "pulseseq": {
    "name": "name of the pulse sequence",
    "components": {
      ...
    },
    "sections": {
      ...
    },
    "sequence": [
      ...
    ]
  },
  "settings": {
    "euler": {
      "alpha": <float>,
      "beta": <float>,
      "gamma": <float>
    },
    "ncores": <int>,
    "acq": "H1",
    "Magnet": {
      "b0": <float>
    },
    "Gyroton": {
      "em_frequency": <float>
    },
    "Probe": {
      "mas_frequency": <float>,
      "temperature": <float>,
      "mas_increment": <float>
    },
    "task": <task name>,
    ...
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


## License

[BSD-3 License](./LICENSE)

