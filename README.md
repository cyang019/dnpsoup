# DNPSOUP

[![pipeline status](https://gitlab.com/chen.yang/dnpsoup/badges/master/pipeline.svg)](https://gitlab.com/chen.yang/dnpsoup/-/commits/master)

> **D**ynamic **N**uclear **P**olarization **S**imulation **O**ptimized with a **U**nified **P**ropagator

Simulates polarization enhancement on nuclei(s) due to the presence of electron(s) and EM radiations, with or without Magic Angle Spinning (MAS).

## Getting Started

### Download

``` bash
git clone --recursive git@gitlab.com:chen.yang/dnpsoup.git
```


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

+ [NOVEL]()

+ [Offres-NOVEL]()

+ [ISE SSE]()

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
Chen Yang, Kong Ooi Tan, Robert G. Griffin


## Acknowledgements
- The authors would like to thank Dr. Sheetal Kumar Jain for initial discussions and encouragement on creating a separate simulation package.
- The authors would also like to thank Dr. Kan-Nian Hu for discussions on cross effect, and Dr. John Wright for computational resources of PSFC at MIT. 
- Chen Yang thanks Professor Leonard Mueller for introducing him to solid-state NMR and product operators. 
- Chen Yang also thanks his family for support, and Mr. Xiang Gao for his suggestions on debugging frontend javascript code.


## License

[BSD-3](./LICENSE)

