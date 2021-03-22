# dnpsoup

> A Dynamic Nuclear Polarization Simulation Program

Simulate polarization enhancement on nuclei(s) due to the presence of electron(s) and electromagnetic radiations, with or without Magic Angle Spinning (MAS).


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

[Sample Input Scripts](./examples/inputs/)

**[dnpsoup_analytics](https://gitlab.com/chen.yang/dnpsoup_analytics)** contains details of simulation results for some of the above simulation examples.

One can use **[dnpsoup_gui](https://github.com/cyang019/dnpsoup_gui)** to generate simulation inputs for convenience.


### License

[BSD-3](./LICENSE)

