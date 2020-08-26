# DNPSoup

> A DNP Simulation Program

Simulate DNP Enhancements of FieldProfiles, Buildup Curves, and Parameter Scannings.


### Download

``` bash
git clone --recursive git@bitbucket.org:cyang019/dnpsoup.git
```


### Build

#### Ubuntu

``` bash
# install dependencies
apt-get update
apt-get install -y g++ git cmake ninja-build libopenblas-dev liblapacke-dev libpthread-stubs0-dev gfortran libatlas-base-dev

# can put this line in .bashrc
export OMP_NUM_THREADS=1

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
./build/dnpsoup_cli/dnpsoup_exec [input json file path] [output file path]
```

#### Examples
[Examples](./examples/inputs/)

