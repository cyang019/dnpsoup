FROM ubuntu:latest

WORKDIR /usr/src/dnpsoup

RUN apt-get update \
    && ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y tzdata \
    && dpkg-reconfigure --frontend noninteractive tzdata \
    && apt-get install -y g++ git cmake ninja-build libopenblas-dev liblapacke-dev libpthread-stubs0-dev gfortran

COPY ./cmake ./cmake
COPY ./configure_dnpsoup.h.in ./
COPY ./CMakeLists.txt ./
COPY ./matrix/matrix_impl ./matrix/matrix_impl
COPY ./tests ./tests
COPY ./dnpsoup_cli ./dnpsoup_cli
COPY ./dnpsoup_impl ./dnpsoup_impl

RUN mkdir build \
    && cd build \
    && cmake .. -GNinja \
    && ninja

CMD ["./build/dnpsoup_cli/dnpsoup_exec", "mounted/input.json", "mounted/output.result"]
