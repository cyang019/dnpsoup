FROM ubuntu:latest

WORKDIR /usr/src/dnpsoup

COPY CMakeLists.txt .
COPY cmake .
COPY matrix .
COPY configure_dnpsoup.h.in .
COPY tests .
COPY dnpsoup_cli .
COPY dnpsoup_impl .

RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN export DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata
RUN apt-get install -y g++ git cmake ninja-build libopenblas-dev liblapacke-dev libpthread-stubs0-dev gfortran
RUN mkdir build
RUN cd build
RUN cmake .. -GNinja
RUN ninja
RUN cd ..

CMD ["./build/dnpsoup_cli/dnpsoup_exec", "input.json", "output.result"]
