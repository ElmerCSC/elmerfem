FROM ubuntu:xenial
RUN apt-get -y update && apt-get install -y gcc cmake g++ gfortran libopenblas-dev
RUN apt-get install -y unzip libopenmpi-dev
RUN apt-get install -y libmumps-dev libparmetis-dev libmetis-dev
ADD ./ ./elmerfem
ADD ci/precaches/opts-ubuntu-xenial.cmake /opts-ubuntu-xenial.cmake
ADD ci/precaches/opts-generic.cmake /opts-generic.cmake
RUN ls
RUN mkdir elmerbuild
WORKDIR elmerbuild
