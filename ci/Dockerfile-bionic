FROM ubuntu:bionic
RUN apt-get -y update && apt-get install -y gcc cmake g++ gfortran libopenblas-dev
RUN apt-get install -y unzip libopenmpi-dev
RUN apt-get install -y libmumps-dev libparmetis-dev openssh-client
ADD ./ ./elmerfem
ADD ci/precaches/opts-ubuntu-bionic.cmake /opts-ubuntu-bionic.cmake
ADD ci/precaches/opts-generic.cmake /opts-generic.cmake
RUN ls
RUN mkdir elmerbuild
WORKDIR elmerbuild
