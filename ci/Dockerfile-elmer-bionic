FROM elmertest/bionic
RUN cmake -C ../opts-ubuntu-bionic.cmake -C ../opts-generic.cmake ../elmerfem
RUN make -j$(nproc)
