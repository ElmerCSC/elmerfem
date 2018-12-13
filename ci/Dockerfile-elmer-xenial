FROM elmertest/xenial
RUN cmake -C ../opts-ubuntu-xenial.cmake -C ../opts-generic.cmake ../elmerfem
RUN make -j$(nproc)
