---------------------------------------------------------------
           Compile and install ElmerPost's plugins:
---------------------------------------------------------------

$ export CFLAGS="-I/usr/include/tcl8.4 -I/usr/include/ffmpeg"
$ ./configure --prefix=$ELMER_POST_HOME
$ make
$ make install

Usually ELMER_POST_HOME is defined as

$ export ELMER_POST_HOME=$ELMER_HOME/share/elmerpost

If the configuration script fails to determine correct
setting, compilation intructions can be found from the
headers of the individual files (*.c).
