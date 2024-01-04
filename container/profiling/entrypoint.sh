#! /bin/bash

export LIBPROFILER=$(find /usr/lib -name "libprofiler.so")
LD_PRELOAD=$LIBPROFILER CPUPROFILE=/opt/main.prof CPUPROFILE_FREQUENCY=100000 /opt/pysvzerod-build/relwithdebinfo/svzerodsolver $@ /opt/output.csv
google-pprof --pdf /opt/pysvzerod-build/relwithdebinfo/svzerodsolver /opt/main.prof > profiling_report.pdf
