LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so CPUPROFILE=main.prof CPUPROFILE_FREQUENCY=100000 @
google-pprof --pdf ./svzerodsolver main.prof > profiling_report.pdf