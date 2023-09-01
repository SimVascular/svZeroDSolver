
src=LPNSolverInterface.cpp main.cpp
obj=$(src:.cpp=.o)

CC = g++

copt = -g
copt = -O
#flags =  -std=c++11 -stdlib=libc++ -w $(inc) $(copt)
flags =  -std=c++11  -w $(inc) $(copt)

svZeroD_interface_test: $(obj) LPNSolverInterface.h
	$(CC) -o $@ $(obj) $(libs) -ldl

.cpp.o: 
	$(CC) $(flags) -c $<

clean:
	rm -f *.o svZeroD_interface_test

