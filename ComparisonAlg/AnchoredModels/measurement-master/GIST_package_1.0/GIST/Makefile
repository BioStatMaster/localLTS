MATLABDIR ?=/home/mgong/MATLAB/R2014a
CXX = /usr/bin/g++ #compile in Macbook
CFLAGS = -Wall -O0 -fPIC -I$(MATLABDIR)/extern/include -I/home/mgong/libs/boost/include
MEX = $(MATLABDIR)/bin/mex
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)
MEX_OPTION +=  -cxx

all: gradientA.$(MEX_EXT)

gradientA.$(MEX_EXT): gradientA.o
		$(MEX) $(MEX_OPTION) -output $@ $^
			
gradientA.o: gradientA.cpp
		$(CXX) $(CFLAGS) -o $@ -c $^	

clean: 
		rm -f *.mex*
			rm -f *.o		 

