# build an executable named mysampling
all:samplingmain.C samplingfunctions.C 
	g++ -c samplingfunctions.C -o samplingfunctions.o -g -Wall -Wextra -pedantic
	g++ -c samplingmain.C -o samplingmain.o -g -Wall -Wextra -pedantic
	g++ samplingfunctions.o   samplingmain.o -o mysampling -O3 -g -Wall -Wextra -pedantic