
CC = g++ -Wall -ggdb
CC = g++ -pg
CC = g++

# without OpenMP
CCFLAGS = -O3 -DNO_OPENMP

# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),yes)
CCFLAGS = -O3 -fopenmp
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb
endif

#LDFLAGS = -static -o
LDFLAGS = -o

PROGS = cd-hit cd-hit-est cd-hit-2d cd-hit-est-2d cd-hit-div

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)

clean:
	rm *.o $(PROGS)

# programs

cd-hit: cdhit-common.o cdhit-utility.o cdhit.o
	$(CC) $(CCFLAGS) cdhit.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit

cd-hit-2d: cdhit-common.o cdhit-utility.o cdhit-2d.o
	$(CC) $(CCFLAGS) cdhit-2d.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-2d

cd-hit-est: cdhit-common.o cdhit-utility.o cdhit-est.o
	$(CC) $(CCFLAGS) cdhit-est.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-est

cd-hit-est-2d: cdhit-common.o cdhit-utility.o cdhit-est-2d.o
	$(CC) $(CCFLAGS) cdhit-est-2d.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-est-2d

cd-hit-div: cdhit-common.o cdhit-utility.o cdhit-div.o
	$(CC) $(CCFLAGS) cdhit-div.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-div

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-utility.c++ -c

cdhit.o: cdhit.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit.c++ -c

cdhit-2d.o: cdhit-2d.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-2d.c++ -c

cdhit-est.o: cdhit-est.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-est.c++ -c

cdhit-est-2d.o: cdhit-est-2d.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-est-2d.c++ -c

cdhit-div.o: cdhit-div.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-div.c++ -c

