########################################################################
# Compiler and external dependences
########################################################################
CC = cc
# FC = ftn
CC = pgcc
# FC = pgf90
CC = icc
# FC = ifort -nofor-main
CC = gcc
# FC = gfortran -fcray-pointer
AR = ar ruc

########################################################################
# Compiling and linking options
########################################################################
#COPTS     = -g -pedantic -Wall 
#COPTS     = -g -cc=icc -config=icc -O3
#COPTS     = -g -config=icc -O3
#COPTS     = -g -O3

COPTS     = $(BOPT)
CINCLUDES = -I./include
CDEFS     = 
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)

# FOPTS     = $(BOPT)
# FINCLUDES = $(CINCLUDES)
# FDEFS     = $(CDEFS)
# FFLAGS    = $(FOPTS) $(FINCLUDES) $(FDEFS)

AMGLIB = lib/libAMG.a

LINKOPTS  = $(COPTS) #-Mnomain
LIBS      = -lm

CLFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++
FLFLAGS    = $(LINKOPTS) $(LIBS)


########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c

CSRCDIR = ./csrc
ELASDIR = ./elas
ELIPDIR = ./ellipticeqn

CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))
CSRC += $(foreach dir,$(ELASDIR),$(wildcard $(ELASDIR)/*.c))
CSRC += $(foreach dir,$(ELIPDIR),$(wildcard $(ELIPDIR)/*.c))

OBJSC := $(patsubst %.c,%.o,$(CSRC))

########################################################################
# List of all programs to be compiled
########################################################################

# Everything
#ALLPROG = $(TESTPROG) ($TESTMATPROG)
ALLPROG = $(AMGLIB) test testmat

# Test for solvers
TESTPROG = $(AMGLIB) main/test.o 

# Test for matrix properties
TESTMATPROG = $(AMGLIB) main/testmat.o


########################################################################
# Link
########################################################################

all: $(ALLPROG)

Default: test

$(AMGLIB): $(OBJSC)
	ranlib $(AMGLIB)

test: 
	$(CC) $(CFLAGS) -c main/test.c  -o main/test.o
	$(CC) main/test.o $(AMGLIB) $(FLFLAGS) -o test.ex

testmat: 
	$(CC) $(CFLAGS) -c main/testmat.c  -o main/testmat.o
	$(CC) main/testmat.o $(AMGLIB) $(FLFLAGS) -o testmat.ex


########################################################################
# Clean up
########################################################################

.PHONY : clean allclean

clean:
	rm -f main/*.o
	rm -f $(CSRCDIR)/*.o
	rm -f $(FSRCDIR)/*.o
	rm -f $(ELASDIR)/*.o
	rm -f $(ELIPDIR)/*.o
	rm -f lib/*.a

allclean:
	make clean
	rm -f *~
	rm -f main/*.o
	rm -f *.ex

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -c $<  -o $@
	$(AR) $(AMGLIB) $@
