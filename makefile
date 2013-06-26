# ======== COMPILER ========

CC	= gcc
#OPT	= -Wall -O3 -DSO3_VERSION=\"0.1\" -DSO3_BUILD=\"`svnversion -n .`\"
OPT	= -Wall -g -DSO3_VERSION=\"0.1\" -DSO3_BUILD=\"`svnversion -n .`\"


# ======== LINKS ========

UNAME := $(shell uname)
PROGDIR = ..

ifeq ($(UNAME), Linux)
  MLAB		= /usr/local/MATLAB/R2011b
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib

  MEXEXT	= mexa64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif
ifeq ($(UNAME), Darwin)
  MLAB		= /Applications/MATLAB_R2011b.app
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib

  MEXEXT	= mexmaci64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif

SO3DIR   = $(PROGDIR)/so3
SO3LIB   = $(SO3DIR)/lib/c
SO3LIBNM = so3
SO3SRC   = $(SO3DIR)/src/c
SO3BIN   = $(SO3DIR)/bin/c
SO3OBJ   = $(SO3SRC)
SO3INC   = $(SO3DIR)/include/c
SO3DOC   = $(SO3DIR)/doc/c

SSHTDIR  = $(PROGDIR)/ssht
SSHTLIB  = $(SSHTDIR)/lib/c
SSHTLIBNM= ssht
SSHTBIN  = $(SSHTDIR)/bin/c
SSHTINC  = $(SSHTDIR)/include/c

ifeq ($(UNAME), Linux)
  FFTWDIR      = $(PROGDIR)/fftw-3.2.2_fPIC
endif
ifeq ($(UNAME), Darwin)
  FFTWDIR      = $(PROGDIR)/fftw
endif

FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3

SSHTSRCMAT	= $(SSHTDIR)/src/matlab
SSHTOBJMAT  	= $(SSHTSRCMAT)
SSHTOBJMEX  	= $(SSHTSRCMAT)

SO3SRCMAT	= $(SO3DIR)/src/matlab
SO3OBJMAT	= $(SO3SRCMAT)
SO3OBJMEX	= $(SO3SRCMAT)


# ======== SOURCE LOCATIONS ========

vpath %.c $(SO3SRC)
vpath %.h $(SO3SRC)
vpath %_mex.c $(SO3SRCMAT)


# ======== FFFLAGS ========

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC) -I$(SO3INC)
ifeq ($(UNAME), Linux)
  # Add -fPIC flag (required for mex build).
  # (Note that fftw must also be built with -fPIC.)
  FFLAGS += -fPIC
endif


# ======== LDFLAGS ========

LDFLAGS = -L$(SO3LIB) -l$(SO3LIBNM) -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWLIBNM) -lm

LDFLAGSMEX = -L$(SO3LIB) -l$(SO3LIBNM) -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWLIBNM)


# ======== OBJECT FILES TO MAKE ========

SO3OBJS = $(SO3OBJ)/so3_sampling.o    \
          $(SO3OBJ)/so3_core.o        \

SO3HEADERS = so3_types.h     \
             so3_error.h     \
	     so3_sampling.h  \
	     so3_core.h

SO3OBJSMAT = $(SO3OBJMAT)/so3_sampling_mex.o  \
             $(SO3OBJMAT)/so3_forward_mex.o   \
             $(SO3OBJMAT)/so3_inverse_mex.o

SO3OBJSMEX = $(SO3OBJMEX)/so3_sampling_mex.$(MEXEXT)  \
             $(SO3OBJMEX)/so3_forward_mex.$(MEXEXT)   \
             $(SO3OBJMEX)/so3_inverse_mex.$(MEXEXT)


# ======== MAKE RULES ========

$(SO3OBJ)/%.o: %.c $(SO3HEADERS)
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

.PHONY: default
default: lib test about

.PHONY: test
test: $(SO3BIN)/so3_test about
$(SO3BIN)/so3_test: $(SO3OBJ)/so3_test.o $(SO3LIB)/lib$(SO3LIBNM).a
	$(CC) $(OPT) $< -o $(SO3BIN)/so3_test $(LDFLAGS) 

.PHONY: about
about: $(SO3BIN)/so3_about
$(SO3BIN)/so3_about: $(SO3OBJ)/so3_about.o 
	$(CC) $(OPT) $< -o $(SO3BIN)/so3_about

.PHONY: runtest
runtest: test
	$(SO3BIN)/so3_test 64 0

.PHONY: all
all: lib test about matlab


# Library

.PHONY: lib
lib: $(SO3LIB)/lib$(SO3LIBNM).a
$(SO3LIB)/lib$(SO3LIBNM).a: $(SO3OBJS)
	ar -r $(SO3LIB)/lib$(SO3LIBNM).a $(SO3OBJS)


# Matlab

$(SO3OBJMAT)/%_mex.o: %_mex.c $(SO3LIB)/lib$(SO3LIBNM).a
	echo "hmmmm"
	#$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC} 

$(SO3OBJMEX)/%_mex.$(MEXEXT): $(SO3OBJMAT)/%_mex.o $(SO3LIB)/lib$(SO3LIBNM).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

.PHONY: matlab
matlab: $(SO3OBJSMEX)


# Documentation 

.PHONY: doc
doc:
	doxygen $(SO3SRC)/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -f $(SO3DOC)/html/*


# Cleaning up

.PHONY: clean
clean:	tidy
	rm -f $(SO3OBJ)/*.o
	rm -f $(SO3LIB)/lib$(SO3LIBNM).a
	rm -f $(SO3BIN)/so3_test
	rm -f $(SO3BIN)/so3_about
	rm -f $(SO3OBJMAT)/*.o
	rm -f $(SO3OBJMEX)/*.$(MEXEXT)

.PHONY: tidy
tidy:
	rm -f *~ 

.PHONY: cleanall
cleanall: clean cleandoc
