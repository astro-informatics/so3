# ======== COMPILER ========


CC	= gcc

#OPT	= -Wall -O3 -fopenmp -DSO3_VERSION=\"0.1\" -DSO3_BUILD=\"`git rev-parse HEAD`\"
OPT	= -Wall -g -fopenmp -DSO3_VERSION=\"1.3.0\" -DSO3_BUILD=\"`git rev-parse HEAD`\"


# ======== LINKS ========

UNAME := $(shell uname)
PROGDIR = ..

MLAB		= ${MATLAB}
MLABINC	= ${MLAB}/extern/include
MLABLIB	= ${MLAB}/extern/lib
# -------------------- 
ifeq ($(UNAME), Linux)
  MEXEXT	= mexa64
endif
ifeq ($(UNAME), Darwin)
  MEXEXT	= mexmaci64
endif
# -------------------- 
MEX 		= ${MLAB}/bin/mex
MEXFLAGS	= -cxx

SO3DIR   = $(PROGDIR)/so3
SO3LIB   = $(SO3DIR)/lib/c
SO3LIBNM = so3
SO3SRC   = $(SO3DIR)/src/c
SO3BIN   = $(SO3DIR)/bin/c
SO3OBJ   = $(SO3SRC)
SO3INC   = $(SO3DIR)/include/c
SO3DOC   = $(SO3DIR)/docs/c

SSHTDIR  = $(PROGDIR)/ssht
SSHTLIB  = $(SSHTDIR)/lib/c
SSHTLIBNM= ssht
SSHTBIN  = $(SSHTDIR)/bin/c
SSHTINC  = $(SSHTDIR)/include/c

FFTWDIR      = $(FFTW)
FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3
#FFTWOMPLIBNM = fftw3_threads

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


#LDFLAGS = -L$(SO3LIB) -l$(SO3LIBNM) -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWOMPLIBNM) -l$(FFTWLIBNM) -lm

#LDFLAGSMEX = -L$(SO3LIB) -l$(SO3LIBNM) -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWOMPLIBNM) -l$(FFTWLIBNM)


# ======== OBJECT FILES TO MAKE ========

SO3OBJS = $(SO3OBJ)/so3_sampling.o    \
          $(SO3OBJ)/so3_core.o        \
          $(SO3OBJ)/so3_adjoint.o

SO3HEADERS = so3_types.h     \
             so3_error.h     \
             so3_sampling.h  \
             so3_core.h 	 \
             so3_adjoint.h

SO3OBJSMAT = $(SO3OBJMAT)/so3_sampling_mex.o \
             $(SO3OBJMAT)/so3_elmn2ind_mex.o \
             $(SO3OBJMAT)/so3_ind2elmn_mex.o \
             $(SO3OBJMAT)/so3_forward_mex.o  \
             $(SO3OBJMAT)/so3_inverse_mex.o  \
             $(SO3OBJMAT)/so3_forward_direct_mex.o  \
             $(SO3OBJMAT)/so3_inverse_direct_mex.o  \
             $(SO3OBJMAT)/so3_forward_adjoint_direct_mex.o  \
             $(SO3OBJMAT)/so3_inverse_adjoint_direct_mex.o


SO3OBJSMEX = $(SO3OBJMEX)/so3_sampling_mex.$(MEXEXT) \
             $(SO3OBJMEX)/so3_elmn2ind_mex.$(MEXEXT) \
             $(SO3OBJMEX)/so3_ind2elmn_mex.$(MEXEXT) \
             $(SO3OBJMEX)/so3_forward_mex.$(MEXEXT)  \
             $(SO3OBJMEX)/so3_inverse_mex.$(MEXEXT)  \
             $(SO3OBJMEX)/so3_forward_direct_mex.$(MEXEXT)  \
             $(SO3OBJMEX)/so3_inverse_direct_mex.$(MEXEXT)  \
             $(SO3OBJMEX)/so3_forward_adjoint_direct_mex.$(MEXEXT)  \
             $(SO3OBJMEX)/so3_inverse_adjoint_direct_mex.$(MEXEXT)



# ======== MAKE RULES ========

$(SO3OBJ)/%.o: %.c $(SO3HEADERS)
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

.PHONY: default
default: lib unittest test about

.PHONY: unittest
unittest: $(SO3BIN)/unittest/so3_unittest
$(SO3BIN)/unittest/so3_unittest: $(SO3OBJ)/unittest/so3_unittest.o $(SO3LIB)/lib$(SO3LIBNM).a
	$(CC) $(OPT) $< -o $(SO3BIN)/unittest/so3_unittest $(LDFLAGS)

.PHONY: rununittest
rununittest: unittest
	$(SO3BIN)/unittest/so3_unittest

.PHONY: test
test: $(SO3BIN)/so3_test about
$(SO3BIN)/so3_test: $(SO3OBJ)/so3_test.o $(SO3LIB)/lib$(SO3LIBNM).a
	$(CC) $(OPT) $< -o $(SO3BIN)/so3_test $(LDFLAGS)

.PHONY: test_csv
test_csv: $(SO3BIN)/so3_test_csv about
$(SO3BIN)/so3_test_csv: $(SO3OBJ)/so3_test_csv.o $(SO3LIB)/lib$(SO3LIBNM).a
	$(CC) $(OPT) $< -o $(SO3BIN)/so3_test_csv $(LDFLAGS)

.PHONY: about
about: $(SO3BIN)/so3_about
$(SO3BIN)/so3_about: $(SO3OBJ)/so3_about.o
	$(CC) $(OPT) $< -o $(SO3BIN)/so3_about

.PHONY: runtest
runtest: test
	$(SO3BIN)/so3_test

.PHONY: all
all: lib unittest test test_csv about matlab


# Library

.PHONY: lib
lib: $(SO3LIB)/lib$(SO3LIBNM).a
$(SO3LIB)/lib$(SO3LIBNM).a: $(SO3OBJS)
	ar -r $(SO3LIB)/lib$(SO3LIBNM).a $(SO3OBJS)


# Matlab

$(SO3OBJMAT)/%_mex.o: %_mex.c $(SO3LIB)/lib$(SO3LIBNM).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC}

$(SO3OBJMEX)/%_mex.$(MEXEXT): $(SO3OBJMAT)/%_mex.o $(SO3LIB)/lib$(SO3LIBNM).a
	$(MEX) $< -output $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

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
	rm -f $(SO3OBJ)/unittest/*.o
	rm -f $(SO3LIB)/lib$(SO3LIBNM).a
	rm -f $(SO3BIN)/so3_test
	rm -f $(SO3BIN)/so3_about
	rm -f $(SO3BIN)/unittest/so3_unittest
	rm -f $(SO3OBJMAT)/*.o
	rm -f $(SO3OBJMEX)/*.$(MEXEXT)

.PHONY: tidy
tidy:
	rm -f *~

.PHONY: cleanall
cleanall: clean cleandoc

