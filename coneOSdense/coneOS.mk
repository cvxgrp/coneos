UNAME = $(shell uname -s)
CC = gcc
CFLAGS = -g -Wall -pedantic -O3 -I.
LDFLAGS = -lm 
OPENFLAG = #-fopenmp

# LINK YOUR INSTALLED BLAS LIBRARY BELOW:
CFLAGS += -I/opt/local/include -I/usr/local/include 
LDFLAGS += -L/opt/local/lib -L/usr/local/lib

# standard blas libs
#LDFLAGS += -llapack -lcblas -lgfortran

# openblas
LDFLAGS += -lopenblas

# ATLAS, multi-thread (if atlas compiled with multi-thread support):
#LDFLAGS += -L/usr/local/atlas/lib -lpthread -lptcblas -lptf77blas -llapack -latlas

# ATLAS, single thread:
#LDFLAGS += -L/usr/local/atlas/lib -lcblas -lf77blas -llapack -latlas

CFLAGS += -DLAPACK_LIB_FOUND

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib
