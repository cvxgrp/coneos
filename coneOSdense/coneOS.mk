UNAME = $(shell uname -s)
CC = gcc
CFLAGS = -g -Wall -pedantic -O3 -I. -Iexternal/ 
LDFLAGS = -lm 
OPENFLAG = #-fopenmp

# LINK YOUR INSTALLED BLAS LIBRARY BELOW:
LAFLAGS = -I/opt/local/include -I/usr/local/include
LALIBS = -L/opt/local/lib -L/usr/local/lib 

# standard blas libs
#LALIBS += -llapack -lcblas -lgfortran

# openblas
LALIBS += -lopenblas

# ATLAS, multi-thread (if atlas compiled with multi-thread support):
#LALIBS += -L/usr/local/atlas/lib -lpthread -lptcblas -lptf77blas -llapack -latlas

# ATLAS, single thread:
#LALIBS += -L/usr/local/atlas/lib -lcblas -lf77blas -llapack -latlas

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib
