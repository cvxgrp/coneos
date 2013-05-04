UNAME = $(shell uname -s)
CC = gcc
CFLAGS = -g -Wall -pedantic -O3 -I.
LDFLAGS = -lm 
OPENFLAG = -fopenmp

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib
