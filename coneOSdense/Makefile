# MAKEFILE for coneOS
# change blas library in coneOS.mk
include coneOS.mk

OBJECTS = coneOS.o util.o cones.o linAlg.o
TARGETS = bin/demo_direct bin/demo_indirect

.PHONY: default 
default: lib/libconeOSdensedir.a bin/demo_direct lib/libconeOSdenseindir.a bin/demo_indirect

coneOS   	: coneOS.c coneOS.h normalize.h
util		: util.c util.h
cones		: cones.c cones.h
linAlg      : linAlg.c linAlg.h

direct/private.o	: direct/private.c
indirect/private.o	: indirect/private.c 	

lib/libconeOSdensedir.a: $(OBJECTS) direct/private.o
	mkdir -p lib
	$(ARCHIVE) lib/libconeOSdensedir.a $^
	- $(RANLIB) lib/libconeOSdensedir.a

lib/libconeOSdenseindir.a: $(OBJECTS) indirect/private.o
	mkdir -p lib
	$(ARCHIVE) lib/libconeOSdenseindir.a $^
	- $(RANLIB) lib/libconeOSdenseindir.a

bin/demo_direct: run_coneOS.c lib/libconeOSdensedir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_dense\"" -o $@ $^ $(LDFLAGS)

bin/demo_indirect: run_coneOS.c lib/libconeOSdenseindir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_dense\"" -o $@ $^ $(LDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) direct/private.o indirect/private.o
	@rm -rf bin/*.dSYM
	@rm -rf matlab/*.mex*
purge: clean 
	@rm -rf bin lib
