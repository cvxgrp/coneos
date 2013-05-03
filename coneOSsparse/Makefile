# MAKEFILE for coneOS
include coneOS.mk

OBJECTS = coneOS.o linAlg.o util.o cones.o cs.o
AMD_SOURCE = $(wildcard direct/amd_*.c)
DIRECT_OBJECTS = direct/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = bin/demo_direct bin/demo_indirect

.PHONY: default 
default: lib/libconeOSdir.a lib/libconeOSindir.a bin/demo_direct bin/demo_indirect

coneOS.o 	: coneOS.h
linAlg.o 	: linAlg.h
util.o		: util.h
cones.o		: cones.c
cs.o		: cs.h

direct/private.o		: direct/private.h
direct/ldl.o			: direct/ldl.h
direct/amd_1.o			: direct/amd_internal.h direct/amd.h
direct/amd_2.o			: direct/amd_internal.h direct/amd.h
direct/amd_aat.o		: direct/amd_internal.h direct/amd.h
direct/amd_control.o	: direct/amd_internal.h direct/amd.h
direct/amd_defaults.o 	: direct/amd_internal.h direct/amd.h
direct/amd_dump.o		: direct/amd_internal.h direct/amd.h
direct/amd_global.o		: direct/amd_internal.h direct/amd.h
direct/amd_info.o		: direct/amd_internal.h direct/amd.h
direct/amd_order.o		: direct/amd_internal.h direct/amd.h
direct/amd_post_tree.o	: direct/amd_internal.h direct/amd.h
direct/amd_postorder.o	: direct/amd_internal.h direct/amd.h
direct/amd_preprocess.o	: direct/amd_internal.h direct/amd.h
direct/amd_valid.o		: direct/amd_internal.h direct/amd.h

indirect/private.o	: indirect/private.c
	$(CC) $(CFLAGS) $(OPENFLAG) -c -o $@ $^ 

lib/libconeOSdir.a: $(OBJECTS) direct/private.o  $(DIRECT_OBJECTS)
	mkdir -p lib
	$(ARCHIVE) lib/libconeOSdir.a $^
	- $(RANLIB) lib/libconeOSdir.a

lib/libconeOSindir.a: $(OBJECTS) indirect/private.o
	mkdir -p lib
	$(ARCHIVE) lib/libconeOSindir.a $^
	- $(RANLIB) lib/libconeOSindir.a

bin/demo_direct: run_coneOS.c lib/libconeOSdir.a
	mkdir -p bin
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/data_coneOS\"" -o $@ $^ $(LDFLAGS) 

bin/demo_indirect: run_coneOS.c lib/libconeOSindir.a
	mkdir -p bin
	$(CC) $(CFLAGS) $(OPENFLAG) -DDEMO_PATH="\"$(CURDIR)/data_coneOS\"" -o $@ $^ $(LDFLAGS) 

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS) $(DIRECT_OBJECTS) direct/private.o indirect/private.o 
	@rm -rf bin/*.dSYM
	@rm -rf matlab/*.mex*
purge: clean 
	@rm -rf bin lib
