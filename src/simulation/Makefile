# ============================================================
# Makefile — Ising model simulations
# ============================================================

SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --no-builtin-rules
.DEFAULT_GOAL := all

.PHONY: all clean \
        ising2d-cell ising2d-vertex \
        isingtr-cell isingtr-vertex

CC ?= gcc

CFLAGS = -std=c11 -O3 -Wall -pedantic -Iranlux-3.4

RANLUX = ranlux-3.4/ranlxs.o ranlux-3.4/ranlxd.o ranlux-3.4/ranlux_common.o
LIBS   = -lm -lstdc++

LATTICE_SIZES = 64 96 128 192 256 384 512 768 1024 1536 2048 3072

# ============================================================
# Main target
# ============================================================

all: $(RANLUX) \
     ising2d-cell \
     ising2d-vertex \
     isingtr-cell \
     isingtr-vertex
	@echo "Build completed."

# ============================================================
# Square lattice
# ============================================================

ising2d-cell: $(RANLUX)
	@echo "Building square-lattice spin-as-cell executables..."
	for L in $(LATTICE_SIZES); do \
	    $(CC) $(CFLAGS) \
	        -DDIM=2 \
	        -DL=$$L \
	        ising.c $(RANLUX) $(LIBS) \
	        -o ising2d-cell-$$L; \
	done

ising2d-vertex: $(RANLUX)
	@echo "Building square-lattice spin-as-vertex executables..."
	for L in $(LATTICE_SIZES); do \
	    $(CC) $(CFLAGS) \
	        -DDIM=2 \
	        -DSPIN_AS_VERTEX \
	        -DL=$$L \
	        ising.c $(RANLUX) $(LIBS) \
	        -o ising2d-vertex-$$L; \
	done

# ============================================================
# Triangular lattice
# ============================================================

isingtr-cell: $(RANLUX)
	@echo "Building triangular-lattice spin-as-cell executables..."
	for L in $(LATTICE_SIZES); do \
	    $(CC) $(CFLAGS) \
	        -DL=$$L \
	        ising.c $(RANLUX) $(LIBS) \
	        -o isingtr-cell-$$L; \
	done

isingtr-vertex: $(RANLUX)
	@echo "Building triangular-lattice spin-as-vertex executables..."
	for L in $(LATTICE_SIZES); do \
	    $(CC) $(CFLAGS) \
	        -DSPIN_AS_VERTEX \
	        -DL=$$L \
	        ising.c $(RANLUX) $(LIBS) \
	        -o isingtr-vertex-$$L; \
	done

# ============================================================
# Dependency compilation rules
# ============================================================

ranlux-3.4/%.o: ranlux-3.4/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

# ============================================================
# Clean
# ============================================================

clean:
	rm -f *.o
	rm -f ising2d-cell-*
	rm -f ising2d-vertex-*
	rm -f isingtr-cell-*
	rm -f isingtr-vertex-*
	rm -f ranlux-3.4/*.o