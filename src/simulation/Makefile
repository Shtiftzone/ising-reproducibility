# ============================================================
# Makefile — Ising model simulations
# ============================================================

SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --no-builtin-rules
.DEFAULT_GOAL := all

.PHONY: all clean ising2d isingtr

CC       ?= gcc

CFLAGS   = -std=c11 -O3 -Wall -pedantic -Iranlux-3.4

RANLUX = ranlux-3.4/ranlxs.o ranlux-3.4/ranlxd.o ranlux-3.4/ranlux_common.o
LIBS   = -lm -lstdc++

LATTICE_SIZES_2D = 64 96 128 192 256 384 512 768 1024 1536 2048 3072

# ============================================================
# Main targets
# ============================================================

all: $(RANLUX) ising2d isingtr
	@echo "Build completed."

ising2d: $(RANLUX)
	@echo "Building 2D Ising executables..."
	for L in $(LATTICE_SIZES_2D); do \
	    $(CC) $(CFLAGS) -DDIM=2 -DL=$$L ising.c $(RANLUX) $(LIBS) -o ising2d-$$L; \
	done

isingtr: $(RANLUX)
	@echo "Building triangular-lattice Ising executables..."
	for L in $(LATTICE_SIZES_2D); do \
	    $(CC) $(CFLAGS) -DL=$$L ising.c $(RANLUX) $(LIBS) -o isingtr-$$L; \
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
	rm -f ising2d-* ising2d4-*
	rm -f isingtr-* isingtr6-*
	rm -f ising3d-* ising3d6-*
	rm -f ising4d-* ising4d8-*
	rm -f ranlux-3.4/*.o