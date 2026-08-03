# ============================================================
#  Makefile — Ising models with hybrid Metropolis and Wolff
# ============================================================

SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --no-builtin-rules
.DEFAULT_GOAL := all
.PHONY: all clean

CC       ?= gcc
CXX      ?= g++

CFLAGS    = -std=c11 -O3 -Wall -pedantic -Iranlux-3.4
CXXFLAGS  = -std=c++14 -O3 -Wall -pedantic -Igudhi/include

# Boost z modu³u
ifdef BOOST_ROOT
CXXFLAGS += -I$(BOOST_ROOT)/include
endif

RANLUX    = ranlux-3.4/ranlxs.o ranlux-3.4/ranlxd.o ranlux-3.4/ranlux_common.o
LIBS      = -lm -lstdc++ -L/mnt/d/PhD/IsingModels/gudhi/lib

# ============================================================
# G³ówne cele
# ============================================================
all: $(RANLUX) betti3d.o
	@echo "== Budujê 2D =="
	set -x; for L in 64 96 128 192 256 384 512 768 1024 1536 2048 3072; do \
	    $(CC) $(CFLAGS) -DDIM=2 -DL=$$L ising.c $(RANLUX) $(LIBS) -o ising2d-$$L; \
	    $(CC) $(CFLAGS) -DDIM=2 -DCELL -DL=$$L ising.c $(RANLUX) $(LIBS) -o ising2d4-$$L; \
	    $(CC) $(CFLAGS) -DL=$$L ising.c $(RANLUX) $(LIBS) -o isingtr-$$L; \
	    $(CC) $(CFLAGS) -DCELL -DL=$$L ising.c $(RANLUX) $(LIBS) -o isingtr6-$$L; \
	done

	@echo "== Budujê 3D =="
	set -x; for L in 16 24 32 48 64 96 128 192 256 384; do \
	    $(CC) $(CFLAGS) -DDIM=3 -DL=$$L ising.c $(RANLUX) betti3d.o $(LIBS) -o ising3d-$$L; \
	    $(CC) $(CFLAGS) -DDIM=3 -DCELL -DL=$$L ising.c $(RANLUX) $(LIBS) -o ising3d6-$$L; \
	done

	@echo "== Budujê 4D =="
	set -x; for L in 8 12 16 24 32 48 60 64 72 78 84 96 108 128; do \
	    $(CC) $(CFLAGS) -DDIM=4 -DL=$$L ising.c $(RANLUX) $(LIBS) -o ising4d-$$L; \
	    $(CC) $(CFLAGS) -DDIM=4 -DCELL -DL=$$L ising.c $(RANLUX) $(LIBS) -o ising4d8-$$L; \
	done

# ============================================================
# Regu³y kompilacji zale¿noœci
# ============================================================
ranlux-3.4/%.o: ranlux-3.4/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

betti3d.o: betti3d.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# ============================================================
# Clean
# ============================================================
clean:
	rm -f *.o ising*-*
	rm -f ranlux-3.4/*.o
