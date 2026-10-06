# Makefile for the 5-round AES distinguishers (Sections 4 and 5 of the paper).
#
#   cp_nz01 : chosen-plaintext (CP) distinguisher, pairings n_z=0 and n_z=1   (Section 4.2)
#   cp_nz2  : chosen-plaintext (CP) distinguisher, pairing  n_z=2             (Section 4.2)
#   acp_all : adaptively-chosen-plaintext (ACP) distinguisher, all pairings   (Section 5.2)
#   cp_nz012: CP, pairings n_z=1 on {0,10} (Yan et al.), n_z=1 on {5,15} and n_z=0 on one
#             structure, with their union (Section 4.3); cp_nz012f is a faster AES-NI build
#             with the same T lines and footer (only the program name in the header line
#             differs)
#
# Build:   make
#
# Usage examples (5-round AES, w = 8):
#   ./cp_nz01 11664 108 108 5 300     # n_z=0 and n_z=1 on one structure ~2^27
#                                     #   args: s0 s1 s2 rounds [trials] [seed]
#                                     #   structure = {0,10}=s0 x {5}=s1 x {15}=s2
#   ./cp_nz2  11664 11664 5 300       # n_z=2 on one structure ~2^27
#                                     #   args: s0 s1 rounds [trials] [seed]
#                                     #   structure = {0,10}=s0 x {5,15}=s1
#   ./acp_all 1 23.32 5 300 2024      # ACP, pairing nz, data complexity 2^Dexp
#                                     #   args: nz Dexp rounds [trials] [seed]
#   ./cp_nz012f 9153 96 96 5 30 42 0  # Section 4.3, three pairings on one structure
#                                     #   args: s0 s1 s2 rounds trials seed [t_start]
#                                     #   (cp_nz012: same T lines and footer, adds [verify] after t_start)
#
# Set rounds=10 for the random-permutation baseline.

CC      = gcc
CFLAGS  = -O3 -march=native
LDLIBS  = -lm

ifeq ($(OS),Windows_NT)
$(error Run this artifact from WSL/Linux. Native Windows builds create .exe artifacts and are intentionally unsupported.)
endif

all: cp_nz01 cp_nz2 acp_all cp_nz012 cp_nz012f

cp_nz01: cp_nz01.c aes.c aes.h
	$(CC) $(CFLAGS) cp_nz01.c aes.c -o $@ $(LDLIBS)

cp_nz2: cp_nz2.c aes.c aes.h
	$(CC) $(CFLAGS) cp_nz2.c aes.c -o $@ $(LDLIBS)

acp_all: acp_all.c aes.c aes.h
	$(CC) $(CFLAGS) acp_all.c aes.c -o $@ $(LDLIBS)

cp_nz012: cp_nz012.c aes.c aes.h
	$(CC) $(CFLAGS) cp_nz012.c aes.c -o $@ $(LDLIBS)

cp_nz012f: cp_nz012f.c aes.c aes.h ../common/rmd_aesni.h
	$(CC) $(CFLAGS) cp_nz012f.c aes.c -o $@ $(LDLIBS)

clean:
	rm -f cp_nz01 cp_nz2 acp_all cp_nz012 cp_nz012f cp_nz012_m8 cp_nz012_m16

.PHONY: all clean
