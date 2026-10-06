CC     := gcc
CFLAGS := -O3 -maes -mpclmul -fopenmp -Wall

ifeq ($(OS),Windows_NT)
$(error Run this artifact from WSL/Linux. Native Windows builds create .exe artifacts and are intentionally unsupported.)
endif

all: count_quartets_nz2 count_quartets_nz1 count_quartets_nz0

count_quartets_nz2: count_quartets_nz2.c ../common/rmd_aesni.h ../common/rmd_random.h
	$(CC) $(CFLAGS) count_quartets_nz2.c -o $@ -lm

count_quartets_nz1: count_quartets_nz1.c ../common/rmd_aesni.h ../common/rmd_random.h
	$(CC) $(CFLAGS) count_quartets_nz1.c -o $@ -lm

count_quartets_nz0: count_quartets_nz0.c ../common/rmd_aesni.h ../common/rmd_random.h
	$(CC) $(CFLAGS) count_quartets_nz0.c -o $@ -lm

clean:
	rm -f count_quartets_nz2 count_quartets_nz1 count_quartets_nz0

.PHONY: all clean
