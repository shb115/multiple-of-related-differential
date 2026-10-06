CC     := gcc
CFLAGS := -O3 -Wall

ifeq ($(OS),Windows_NT)
$(error Run this artifact from WSL/Linux. Native Windows builds create .exe artifacts and are intentionally unsupported.)
endif

all: count_quartets_nz2 count_quartets_nz1 count_quartets_nz0

count_quartets_nz2: count_quartets_nz2.c small_aes.c small_aes.h ../common/rmd_random.h
	$(CC) $(CFLAGS) count_quartets_nz2.c small_aes.c -o $@

count_quartets_nz1: count_quartets_nz1.c small_aes.c small_aes.h ../common/rmd_random.h
	$(CC) $(CFLAGS) count_quartets_nz1.c small_aes.c -o $@

count_quartets_nz0: count_quartets_nz0.c small_aes.c small_aes.h ../common/rmd_random.h
	$(CC) $(CFLAGS) count_quartets_nz0.c small_aes.c -o $@

clean:
	rm -f count_quartets_nz2 count_quartets_nz1 count_quartets_nz0

.PHONY: all clean
