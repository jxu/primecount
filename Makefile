CC = gcc
CCFLAGS = -Wall -Wextra -lm -fopenmp

.PHONY: primecount tests
# Target-specific variable values
release: CCFLAGS += -O3
release: primecount

debug: CCFLAGS += -Og -g
debug: tests

primecount: main.c
	$(CC) $(CCFLAGS) -o $@ $<

tests: tests.c
	$(CC) $(CCFLAGS) -o $@ $<
