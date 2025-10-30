CC = gcc
CCFLAGS = -Wall -Wextra -fopenmp
LDFLAGS = -lm

# Target-specific variable values
release: CCFLAGS += -O3
release: primecount

debug: CCFLAGS += -Og -g
debug: primecount

primecount: primecount.c
	$(CC) $(CCFLAGS) -o $@ $< $(LDFLAGS)

