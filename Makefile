CC = gcc
CCFLAGS = -Wall -Wextra
LDFLAGS = -lm

# Target-specific variable values
# valgrind reports openmp as still reachable leak
release: CCFLAGS += -O3 -fopenmp
release: primecount

debug: CCFLAGS += -Og -g -fsanitize=undefined
debug: primecount

primecount: primecount.c
	$(CC) $(CCFLAGS) -o $@ $< $(LDFLAGS)

