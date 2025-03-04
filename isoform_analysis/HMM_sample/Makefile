CC = gcc
CFLAGS = -Wall -Wextra -g

OBJS = main.o HMM_algo.o

all: conventional_HMM

conventional_HMM: $(OBJS)
	$(CC) $(CFLAGS) -o conventional_HMM $(OBJS) -lm
	@chmod +x conventional_HMM

main.o: main.c HMM_algo.c HMM_model.h
	@echo "Compling main"
	$(CC) $(CFLAGS) -c main.c

HMM_sample.o: HMM_algo.c
	@echo "Compling HMM algorithms"
	$(CC) $(CFLAGS) -c HMM_algo.c

clean: 
	@echo "Removing everything but the source file"
	rm -f conventional_HMM *.o

.PHONY: all clean