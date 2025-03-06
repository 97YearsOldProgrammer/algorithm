CC = gcc
CFLAGS = -Wall -Wextra -g

OBJS = main.o parser.o model.o

all: EDHMM

EDHMM: $(OBJS)
	$(CC) $(CFLAGS) -o EDHMM $(OBJS) -lm
	@chmod +x EDHMM

main.o: main.c parser.c model.c model.h
	@echo "Compling main"
	$(CC) $(CFLAGS) -o main.c

parser.o: parser.c model.h
	@echo "Compling function and data structure for computation"
	$(CC) $(CFLAGS) -o parser.c

model.o: model.c model.h
	@echo "Compling model"
	$(CC) $(CFLAGS) -o model.c

clean: 
	@echo "Removing everything but the source file"
	rm -f EDHMM *.o

.PHONY: all clean