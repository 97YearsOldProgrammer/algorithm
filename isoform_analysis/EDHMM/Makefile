CC = gcc
CFLAGS = -Wall -Wextra -g

OBJS = main.o parser.o model.o

all: EDHMM

EDHMM: $(OBJS)
	$(CC) $(CFLAGS) -o EDHMM $(OBJS) -lm
	@chmod +x EDHMM

main.o: main.c parser.c model.c model.h
	@echo "Compling main"
	$(CC) $(CFLAGS) -c main.c

parser.o: parser.c model.h
	@echo "Compling function and data structure for computation"
	$(CC) $(CFLAGS) -c parser.c

model.o: model.c model.h
	@echo "Compling model"
	$(CC) $(CFLAGS) -c model.c

clean: 
	@echo "Removing everything but the source file"
	rm -f EDHMM *.o

test: EDHMM
	./EDHMM seq
	
.PHONY: all clean