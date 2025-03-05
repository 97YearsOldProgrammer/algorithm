CC = gcc
CFLAGS = -Wall -Wextra -g

OBJS = main.o donor_acceptor_parser.o model.o

all: EDHMM

EDHMM: $(OBJS)
	$(CC) $(CFLAGS) -o EDHMM $(OBJS) -lm
	@chmod +x EDHMM

main.o: main.c donor_acceptor_parser.c model.c model.h
	@echo "Compling main"
	$(CC) $(CFLAGS) -o main.c

donor_acceptor_parser.o: donor_acceptor_parser.c model.h
	@echo "Compling emission probability"
	$(CC) $(CFLAGS) -c donor_acceptor_parser.c

clean: 
	@echo "Removing everything but the source file"
	rm -f EDHMM *.o

.PHONY: all clean