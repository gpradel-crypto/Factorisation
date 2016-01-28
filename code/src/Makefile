#Variables
EXE = factorisation
SRC = factorisation.c
OBJS = $(SRC:.c=.o)

# Compilation flags
CFLAGS = -Wall -Wextra -std=c99 -g -fprofile-arcs -ftest-coverage
LDFLAGS = -L/usr/local/lib/ -lgmp -lm

# Special rules and targets
.PHONY : clean help all superclean

# Rules and targets
all : $(EXE)

$(EXE) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJS) : $(SRC) factorisation.h
	$(CC) $(CFLAGS) -c $(SRC)

help :
	@echo "Available targets :"
	@echo "all : compile factorisation"
	@echo "help : display this help"
	@echo "clean : delete all temporary files"
	@echo "softclean : delete all temporary files except the executable"

softclean :
	rm -rf $(OBJS)
	rm -rf *~
	rm -rf \#*\#
	rm -rf *.h.gch
	rm -rf *.gcno
	rm -rf *.gcda

clean : softclean
	rm -rf $(EXE)		
