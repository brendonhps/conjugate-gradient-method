  PROG   = cgSolver
  MODULOS   = utils\
  			  gradienteconj\

           CC = gcc -std=c11 -O3 -g -mavx -march=native -ffast-math

         OBJS = $(addsuffix .o,$(MODULOS))

       LIKWID = /usr/local#/home/soft/likwid
       LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID)/include
       LIKWID_LIBS = -L$(LIKWID)/lib -llikwid

       CFLAGS = $(LIKWID_FLAGS)
       LFLAGS = $(LIKWID_LIBS) -lm

.PHONY: all clean limpa purge faxina distclean

%.o: %.c %.h
	$(CC) $(LFLAGS) $(CFLAGS) -c $<

all: $(PROG)

$(PROG): $(PROG).o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	@echo "Limpando ...."
	@rm -f *~ *.bak *.tmp

purge distclean:   clean
	@echo "Faxina ...."
	@rm -f  $(PROG) *.o core a.out
	@rm -f *.png marker.out *.log
