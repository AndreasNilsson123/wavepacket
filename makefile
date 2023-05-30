#
# Sample makefile for the compilation of program "wavepacket"
# on a parallel computer
#
# C.M. Dion, A. Hashemloo, and G. Rahali
# April 2013
#

CC = gcc

CFLAGS = -O3 -std=c99 

LDLIBS = -lm -lgsl -lgslcblas

sim_bonus:  sim.o wavepacket_bonus.o 
	$(CC) $(CFLAGS) -o sim_bonus $^ $(LDLIBS)

.PHONY: clean
clean:	
	rm -f *.o sim_bonus