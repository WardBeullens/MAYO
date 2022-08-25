CC=gcc
CFLAGS=-I ../XKCP/bin/AVX2/
LFLAGS=-L ../XKCP/bin/AVX2/ -lXKCP -lcrypto

IMPLEMENTATION_SOURCE = MAYO.c arithmetic.c test.c  
IMPLEMENTATION_HEADERS= MAYO.h arithmetic.h MAYO_params.h keccaklib  

test: $(IMPLEMENTATION_SOURCE) $(IMPLEMENTATION_HEADERS)
	gcc -o test $(IMPLEMENTATION_SOURCE) $(CFLAGS) $(LFLAGS) -std=c11 -O3 -g -march=native -funroll-loops -Wall

keccaklib: 
	(cd ../XKCP; make AVX2/libXKCP.a)

.PHONY: clean
clean:
	rm -f test >/dev/null
