CC = gcc
CFLAGS = -g -O0 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDFLAGS = -lpthread -lm
ALL = Classifier

all: $(ALL)

libfastk.c: gene_core.c
libfastk.h: gene_core.h

Classifier: Classifier.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h bessel.c bessel.h
	$(CC) $(CFLAGS) -o $@ Classifier.c libfastk.c DB.c QV.c bessel.c $(LDFLAGS)

clean:
	rm -f $(ALL)
