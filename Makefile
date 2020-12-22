CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDFLAGS = -lpthread -lm
ALL = EmerRate Classifier

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

EmerRate: EmerRate.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o $@ EmerRate.c libfastk.c $(LDFLAGS)

Classifier: Classifier.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ Classifier.c libfastk.c DB.c QV.c $(LDFLAGS)

clean:
	rm -f $(ALL)
