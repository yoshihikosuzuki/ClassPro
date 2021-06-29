CC = gcc
CFLAGS = -g -O0 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDFLAGS = -lpthread -lm
ALL = ClassPro

all: $(ALL)

libfastk.c: gene_core.c
libfastk.h: gene_core.h

ClassPro: ClassPro.c ClassPro.h class.c io.c const.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h bessel.c bessel.h
	$(CC) $(CFLAGS) -o $@ ClassPro.c class.c io.c const.c libfastk.c DB.c QV.c bessel.c $(LDFLAGS)

clean:
	rm -f $(ALL)
