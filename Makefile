CC = gcc
CFLAGS = -g -O0 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing -Wno-unused-function
LDFLAGS = -lpthread -lm -lz
ALL = ClassPro ClassGS prof2class

all: $(ALL)

libfastk.c: gene_core.c
libfastk.h: gene_core.h

ClassPro: ClassPro.c ClassPro.h class.c io.c const.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h bessel.c bessel.h kseq.h
	$(CC) $(CFLAGS) -o $@ ClassPro.c class.c io.c const.c libfastk.c DB.c QV.c bessel.c $(LDFLAGS)

ClassGS: ClassGS.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ ClassGS.c libfastk.c DB.c QV.c $(LDFLAGS)

prof2class: prof2class.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ prof2class.c libfastk.c DB.c QV.c $(LDFLAGS)

clean:
	rm -f $(ALL)
