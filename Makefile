CC = gcc
CFLAGS = -g -O0 -Wall -Wextra -Wno-unused-function
LDFLAGS = -lm -lz -lpthread

ALL = ClassPro ClassGS prof2class

all: $(ALL)

ClassPro: ClassPro.c ClassPro.h const.c io.c emodel.c hist.c context.c wall.c class.c bessel.c bessel.h prob.c prob.h kseq.h libfastk.c libfastk.h gene_core.c gene_core.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ ClassPro.c libfastk.c DB.c QV.c $(LDFLAGS)

ClassGS: ClassGS.c libfastk.c libfastk.h gene_core.c gene_core.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ ClassGS.c libfastk.c DB.c QV.c $(LDFLAGS)

prof2class: prof2class.c libfastk.c libfastk.h gene_core.c gene_core.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ prof2class.c libfastk.c DB.c QV.c $(LDFLAGS)

clean:
	rm -f $(ALL)
