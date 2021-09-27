CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-function
LDFLAGS = -lm -lz -lpthread

ALL = ClassPro ClassGS prof2class class2acc class2cns
GENE_LIBS = libfastk.c libfastk.h gene_core.c gene_core.h DB.c DB.h QV.c QV.h
GENE_LIBS_C = libfastk.c DB.c QV.c

all: $(ALL)

ClassPro: ClassPro.c ClassPro.h const.c io.c emodel.c hist.c context.c wall.c class_rel.c class_unrel.c bessel.c bessel.h prob.h util.h kseq.h $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ ClassPro.c $(GENE_LIBS_C) $(LDFLAGS)

ClassGS: ClassGS.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ ClassGS.c $(GENE_LIBS_C) $(LDFLAGS)

prof2class: prof2class.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ prof2class.c $(GENE_LIBS_C) $(LDFLAGS)

class2acc: class2acc.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ class2acc.c $(GENE_LIBS_C) $(LDFLAGS)

class2cns: class2cns.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ class2cns.c $(GENE_LIBS_C) $(LDFLAGS)

clean:
	rm -f $(ALL)
