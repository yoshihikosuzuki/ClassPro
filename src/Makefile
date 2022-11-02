CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wno-unused-function
LIBS = -Igsl/include -Lgsl/lib -lgsl -lgslcblas -lm -lz -lpthread

ALL = ClassPro ClassGS prof2class class2acc class2cns RepDepAnno

GENE_LIBS = libfastk.c libfastk.h gene_core.c gene_core.h DB.c DB.h QV.c QV.h
GENE_LIBS_C = libfastk.c DB.c QV.c

GSL_LIBS = gsl/include/gsl/gsl_multifit.h
GSL_INSTALL := $(shell pwd)/gsl/

all: $(ALL)

gsl-2.7: gsl-2.7.tar.gz
	tar xzvf $<

$(GSL_LIBS): gsl-2.7
	mkdir -p gsl; cd $<; ./configure --prefix=$(GSL_INSTALL) --enable-shared=no; make; make install; cd ..

ClassPro: ClassPro.c ClassPro.h const.c io.c util.c prob.c hist.c context.c wall.c class_rel.c class_unrel.c seed.c kseq.h bessel.c bessel.h kdq.h nthash.h benchmark.c benchmark.h $(GENE_LIBS) $(GSL_LIBS)
	$(CC) $(CFLAGS) -o $@ ClassPro.c benchmark.c $(GENE_LIBS_C) $(LIBS)

ClassGS: ClassGS.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ ClassGS.c $(GENE_LIBS_C) $(LIBS)

prof2class: prof2class.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ prof2class.c $(GENE_LIBS_C) $(LIBS)

class2acc: class2acc.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ class2acc.c $(GENE_LIBS_C) $(LIBS)

class2cns: class2cns.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ class2cns.c $(GENE_LIBS_C) $(LIBS)

RepDepAnno: RepDepAnno.c $(GENE_LIBS)
	$(CC) $(CFLAGS) -o $@ RepDepAnno.c $(GENE_LIBS_C) $(LIBS)

clean:
	rm -rf gsl gsl-2.7 $(ALL)
