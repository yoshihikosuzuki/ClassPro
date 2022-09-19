#include <sys/resource.h>
#include <time.h>
#include "gene_core.h"
#include "benchmark.h"

static struct rusage   Itime;
static struct timespec Iwall;

static struct rusage   Mtime;
static struct timespec Mwall;

void startTime()
{ getrusage(RUSAGE_SELF,&Itime);
  clock_gettime(CLOCK_MONOTONIC,&Iwall);
  Mtime = Itime;
  Mwall = Iwall;
}

void timeTo(FILE *f, int all)
{ struct rusage    now;
  struct timespec  today;
  struct rusage   *t;
  struct timespec *w;
  int usecs, umics;
  int ssecs, smics;
  int tsecs, tmics;
  int64 mem;

  getrusage(RUSAGE_SELF, &now);
  clock_gettime(CLOCK_MONOTONIC,&today);

  if (all)
    { t = &Itime;
      w = &Iwall;
      fprintf (f,"Total Resources:");
    }
  else
    { t = &Mtime;
      w = &Mwall;
      fprintf (f,"Resources for phase:");
    }


  usecs = now.ru_utime.tv_sec  - t->ru_utime.tv_sec;
  umics = now.ru_utime.tv_usec - t->ru_utime.tv_usec;
  if (umics < 0)
    { umics += 1000000;
      usecs -= 1;
    }
  if (usecs >= 60)
    fprintf (f,"  %d:%02d.%03d (m:s.ms) user",usecs/60,usecs%60,umics/1000);
  else if (usecs > 0)
    fprintf (f,"  %d.%03d (s.ms) user",usecs,umics/1000);
  else
    fprintf (f,"  %d mics user",umics);

  ssecs = now.ru_stime.tv_sec  - t->ru_stime.tv_sec;
  smics = now.ru_stime.tv_usec - t->ru_stime.tv_usec;
  if (smics < 0)
    { smics += 1000000;
      ssecs -= 1;
    }
  if (ssecs >= 60)
    fprintf (f,"  %d:%02d.%03d (m:s.ms) sys",ssecs/60,ssecs%60,smics/1000);
  else if (ssecs > 0)
    fprintf (f,"  %d.%03d (s.ms) sys",ssecs,smics/1000);
  else
    fprintf (f,"  %d mics sys",smics);

  tsecs = today.tv_sec  - w->tv_sec;
  tmics = today.tv_nsec/1000 - w->tv_nsec/1000;
  if (tmics < 0)
    { tmics += 1000000;
      tsecs -= 1;
    }
  if (tsecs >= 60)
    fprintf (f,"  %d:%02d.%03d (m:s.ms) wall",tsecs/60,tsecs%60,tmics/1000);
  else if (tsecs > 0)
    fprintf (f,"  %d.%03d (s.ms) wall",tsecs,tmics/1000);
  else
    fprintf (f,"  %d mics wall",tmics);

  fprintf(f,"  %.1f%%",(100.*(usecs+ssecs) + (umics+smics)/10000.)/(tsecs+tmics/1000000.));

  // if (all)
    { mem = now.ru_maxrss/1000000;
      fprintf(f,"  ");
      Print_Number(mem,0,f);
      fprintf(f," MB max rss");
    }

  fprintf(f,"\n");

  Mtime = now;
  Mwall = today;
}
