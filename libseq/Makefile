#
#  This file is part of MUMPS 5.5.1, released
#  on Tue Jul 12 13:17:24 UTC 2022
#
all: staticlibmpiseq

.PHONY: all staticlibmpiseq sharedlibmpiseq clean

topdir = ..
include $(topdir)/Makefile.inc

staticlibmpiseq: libmpiseq$(PLAT)$(LIBEXT)

libmpiseq$(PLAT)$(LIBEXT): mpi.o mpic.o elapse.o
	$(AR)$@ mpi.o mpic.o elapse.o
	$(RANLIB) $@

sharedlibmpiseq:
	$(MAKE) FPIC="$(FPIC_OPT)" libmpiseq$(PLAT)$(LIBEXT_SHARED)

libmpiseq$(PLAT)$(LIBEXT_SHARED): mpi.o mpic.o elapse.o
	$(FC) -shared mpi.o mpic.o elapse.o -o libmpiseq$(PLAT)$(LIBEXT_SHARED)

.f.o:
	$(FC) $(OPTF) $(FPIC)              -c $*.f $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(FPIC) $(CDEFS) -I. -c $*.c $(OUTC)$*.o

clean:
	$(RM) *.o *$(LIBEXT)
