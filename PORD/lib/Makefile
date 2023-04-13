
# To compile directly, uncomment the line below.
# include ../Make.in
# # include $(BUILDDIR)/Makefile.inc
#
# Otherwise, adequate variables for CC, CFLAGS, AR and
# RANLIB must be passed to make.
#

INCLUDES = -I../include

COPTIONS = $(INCLUDES) $(CFLAGS) $(OPTFLAGS) $(FPIC)

OBJS = graph.o gbipart.o gbisect.o ddcreate.o ddbisect.o nestdiss.o \
       multisector.o gelim.o bucket.o tree.o \
       symbfac.o interface.o sort.o minpriority.o 

# Note: numfac.c read.c mapping.c triangular.c matrix.c kernel.c
# were not direcly used by MUMPS and have been removed from the
# original SPACE package.

# OUTC = -o

.c.o:
	$(CC) $(COPTIONS) -c $*.c $(OUTC)$*.o

libpord$(PLAT)$(LIBEXT):$(OBJS)
	$(AR)$@ $(OBJS)
	$(RANLIB) $@

libpord$(PLAT)$(LIBEXT_SHARED):$(OBJS)
	$(CC) -shared $(OBJS) -o libpord.so

clean:
	rm -f *.o

realclean:
	rm -f *.o libpord$(PLAT)$(LIBEXT) libpord$(PLAT)$(LIBEXT_SHARED)
