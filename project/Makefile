
#  Lu Chen 2023/4/13
#
topdir = ..
libdir = $(topdir)/lib

.PHONY: z clean
.SECONDEXPANSION:

z:	ztest1

include $(topdir)/Makefile.inc

LIBZMUMPS = -L$(libdir) -lzmumps$(PLAT) -lmumps_common$(PLAT)
	 
ztest1:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest1.o $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 


.SUFFIXES: .o .f90
	
.f90.o:
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.f90 $(OUTF)$*.o

clean:
	$(RM) *.o [sdcz]ztest1
