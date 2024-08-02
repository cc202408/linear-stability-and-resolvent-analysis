
#  Lu Chen 2023/4/13
#
topdir = ..
libdir = $(topdir)/lib

SuperLUroot	= /usr/local
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu.a
SUPERLUHEADER  = /home/admin191/tools/superlu-6.0.0/SRC

.PHONY: z clean
.SECONDEXPANSION:

z:	ztest1 ztest1_parp ztest2 ztest3 ztestHEL_parp ztestTA_parp ztest4 ztestP ztestTA ztestTAst ztestTA_control ztest5_parp

include $(topdir)/Makefile.inc

LIBZMUMPS = -L$(libdir) -lzmumps$(PLAT) -lmumps_common$(PLAT)
	 
ztest1:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest1.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 

ztest2:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest2.o  $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 
	 
ztest3:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest3.o  $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 

ztest4:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest4.o  $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 
	 
ztestP:  $$@.o
	 $(FL) -o $@ $(OPTL) ztestP.o  $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 	 
	 
ztestTA:  $$@.o
	 $(FL) -o $@ $(OPTL) ztestTA.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)
	 	 	 
ztestTAst:  $$@.o
	 $(FL) -o $@ $(OPTL) ztestTAst.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)	 	 	 
	 	 	 
ztest1_parp:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest1_parp.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 	

ztestHEL_parp:  $$@.o
	 $(FL) -o $@ $(OPTL) ztestHEL_parp.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS) 		 	 
ztestTA_parp:  $$@.o
	 $(FL) -o $@ $(OPTL) ztestTA_parp.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

ztestTA_control:  $$@.o
	 $(FL) -o $@ $(OPTL) ztestTA_control.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

ztest5_parp:  $$@.o
	 $(FL) -o $@ $(OPTL) ztest5_parp.o c_fortran_zgssv.o st_to_ccs.o $(SUPERLULIB) $(LIBARPACK) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)
	 	  
.SUFFIXES: .o .f90 
	
.f90.o:
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.f90 $(OUTF)$*.o

clean:
	$(RM) *.o [sdcz]ztest1 [sdcz]ztest1_parp [sdcz]ztest2 [sdcz]ztest3 [sdcz]ztestHEL_parp [sdcz]ztestTA_parp [sdcz]ztest4 [sdcz]ztestP [sdcz]ztestTA [sdcz]ztestTA_control [sdcz]ztest5_parp
