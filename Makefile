topdir = .
libdir = $(topdir)/lib

default: d

.PHONY: default clean
.PHONY: all s d c z prerequisites libseqneeded
.PHONY: allshared sshared dshared cshared zshared prerequisitesshared libseqneededsharedlibseq sharedlibseq

all: prerequisites
	cd src; $(MAKE) all
	cd examples; $(MAKE) all

s: prerequisites
	cd src; $(MAKE) s
	cd examples; $(MAKE) s

d: prerequisites
	cd src; $(MAKE) d
	cd examples; $(MAKE) d

c: prerequisites
	cd src; $(MAKE) c
	cd examples; $(MAKE) c

z: prerequisites
	cd src; $(MAKE) z
	cd examples; $(MAKE) z


allshared: prerequisitesshared
	cd src; $(MAKE) allshared
	cd examples; $(MAKE) all

sshared: prerequisitesshared
	cd src; $(MAKE) sshared
	cd examples; $(MAKE) s

dshared: prerequisitesshared
	cd src; $(MAKE) dshared
	cd examples; $(MAKE) d

cshared: prerequisitesshared
	cd src; $(MAKE) cshared
	cd examples; $(MAKE) c

zshared: prerequisitesshared
	cd src; $(MAKE) zshared
	cd examples; $(MAKE) z


# Is Makefile.inc available ?
Makefile.inc:
	@echo "######################################################################"
	@echo "# BEFORE COMPILING MUMPS, YOU MUST HAVE AN APPROPRIATE Makefile.inc"
	@echo "# FILE AVAILABLE. PLEASE CHECK THE DIRECTORY ./Make.inc FOR EXAMPLES"
	@echo "# OF Makefile.inc FILES, AND USE Make.inc/Makefile.inc.generic IF YOU"
	@echo "# NEED TO BUILD A NEW ONE. SEE ALSO THE README AND INSTALL FILES"
	@echo "######################################################################"
	@exit 1

include Makefile.inc

prerequisites: Makefile.inc $(LIBSEQNEEDED) $(libdir)/libpord$(PLAT)$(LIBEXT)

prerequisitesshared: Makefile.inc $(LIBSEQNEEDED)sharedlibseq $(libdir)/libpord$(PLAT)$(LIBEXT_SHARED)

# Rules for fake MPI library used to avoid using MPI:
#
# If $(LIBSEQNEEDED) is empty, prerequisitesshared includes a dependenecy on
# the sharedlibseq suffix dependency which we always satisfy

sharedlibseq:

libseqneeded:
	(cd libseq; $(MAKE))
	(cp libseq/lib* $(libdir))

libseqneededsharedlibseq:
	(cd libseq; $(MAKE) sharedlibmpiseq)
	(cp libseq/lib* $(libdir))

# Build the libpord.a library and copy it into $(topdir)/lib
$(libdir)/libpord$(PLAT)$(LIBEXT):
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); \
	  $(MAKE) CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" RANLIB="$(RANLIB)" OUTC="$(OUTC)" LIBEXT="$(LIBEXT)" LIBEXT_SHARED="$(LIBEXT_SHARED)" libpord$(LIBEXT); \
	fi;
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cp $(LPORDDIR)/libpord$(LIBEXT) $@; \
	fi;

$(libdir)/libpord$(PLAT)$(LIBEXT_SHARED):
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); \
	  $(MAKE) PLAT="$(PLAT)" FPIC="$(FPIC_OPT)" CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" RANLIB="$(RANLIB)" OUTC="$(OUTC)" LIBEXT="$(LIBEXT)" LIBEXT_SHARED="$(LIBEXT_SHARED)" libpord$(LIBEXT_SHARED); \
	fi;
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cp $(LPORDDIR)/libpord$(LIBEXT_SHARED) $@; \
	fi;




clean:
	(cd src; $(MAKE) clean)
	(cd examples; $(MAKE) clean)
	(cd $(libdir); $(RM) lib*$(PLAT)$(LIBEXT) lib*$(PLAT)$(LIBEXT_SHARED))
	(cd libseq; $(MAKE) clean)
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); $(MAKE) CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" RANLIB="$(RANLIB)" OUTC="$(OUTC)" LIBEXT="$(LIBEXT)" LIBEXT_SHARED="$(LIBEXT_SHARED)" PLAT="$(PLAT)" realclean; \
        fi;


