# make file for cpgis, etc., etc.
#
# Configure the following to suit your setup..
#
# edit these to reflect actual library and include locations..
HDRLOC = ../include/

#
# define compile and load commands here.
CFLGS = -O3
CC = gcc
CCLDFLGS =
CINCLUDES = -I../include
ZLIB = -lz
CCLIBS = -lm $(ZLIB)
LCLLIBPATH = /usr/local/lib
LCLLIB = -L $(LCLLIBPATH)
CDEFINES =
CCMP = $(CC) -c $(CFLGS) $(CDEFINES) -o $@ $(CINCLUDES) $*.c
CCLD = $(CC) $(CCLDFLGS) -o $@ $@.o $(SYSLIBS)
OBJS = bas_fns.o sqfl_fns.o sqmat_fns.o wlu_fns.o rbs_fns.o mrg_fns.o
C_OBJS = cpgis.o $(OBJS)
BAMHDRS = $(HDRLOC)bam_fns.h
BAMOBJS = bam_fns.o
FSMOBJS = $(OBJS) $(BAMOBJS) fsm_ops.o
DBPARSOBJS = $(OBJS) dbpars.o sqtrans.o
MATHSOBJS = $(OBJS) cmaths.o
MATHFSMOBJS = $(MATHSOBJS) fsm_ops.o $(BAMOBJS)

# dependency rules

EXES = rmapbsbed2cpg mkrrgenome scan_cpg_depth diffmeth identgeneloc
#EXES = rmapbsbed2cpg bin_cnts mkrrgenome diffmeth identgeneloc

all: $(EXES) cleanadaptors

nogzbuffer: $(all)
	    $(MAKE) $(all) CDEFINES="-D NO_GZBUFFER"

nozlib: $(all)
	$(MAKE) $(all) CDEFINES="-D NO_ZLIB" ZLIB="" BAMOBJS="" BAMHDRS=""

rmapbsbed2cpg: rmapbsbed2cpg.o $(OBJS)
	       $(CCLD) $(OBJS)

rmapbsbed2cpg.o: rmapbsbed2cpg.c rmapbsbed2cpg.h mrg_fns.h rbs_fns.h
		 $(CCMP)

bin_cnts: bin_cnts.o $(FSMOBJS)
	  $(CCLD) $(LCLLIB) $(FSMOBJS) $(CCLIBS)

bin_cnts.o: bin_cnts.c bin_cnts.h
	    $(CCMP)

mkrrgenome: mkrrgenome.o $(FSMOBJS)
	    $(CCLD) $(FSMOBJS) $(CCLIBS)

mkrrgenome.o: mkrrgenome.c bin_cnts.h rmapbsbed2cpg.h
	      $(CCMP)

diffmeth: diffmeth.o $(MATHFSMOBJS) $(BAMOBJS)
	  $(CCLD) $(LCLLIB) $(MATHFSMOBJS) $(CCLIBS)

diffmeth.o: diffmeth.c diffmeth.h $(BAMHDRS)
	    $(CCMP)

diffmeth_nozlib: diffmeth_nozlib.o $(MATHFSMOBJS)
	  $(CCLD) $(MATHFSMOBJS)

diffmeth_nozlib.o: diffmeth.c diffmeth.h
	    $(CC) -c $(CFLGS) -o diffmeth_nozlib.o $(CINCLUDES) diffmeth.c

scan_cpg_depth: scan_cpg_depth.o $(FSMOBJS)
		$(CCLD) $(FSMOBJS) $(CCLIBS)

scan_cpg_depth.o: scan_cpg_depth.c
		  $(CCMP)

cleanadaptors: cleanadaptors.o $(FSMOBJS)
	export LD_RUN_PATH=$(LCLLIBPATH); $(CCLD) -L$(LCLLIBPATH) $(FSMOBJS) $(CCLIBS)

cleanadaptors.o: cleanadaptors.c
	    $(CCMP)

identgeneloc: identgeneloc.o $(DBPARSOBJS)
	$(CCLD) $(DBPARSOBJS) $(CCLIBS)

identgeneloc.o: identgeneloc.c
		$(CCMP)

bas_fns.o: bas_fns.c $(HDRLOC)bas_fns.h
	$(CCMP)

sqfl_fns.o: sqfl_fns.c $(HDRLOC)sqfl_fns.h
	$(CCMP)

wlu_fns.o: wlu_fns.c $(HDRLOC)wlu_fns.h
	$(CCMP)

sqmat_fns.o: sqmat_fns.c $(HDRLOC)sqmat_fns.h
	$(CCMP)

fsm_ops.o: fsm_ops.c $(HDRLOC)fsm_ops.h
	$(CCMP)

rbs_fns.o: rbs_fns.c rbs_fns.h
	$(CCMP)

mrg_fns.o: mrg_fns.c mrg_fns.h
	$(CCMP)

cmaths.o: cmaths.c $(HDRLOC)cmaths.h
	  $(CCMP)

dbpars.o: dbpars.c $(HDRLOC)dbpars.h
	  $(CCMP)

sqtrans.o: sqtrans.c $(HDRLOC)sqtrans.h
	   $(CCMP)

bam_fns.o: bam_fns.c $(HDRLOC)bam_fns.h
	   $(CCMP)

clean:
	-rm *.o

CLEAN: clean
	-rm $(EXES)

debug:
	$(MAKE) all CFLGS=-g

mallocdebug: clean
	$(MAKE) all CFLGS="-g -DMALLOCDBG"
