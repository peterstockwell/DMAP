# DMAP

	  Programs and Scripts for Bisulphite Sequence Data

Peter A. Stockwell
Dept of Pathology,
University of Otago,
Dunedin, New Zealand.

peter.stockwell@otago.ac.nz

----------------------------------------------------------------------


rmapbsbed2cpg
bin_cnts
scan_cpg_depth
mkrrgenome
cleanadaptors
diffmeth
identgeneloc

illum2fasta.awk
mk4to1lines.awk
bismmethex2list.awk
tidyrrnams.awk

These programs have been written in the course of work on differential
CpG methylation of the human genome and represent a series of tools
for preparing, modifying and analysing these data.  The work has
particularly focussed on reduced representation (RR) bisulphite
sequencing (RRBS) as in Meissner, et al., 2008, Nature, 454, 766-770,
Smith, et al., 2009) Methods, 48, 226-232 and Gu, et al., 2011, Nature
Protocols, 6, 468-481.  This software should, none-the-less, have
wider application than that.

This is research software so that it is not necessarily easy to use,
although it should work correctly as intended.  It works from a Unix
(MacOS X) or Linux command line interface: the notes below describe
the functionality and format of intermediate files where appropriate.
This code has been generated and tested on a MacOS X system (10.6 to
10.10) using gcc v4.2.1 to v4.9.3and on various Linux platforms
(RedHat, Centos, Fedora, Ubuntu) but is written to compile and run on
any appropriate C compiler and environment.  The size of files and
data required for this work will generally require a 64 bit
environment.  Awk scripts have been developed for the Gnu AWK (gawk)
distributed with MacOS X but, again, they should run in other
comparable environments.

Some of the SW is mentioned in:

Chatterjee, A., Stockwell, P.A., Rodger, E.J. and Morison, I,M.,
Comparison of alignment software for genome-wide bisulphite sequence
data, Nucleic Acids Research, doi:10.1093/nar/gks150, (2012).

and

Chatterjee, A., Rodger, E.J., Stockwell, P.A., Weeks R.J., Morison,
I.M*. “Technical considerations for reduced representation bisulfite
sequencing with multiplexed libraries” Journal of Biomedicine and
Biotechnology, (2012).

and is extensively described in

Stockwell, P.A., Chatterjee, A., Rodger, E.J. and Morison, I.M. "DMAP:
Differential Methylation Analysis Package for RRBS and WGBS data"
Bioinformatics (2014) DOI: 10.1093/bioinformatics/btu126.

----------------------------------------------------------------------
##		 Unpacking and building instructions:

Prerequisites: cleanadaptors V1.22 and later require zlib, preferably
1.2.5 or later, for data compression.  This is available from
www.zlib.net and should be built and installed in accordance with
instructions.

diffmeth 1.60 and later also require zlib in order to process bam
files.  An option exists to build diffmeth without zlib, but then
input is limited to sam files or CpG position lists.

Distribution: source code is distributed as a compressed tar archive
or as a git clone.  There are two options:

1. Downloading DMAP-main.zip from github.com/peterstockwell/DMAP will
generate a local file which can be unpacked with:

unzip DMAP-master.zip

generating a directory DMAP-master containing three directories:

include - containing some of the generic header files
data - containing contam.fa - the complete cleanadaptors data set
src - containing the rest of the source code and Makefile.

2. Using 'git clone https://github.com/peterstockwell/DMAP' which will
download the DMAP directory:

In either case the set can be built with the following commands:

cd DMAP-master/src

make

or

cd DMAP/src

make

Some warnings may appear but these can be ignored.  The executables
will be in the src directory - no install targets are provided, the
completed executables (rmapbsbed2cpg bin_cnts scan_cpg_depth
mkrrgenome cleanadaptors) should be manually copied to an appropriate
directory (e.g. /usr/local/bin).  The awk scripts are in the src
directory and may similarly be copied somewhere useful.

If the diffmeth compile fails through the lack of zlib, then this
should be obtained and installed (www.zlib.net).  If that is not
possible then you can build the program without zlib, but it will not
be able to read bam files.  To do this:

make diffmeth_nozlib

and rename the diffmeth_nozlib executable to diffmeth if required.

Documentation for the set is present in the file progs_doc.txt which
resides in the top level of the directories.

A user guide for the package along with other software is present in
RRBS_process_guideMkII.pdf
