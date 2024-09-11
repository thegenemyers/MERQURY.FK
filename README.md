# MerquryFK & KatFK: Fast & Simple
  
<font size ="4">**_Authors:  Gene Myers & Arang Rhie_**<br>
**_First:   Feb 24, 2021_**<br>
**_Current: Aug 11, 2021_**</font>

- [Introduction](#introduction)
  - [HAPmaker](#HAPmaker)
  - [CNplot](#CNplot)
  - [ASMplot](#ASMplot)
  - [HAPplot](#HAPplot)
  - [MerquryFK](#MerquryFK)
  - [KatComp](#KatComp)
  - [KatGC](#KatGC)
  - [PloidyPlot](#PloidyPlot)

## Introduction

The original [Merqury](https://github.com/marbl/merqury) is a collection of R, Java, and shell scripts for producing k-mer analysis plots of genomic sequence data and assemblies with **meryl** as its core k-mer counter infra-structure.
**MerquryFK** replaces meryl with the **FastK** k-mer counter suite to considerably speed up analyses.
Moreover, all the R, Java, and shell scripts have been refactored into a typical collection of UNIX command line tools that the user will hopefully experience as easier to comprehend and invoke.  It does still require that R be installed and that the R-packages ```argparse```, ```ggplot2```, ```scales```, ```viridis```, and ```cowplot``` are present.  In addition, we have realized some analyses, **KatComp** and **KatGC**, that one finds
only in the somewhat similar [KAT](https://github.com/TGAC/KAT) k-mer suite developed at the Earlham Institute.
Lastly, we include in this collection, **PloidyPlot** which is an improved version of the
ploidy plotting tool [SmudgePlot](https://github.com/KamilSJaron/smudgeplot).

There are some general conventions for our tools programmed for your convenience.
First, suffix extensions need not be given for arguments of a known type.  For example,
if an argument is a fasta or fastq with root name "foo" without extensions, then
our commands will look for ```foo.fasta, foo.fa, foo.fastq, and foo.fq``` if you specify
```foo``` as the argument.  Second, option arguments (those that begin with a '-') can
be in any order and in any position relative to the non-optional primary arguments (which must
be given in the order specified).  We find this pretty convenient when for example you
have typed out an entire CNplot command (2. below) but forgot that you wanted .pdf's.
All you do is append -pdf to what you've already typed and then hit return.  So for example,
```CNplot -w4 -h3 Assembly -ls Reads -pdf``` is acceptable input.

For the tools that take a FastK k-mer table as an input, we use the syntax \<name>[.ktab],
to describe it on the command line indicating that the .ktab extension is optional as
per the convention above.  Regardless of whether the extension is given, it is expected
that the associated histogram file \<name>.hist is also present (this file is always
produced by a run of FastK that produces a k-mer table).  Also note carefully that these
tables must be produced with the option -t or -t1 set so that all k-mers that occur 1 or more
times in the underlying data set are in the table.

<br>

<a name="HAPmaker"></a>

```
1. HAPmaker [-v] [-T<int(4).] <mat>[.ktab] <pat>[.ktab] <child>[.ktab]
```

Prior to running either HAPplot or MerquryFK in trio mode, one must produce a table
of the **hap-mers** for both the mother and father given FastK k-mer tables of the
maternal, paternal, and child sequence data sets.
The hap-mer k-mers are those in the given parent that are specific to that parent, inherited
by the child, and reliable in that they are unlikely to be error k-mers.
The names of the resulting k-mer tables produced are \<mat>.hap.ktab and \<pat>.hap.ktab.

The -T option can be used to control the number of threads used, and -v turns on verbose
reporting to the standard error.

<br>

<a name="CNplot"></a>

```
2. CNplot [-w<double(6.0)>] [-h<double(4.5)>]
          [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
          [-vk] [-lfs] [-pdf] [-z] [-T<int(4)>] [-P<dir(/tmp)>]
          <reads>[.ktab] <asm:dna> <out>[.cni]
or
          
   CNplot [-w<double(6.0)>] [-h<double(4.5)>]
          [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
          [-v] [-lfs] [-pdf] [-z] <out>[.cni]
```

Given a k-mer table, produced by FastK, for
a read data set, \<reads>,
and an assembly, \<asm>, of the same genome, *CNplot* produces
copy-number spectrum plots for the pair.
The type 'dna' of \<asm> is any dna file format accepted by [FastK](https://github.com/thegenemyers/FASTK) (e.g. .fasta, .cram, .fa.gz, ...)

The width and height in inches of the plots are controlled by the -w and -h options and
by default they are 6 x 4.5 inches.

The maximum x-coordinate (k-mer frequency) and y-coordinate (k-mer count) can be set
in either an absolute or relative fashion.  If absolute one specifes the maxima as integer
arguments to the -X and -Y options.  CN spectra typically have peaks away from the origin.
The x and y coordinates of the highest of these peaks, say x\*,y\* provide a relative
landmark for setting the axis limits and with the -x and -y options one can set the maximums
as a multiple of x\* or y\*.  By default, CNplot sets the axis maximums to x\*&#183;2.1 and y\*&#183;1.1.

CNplot can produce any of a line plot, a filled plot, or a so-called stacked plot where the individual histograms accumulate.  By default it produces all three, but the user can select any subset with the options -l (line), -f (fill), and -s (stack).  By default the plots are .png's
but one can request .pdf's with the -pdf option.

The root path name for the output plots is given by the final \<out> argument.  A suffix of .ln, .fl, or .st is added for the line, fill, and stack plots, respectively, followed
by either .png or .pdf.

If the -z option is set, then CNplot plots at 0, the # of k-mers in \<asm> - \<reads>
broken down into those that are unique or those that are not.

The -T option can be used to control the number of threads used, -v turns on verbose
reporting to the standard error, and -P is passed through to the calls to FastK so you can
specify the temp directory if needed.

If the -k option is set then in addition to producing the requested .png or .pdf files with
root path \<out\>, CNplot also creates \<out\>.cni that contains the processed data that is given to an R-script to produce the plots.  This allows
users to call CNplot again, but with just \<out\>[.cni] as the primary argument, and the
plots will be produced again under control of the parameters given.  This saves having to
wait for the plotting data to have to be computed again, and permits one to repeatedly
plot until they are satisfied with the size, scaling, and form.

The use is responsible for deleting \<out\>.cni when desired.

<br>

<a name="ASMplot"></a>

```
3. ASMplot [-w<double(6.0)>] [-h<double(4.5)>]
           [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
           [-vk] [-lfs] [-pdf] [-z] [-T<int(4)>] [-P<dir(/tmp)>]
           <reads>[.ktab] <asm1:dna> [<asm2:dna>] <out>[.asmi]
or

   ASMplot [-w<double(6.0)>] [-h<double(4.5)>]
           [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
           [-v] [-lfs] [-pdf] [-z] <out>[.asmi]
```

ASMplot has the same optional parameters with the same meaning as CNplot.  What is
different is that this program looks at the spectra of the k-mers that (a) are in neither
asm1 or asm2, (b) in asm1 but not asm2, (c) in asm2 but not asm1, and (d) in both asm1 and
asm2.  If asm2 is missing, then it looks at the spectra of the read k-mers that are and
are not in asm1.  The legend is appropriately labeled.

The -k option works similarly to CNplot save that the short-cut file is \<out\>.asmi.

<br>

<a name="HAPplot"></a>

```
4. HAPplot [-vk] [-w<double(6.0)>] [-h<double(4.5)>] [-pdf] [-T<int(4)>] [-P<dir(/tmp)>]
           <mat>[.hap[.ktab]] <pat>[.hap[.ktab] <asm1:dna> [<asm2:dna>] <out>[.hpi]
or
           
   HAPplot [-v] [-w<double(6.0)>] [-h<double(4.5)>] [-pdf] <out>[.hpi]
```

HAPplot has the relevant optional parameters of CNplot with the same meaning.
It produces a haplotype blob plot of assembly \<asm1> and if present \<asm2> given 
hap-mer tables \<mat>.hap.ktab, and \<pat>.hap.ktab produced earlier by [HAPmaker](#HAPmaker).
Each assembly contig is plotted as a blob where its size is porportional
to the contig's length in bases, and its' position is (x,y) where x = # of maternal hap-mers,
and y = # of paternal hap-mers, in the contig.

The -k option works similarly to CNplot save that the short-cut file is \<out\>.hpi.

<br>

<a name="MerquryFK"></a>

```
5. MerquryFK [-w<double(6.0)>] [-h<double(4.5)>]
             [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
             [-vk] [-lfs] [-pdf] [-z] [-T<int(4)>] [-P<dir(/tmp)>]
                <read>[.ktab] [ <mat>[.hap[.ktab]] <pat>[.hap[.ktab]] ]
                <asm1:dna> [<asm2:dna>] <out>
```

MerquryFK runs all the analyses performed by the original Merqury suite where it will
run the trio analyses if hap-mer tables (produced by [HAPmaker](#HAPmaker) above) are supplied for the mother and father read data sets, and will assume a single unphased assembly if only \<asm1> is given, or two phased, haplotype assemblies if \<asm2> is also given.  The assemblies are assumed to be in a dna-sequence file format acceptable to [FastK](https://github.com/thegenemyers/FASTK) (i.e. fasta or fastq, compressed or not, and cram, bam, sam, or a Dazzler DB).

The options on the first 3 lines of the command description are the same as for CNplot
and have the same interpretation.  Any sub-calls to CNplot, ASMplot, and HAPplot will use
the option setting if not the defaults.  So note that if the -k option is set then various
.cni, .asmi, and .hpi files will be generated with root names as descripted below.

The primary output argument -- \<out> -- is the root path name for all the output files
produced by MerquryFK.  Specifically, it will produce the following where \<asm> is the
root name of an input assembly file:

* **\<out>.\<asm>.spectra-cn.(ln+fl+st).(pdf | png)**: cn-spectra plots of \<asm> versus \<reads>

* **\<out>.\<asm>.qv**: error and qv table for each scaffold of \<asm>

* **\<out>.\<asm>_only.bed**: a .bed file of the locations where \<asm> has k-mer's not supported by those in \<read>.

* **\<out>.spectra-cn.(ln+fl+st).(pdf | png)**: cn-spectra plots of \<asm1> U \<asm2> versus \<reads> (if 2 assemblies)

* **\<out>.spectra-asm.(ln+fl+st).(pdf | png)**: assembly spectra plots of the assemblies versus \<reads>

* **\<out>.qv**: error and qv of each assembly as a whole

* **\<out>.completeness.stats**: coverage of solid read k-mers by the assemblies and their union (if two are given).

When run in trio mode it further outputs the following files where \<hap> is either \<mat> or \<pat>:

* **\<out>.\<asm>.\<hap>.spectra-cn.(ln+fl+st).(pdf | png)***: cn-spectra plots of \<asm> versus \<hap>

* **\<out>.hapmers.blob.(pdf | png)**: haplotype phased blob plot of the assemblys' contigs colored by assembly

* **\<out>.\<asm>.phased_block.bed**: putative phased intervals of an assembly determined with the hapmer tables.

* **\<out>.\<asm>.phased_block.blob.(pdf | png)**: phased blob plot of an assembly with the putative blocks colored by phase

* **\<out>.\<asm>.phased_block.stats**: short summary of the block partitioning of \<asm>

* **\<out>.\<asm>.block.N.(pdf | png)**: length histogram of the phase-colored blocks of \<asm> in order of largest to smallest scaled as a function of N# stat.

* **\<out>.\<asm>.continuity.N.(pdf | png)**: N-plots of haplotype blocks and contigs and if present scaffolds.

<br>

<a name="KatComp"></a>

```
6. KatComp [-w<double(6.0)>] [-h<double(4.5)>]
           [-[xX]<number(x2.1)>] [-[yY]<number(y2.1)>]
           [-lfs] [-pdf] [-T<int(4)>]
           <source1[.ktab] <source2>[.ktab] <out>
```

Given two FastK k-mer tables, `<source1>` and `<source2>`, with the same k-mer size, *KatComp* produces a 3D heat map or contour map of the product of the two k-mer
spectra.  The controlling options are almost identical to those of CNplot, save that
(1) the -z option is not relevant, and (2) the meaning of the plot type options,
-l, -f, and -s -- are as follows:

The -l option produces a contour **line** plot of count iso-lines.   The -f option produces
a **filled** heat map of the counts.  The -s option produces a heat map with a contour plot
**stacked** on top of it.

Another difference with CNplot, is that the y-axis denotes the frequency of k-mers in the second data set, rather than the count of k-mers in the lone data set.
Also the -k option is not provided as the compute time prefacing the plot is reasonable.

<br>

<a name="KatGC"></a>

```
7. KatGC [-w<double(6.0)>] [-h<double(4.5)>]
         [-[xX]<number(x2.1)>] [-lfs] [-pdf] [-T<int(4)>]
         <source>[.ktab] <out>
```

Given a k-mer table, `<source>` produced by FastK, *KatGC* produces a 3D heat map
or contour map of the frequency of a k-mer versus its' GC content.
The controlling options are almost identical to those of CNplot, save that
(1) -z and -[yY] are not relevant, and (2) the meaning of the plot type options,
-l, -f, and -s -- are as follows:

The -l option produces a contour **line** plot of count iso-lines.   The -f option produces
a **filled** heat map of the counts.  The -s option produces a heat map with a contour plot
**stacked** on top of it.

Another difference with CNplot, is that the y-axis denotes the % GC content of k-mers, rather than the count of k-mers in the lone data set.
Also the -k option is not provided as the compute time prefacing the plot is reasonable.

<br>

<a name="PloidyPlot"></a>

```
8. PloidyPlot [-w<double(6.0)>] [-h<double(4.5)>]
              [-vk] [-lfs] [-pdf] [-T<int(4)>] [-P<dir(/tmp)>]
              <source>[.ktab] <out>[.smu]
or

   PloidyPlot [-w<double(6.0)>] [-h<double(4.5)>]
              [-vk] [-lfs] [-pdf] <out>[.smu]
```

This is an improved version of [SmudgePlot](https://github.com/KamilSJaron/smudgeplot)
that produces "cleaner" smudges by avoiding false het-mer signals.  In addition to displaying smudges the program tries to estimate the ploidy of the genome albeit this
is not guaranteed to always be correct.

Almost all the options (-w, -h, -v, -lfs, -pdf, -T) control the output file name and type, display, and number of threads used exactly as described for [CNplot](#CNplot).

Any k-mer with a count of less than -e in the input FastK table \<source> is considered
an error in the analysis.  The analysis is ultimately run over a *symmetric** k-mer table of all k-mers with count 4 or more.  If the supplied table does not meet these specifications,
then the program takes additional compute time to make it so, but if in a preprocessing
step, you use [Logex](https://github.com/thegenemyers/FASTK/#Logex) and [Symmex](https://github.com/thegenemyers/FASTK/#Symmex) to make the table conform to the internal
requirements than this time is saved.  If PloidyPlot is forced to make a trimmed, symmetric table for an input, then the call to Symmex will use the temp directory specified by the -P option (if not the default "/tmp").

Even if the input table is symmetric and trimmed, the bulk of the
time taken by PloidyPlot is in accumulating count statistics of het-mer pairs.  If the
-k option is given then the table of het-mer pair statistics is **k**ept, and stored
in a file named `<out>.smu`.  This saved table can then be given as the sole required
argument to PloidyPlot in subsequent calls, so that the time intensive counting step is avoided.