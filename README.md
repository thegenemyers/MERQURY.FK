# MerquryFK: Fast & Simple
  
<font size ="4">**_Authors:  Gene Myers & Arang Rhie_**<br>
**_First:   Feb 24, 2021_**<br>
**_Current: Feb 24, 2021_**</font>

- [Command Line](#command-line)
  - [CNplot](#CNplot)
  - [ASMplot](#ASMplot)
  - [CNspectra](#CNspectra)
  - [KatComp](#KatComp)
  - [KatGC](#KatGC)


## Command Line

The original **Merqury** is a collection of R and shell scripts for producing k-mer analysis plots of genomic sequence data and assemblies with **meryl** as its core k-mer counter infra-structure.
**MerquryFK** replaces meryl with the **FastK** k-mer counter suite to considerably speed up analyses.
Moreover, all the R and shell scripts have been refactored into a typical collection of UNIX command line tools that the user will hopefully experience as easier to comprehend and invoke.

There are some general conventions for our tools programmed for your convenience.
First, suffix extensions need not be given for arguments of known type.  For example,
if an argument is a fasta or fastq with root name "foo" without extensions, then
our commands will look for ```foo.fasta, foo.fa, foo.fastq, and foo.fq``` if you specify
```foo``` as the argument.  Second, option arguments (those that begin with a '-') can
be in any order and in any position relaltive to the non-optional primary arguments (which must
be given in the order specified).  We find this pretty convenient when for example you
have typed out an entire CNplot command (1. below) but forgot that you wanted .pdf's.
All you do is append -pdf to what you've already typed and then hit return.  So for example,
```CNplot -w4 -h3 Assembly -ls Reads -pdf``` is acceptable input.

<a name="CNplot"></a>

```
1. CNplot [-w<double(6.0)>] [-h<double(4.5)>]
          [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
          [-lfs] [-pdf] [-z] [-T<int(4)>]
          [-o<output>] <asm>[.ktab] <reads>[.ktab]
```

Given k-mer tables, produced by FastK, for an assembly, \<asm>,
and a read data set, \<reads>, of the same genome, *CNplot* produces
copy-number spectrum plots for the pair.

The width and height in inches of the plots are controlled by the -w and -h options and
by default they are 6 x 4.5 inches.

The maximum x-coordinate (k-mer frequency) and y-coordinate (k-mer count) can be set
in either an absolute or relative fashion.  If absolute one specifes the maxima as integer
arguments to the -X and -Y options.  CN spectrum typically have peaks away from the origin.
The x and y coordinates of the highest of these peaks, say x\*,y\* provide a relative
landmark for setting the axis limits and with the -x and -y options one can set the maximums
as a multiple of x\* or y\*.  By default, CNplot sets the axis maximums to x\*&#183;2.1 and y\*&#183;1.1.

CNplot can produce any of a line plot, a filled plot, or a so-called stacked plot where the individual histograms accumulate.  By default it produces all three, but the user can select any subset with the options -l (line), -f (fill), and -s (stack).  By default the plots are .png's
but one can request .pdf's with the -pdf option.

The root path name for the output plots can be set with the -o option.  A suffix of .ln, .fl, or
.st is added for the line, fill, and stack plots, respectively, followed by either .png or .pdf.
If the -o option is not set, then the root path name is that of the \<asm> argument.

If the -z option is set, then CNplot plots at 0, the # of k-mers in \<asm> - \<reads> broken down
into those that are unique or those that are not.

The -T option controls the number of threads used by FastK's "Logex" which is the dominant
computational cost for CNplot.

*Implement and describe -c ?*


<a name="ASMplot"></a>


```
2. ASMplot [-w<double(6.0)>] [-h<double(4.5)>]
           [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
           [-lfs] [-pdf] [-z] [-T<int(4)>]
           [-o<output>] <asm1>[.ktab] [<asm2>[.ktab]] <reads>[.ktab]
```

ASMplot has the same optional parameters with the same meaning as CNplot.  What is
different is that this program looks at the spectra of the reads that (a) are in neither
asm1 or asm2, (b) in asm1 but not asm2, (c) in asm2 but not asm1, and (d) in both asm1 and
asm2.  If asm2 is missing, then it looks at the spectra of the reads that are and are not
in asm1.  The legend is appropriately labeled.

<a name="CNspectra"></a>

```
3. CNspectra [-v] [-lfs] [-pdf] [-T<int(4)>] <read> <asm1> [<asm2>] <out>
```

CNspectra produces copy-number spectra plots for each assembly, a copy-number spectra plot of the union of both assemblies (if two are present), an assembly spectra plot, and tables of qv and completeness statistics, as well as a .bed file of potential error locations in each supplied assembly.

The primary input arguments -- \<read>, \<asm1>, and \<asm2> (if present) -- are expected to be the root path names of FastK tables, histograms, and profiles.  Specifically, CNspectra expects
to find:

* **\<read>.hist & .ktab** produced by <code>FastK -t1 -k\<K> [...] \<read_data> -N\<read></code>

* **\<asm>.ktab & .prof** produced by <code>FastK -t1 -p -k\<K> [...] \<assembly> -N\<asm></code>

* **\<asm>.\<read>.prof** produced by <code>FastK -p:<read> -k\<K> [...] \<assembly> -N\<asm></code>

where the k-mer size \<K> is the same for all calls.

The primary output argument -- \<out> -- is the root path name for all the output files
produced by CNspectra.  Specifically, it *can* produce:

* **\<out>.\<asm>.spectra-cn.***: cn-spectra plots of \<asm>

* **\<out>.\<asm>.qv**: error and qv table for each scaffold of \<asm>

* **\<out>.\<asm>_only.bed**: a .bed file of the locations where the assembly has k-mer's not supported by the read data set.

* **\<out>.spectra-cn.***: cn-spectra plots of the union of \<asm1> and \<asm2>

* **\<out>.spectra-asm.***: assembly spectra plots of the assemblies

* **\<out>.qv**: error and qv of each assembly as a whole

* **\<out>.completeness-stats**: coverage of solid read k-mers by the assemblies and their union (if two are given).

One can select verbose output with -v, .pdf plots versus .png's with -pdf, and which
type of plots -- line, fill, or stacked -- with a combination of the flags -lfs.
If no plot types are set, then all 3 are produced.  Finally, the -T option controls
the number of threads used in those bits of CNspectra that are threaded.

CNspectra uses the default dimensions and scaling parameters of CNplot and ASMplot when
producing plots and always with the -z option set.  These settings can be reset by redefining
a easily identifiable set of defined constants at the top of CNspectra.c


<a name="KatComp"></a>

```
4. KatComp [-w<double(6.0)>] [-h<double(4.5)>]
           [-[xX]<number(x2.1)>] [-[yY]<number(y2.1)>]
           [-lfs] [-pdf] [-T<int(4)>]
           [-o<output>] <source1[.ktab] <source2>[.ktab]
```

Given two FastK k-mer tables, `<source1>` and `<source2>`, with the same k-mer size, *KatComp* produces a 3D heat map or contour map of the product of the two k-mer
spectra.  The controlling options are almost identical to those of CNplot, save that
(1) the -z option is not relevant, and (2) the meaning of the plot type options,
-l, -f, and -s -- are as follows:

The -l option produces a contour **line** plot of count iso-lines.   The -f option produces
a **filled** heat map of the counts.  The -s option produces a heap map with a contour plot
**stacked** on top of it.

Another difference with CNplot, is that the y-axis denotes the frequency of k-mers in the second data set, rather than the count of k-mers in the lone data set.

<a name="KatGC"></a>

```
5. KatGC [-w<double(6.0)>] [-h<double(4.5)>]
         [-[xX]<number(x2.1)>] [-lfs] [-pdf] [-T<int(4)>]
         [-o<output>] <source>[.ktab]
```

Given a k-mer table, `<source>`, produced by FastK, *KatGC* produces a 3D heat map
or contour map of the frequency of a k-mer versus its' GC content.
The controlling options are almost identical to those of CNplot, save that
(1) -z and -[yY] are not relevant, and (2) the meaning of the plot type options,
-l, -f, and -s -- are as follows:

The -l option produces a contour **line** plot of count iso-lines.   The -f option produces
a **filled** heat map of the counts.  The -s option produces a heap map with a contour plot
**stacked** on top of it.