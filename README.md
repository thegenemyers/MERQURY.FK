# MerquryFK: Fast & Simple
  
<font size ="4">**_Authors:  Gene Myers & Arang Rhie_**<br>
**_First:   Feb 24, 2021_**<br>
**_Current: Feb 24, 2021_**</font>

- [Command Line](#command-line)
  - [CNplot](#CNplot)
  - [ASMplot](#ASMplot)
  - [HAPplot](#HAPplot)


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
All you do is append -P to what you've already typed and then hit return.  So for example,
```CNplot -w4 -h3 Assembly -ls Reads -P``` is acceptable input.

<a name="CNplot"></a>

```
1. CNplot [-w<double(6.0)>] [-h<double(4.5)>] [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
     [-lfs] [-P] [-z] [-o<output>] <asm>[.ktab] <reads>[.ktab]
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
but one can request .pdf's with the -P option.

The root path name for the output plots can be set with the -o option.  A suffix of .ln, .fl, or
.st is added for the line, fill, and stack plots, respectively, followed by either .png or .pdf.
If the -o option is not set, then the root path name is that of the \<asm> argument.

*Implement and describe -z and -c*


<a name="ASMplot"></a>


```
2. ASMplot [-w<double(6.0)>] [-h<double(4.5)>] [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
               [-lfs] [-P] [-z] [-o<output>] <asm1>[.ktab] [<asm2>[.ktab]] <reads>[.ktab]
```

ASMplot has the same optional parameters with the same meaning as CNplot.  What is
different is that this program looks at the spectra of the reads that (a) are in neither
asm1 or asm2, (b) in asm1 but not asm2, (c) in asm2 but not asm1, and (d) in both asm1 and
asm2.  If asm2 is missing, then it looks at the spectra of the reads that are and are not
in asm1.  The legend is appropriately labeled.

<a name="HAPplot"></a>

```
3. HAPplot [-w<double(6.0)>] [-h<double(4.5)>] [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]
               [-lfs] [-P] [-z] [-o<output>] <hap1>[.ktab] <hap2>[.ktab] <reads>[.ktab]
```

HAPplot has the same optional parameters with the same meaning as ASMplot.  *The only difference is in the legend and title ???*