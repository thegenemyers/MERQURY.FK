static char *cn_plot_script =
"#!/usr/bin/env Rscript  \n\
  \n\
require(\"argparse\")  \n\
require(\"ggplot2\")  \n\
require(\"scales\")  \n\
  \n\
parser <- ArgumentParser(description = \"Make spectra-cn plots. Line, filled, and stacked spectra-cn plots will be generated.\")  \n\
parser$add_argument(\"-f\", \"--file\", type=\"character\", help=\".spectra-cn.hist file (required)\", default=NULL)  \n\
parser$add_argument(\"-o\", \"--output\", type=\"character\", help=\"output prefix (required)\")  \n\
parser$add_argument(\"-z\", \"--zero-hist\", type=\"character\", default=\"\", help=\".only.hist file (optional, assembly only counts)\")  \n\
parser$add_argument(\"-l\", \"--cutoff\", type=\"character\", default=\"\", help=\"cutoff.txt file (optional, solid k-mer cutoffs)\")  \n\
parser$add_argument(\"-x\", \"--xdim\", type=\"double\", default=6, help=\"width of plot [default %(default)s]\")  \n\
parser$add_argument(\"-y\", \"--ydim\", type=\"double\", default=5, help=\"height of plot [default %(default)s]\")  \n\
parser$add_argument(\"-m\", \"--xmax\", type=\"integer\", default=0, help=\"maximum limit for k-mer multiplicity [default (x where y=peak) * 2.1]\")  \n\
parser$add_argument(\"-n\", \"--ymax\", type=\"integer\", default=0, help=\"maximum limit for k-mer count [default (y where y=peak) * 1.1]\")  \n\
parser$add_argument(\"-t\", \"--type\", type=\"character\", default=\"all\", help=\"available types: line, fill, stack, or all. [default %(default)s]\")  \n\
parser$add_argument(\"-p\", \"--pdf\", dest='pdf', default=FALSE, action='store_true', help=\"get output in .pdf. [default .png]\")  \n\
args <- parser$parse_args()  \n\
  \n\
gray = \"black\"  \n\
red = \"#E41A1C\"  \n\
blue = \"#377EB8\" # light blue = \"#56B4E9\"  \n\
green = \"#4DAF4A\"  \n\
purple = \"#984EA3\"  # purple = \"#CC79A7\"  \n\
orange = \"#FF7F00\"  # orange = \"#E69F00\"  \n\
yellow = \"#FFFF33\"  \n\
  \n\
merqury_col = c(gray, red, blue, green, purple, orange)  \n\
merqury_brw <- function(dat, direction=1) {  \n\
  merqury_colors=merqury_col[1:length(unique(dat))]  \n\
  if (direction == -1) {  \n\
    merqury_colors=rev(merqury_colors)  \n\
  }  \n\
  merqury_colors  \n\
}  \n\
  \n\
ALPHA=0.4  \n\
LINE_SIZE=0.3  \n\
  \n\
fancy_scientific <- function(d) {  \n\
  # turn in to character string in scientific notation  \n\
  d <- format(d, scientific = TRUE)  \n\
  # quote the part before the exponent to keep all the digits and turn the 'e+' into 10^ format  \n\
  d <- gsub(\"^(.*)e\\\\+\", \"'\\\\1'%*%10^\", d)  \n\
  # convert 0x10^00 to 0  \n\
  d <- gsub(\"\\\\'0[\\\\.0]*\\\\'(.*)\", \"'0'\", d)  \n\
  # return this as an expression  \n\
  parse(text=d)  \n\
}  \n\
  \n\
plot_zero_line <- function(zero) {  \n\
  if (!is.null(zero)) {  \n\
    if (length(zero[,1]) == 1) {  \n\
      scale_fill_manual(values = c(red), name=\"k-mer\")  \n\
    } else if (length(zero[,1]) == 2) {  \n\
      scale_fill_manual(values = c(blue, red), name=\"k-mer\")  \n\
    } else if (length(zero[,1]) == 3) {  \n\
      scale_fill_manual(values = c(purple, blue, red), name=\"k-mer\")  \n\
    } else {  \n\
      scale_fill_manual(values = merqury_brw(zero[,1]), name=\"k-mer\")  \n\
    }  \n\
  }  \n\
}  \n\
  \n\
plot_cutoff <- function(cutoff) {  \n\
  if (!is.null(cutoff)) {  \n\
    geom_vline(data = cutoff, aes(xintercept = cutoff[,2], colour = cutoff[,1]), show.legend = FALSE, linetype=\"dashed\", size=LINE_SIZE)  \n\
  }  \n\
}  \n\
  \n\
plot_zero_fill <- function(zero) {  \n\
  if (!is.null(zero)) {  \n\
    geom_bar(data=zero, aes(x=zero[,2], y=zero[,3], fill=zero[,1], colour=zero[,1], group=zero[,1]),  \n\
      position=\"stack\", stat=\"identity\", show.legend = FALSE, width = 2, alpha=ALPHA)  \n\
  }  \n\
}  \n\
  \n\
plot_zero_stack <- function(zero) {  \n\
  if (!is.null(zero)) {  \n\
    geom_bar(data=zero, aes(x=zero[,2], y=zero[,3], fill=zero[,1], colour=zero[,1]),  \n\
      position=\"stack\", stat=\"identity\", show.legend = FALSE, width = 2, alpha=ALPHA)  \n\
  }  \n\
}  \n\
  \n\
format_theme <- function() {  \n\
    theme(legend.text = element_text(size=11),  \n\
          legend.position = c(0.95,0.95),  # Modify this if the legend is covering your favorite circle  \n\
          legend.background = element_rect(size=0.1, linetype=\"solid\", colour =\"black\"),  \n\
          legend.box.just = \"right\",  \n\
          legend.justification = c(\"right\", \"top\"),  \n\
          axis.title=element_text(size=14,face=\"bold\"),  \n\
          axis.text=element_text(size=12))  \n\
}  \n\
  \n\
plot_line <- function(dat, name, x_max, y_max, zero, cutoff) {  \n\
  ggplot(data=dat, aes(x=kmer_multiplicity, y=Count, group=dat[,1], colour=dat[,1])) +  \n\
    geom_line() +  \n\
    scale_color_manual(values = merqury_brw(dat[,1]), name=\"k-mer\") +  \n\
    plot_zero_line(zero=zero) +  \n\
    plot_zero_fill(zero=zero) +  \n\
    plot_cutoff(cutoff) +  \n\
    theme_bw() +  \n\
    format_theme() +  \n\
    scale_y_continuous(labels=fancy_scientific) +  \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))  \n\
}  \n\
  \n\
plot_fill <- function(dat, name, x_max, y_max, zero, cutoff) {  \n\
  ggplot(data=dat, aes(x=kmer_multiplicity, y=Count)) +  \n\
    geom_ribbon(aes(ymin=0, ymax=pmax(Count,0), fill=dat[,1], colour=dat[,1]), alpha=ALPHA, linetype=1) +  \n\
    plot_zero_fill(zero=zero) +  \n\
    plot_cutoff(cutoff) +  \n\
    theme_bw() +  \n\
    format_theme() +  \n\
    scale_color_manual(values = merqury_brw(dat[,1]), name=\"k-mer\") +  \n\
    scale_fill_manual(values = merqury_brw(dat[,1]), name=\"k-mer\") +  \n\
    scale_y_continuous(labels=fancy_scientific) +  \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))  \n\
}  \n\
  \n\
plot_stack <- function(dat, name, x_max, y_max, zero, cutoff) {  \n\
  dat[,1]=factor(dat[,1], levels=rev(levels(dat[,1]))) #[c(4,3,2,1)] reverse the order to stack from read-only  \n\
  ggplot(data=dat, aes(x=kmer_multiplicity, y=Count, fill=dat[,1], colour=dat[,1])) +  \n\
    geom_area(size=LINE_SIZE , alpha=ALPHA) +  \n\
    plot_zero_stack(zero=zero) +  \n\
    plot_cutoff(cutoff) +  \n\
    theme_bw() +  \n\
    format_theme() +  \n\
    scale_color_manual(values = merqury_brw(dat[,1], direction=1), name=\"k-mer\", breaks=rev(levels(dat[,1]))) +  \n\
    scale_fill_manual(values = merqury_brw(dat[,1], direction=1), name=\"k-mer\", breaks=rev(levels(dat[,1]))) +  \n\
    scale_y_continuous(labels=fancy_scientific) +  \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))  \n\
}  \n\
  \n\
save_plot <- function(name, type, outformat, h, w) {  \n\
  ggsave(file = paste(name, type, outformat, sep = \".\"), height = h, width = w)  \n\
}  \n\
  \n\
spectra_cn_plot  <-  function(hist, name, zero=\"\", cutoff=\"\", w=6, h=4.5, x_max, y_max, type, pdf=FALSE) {  \n\
  # Read hist  \n\
  dat=read.table(hist, header=TRUE)  \n\
  dat[,1]=factor(dat[,1], levels=unique(dat[,1]), ordered=TRUE) # Lock in the order  \n\
    \n\
  # Read asm-only  \n\
  dat_0 = NULL  \n\
  if (zero != \"\") {  \n\
    dat_0=read.table(zero, header=FALSE)  \n\
    dat_0[,1]=as.factor(dat_0[,1])  \n\
    dat_0[,1]=factor(dat_0[,1], levels=rev(unique(dat_0[,1])), ordered=TRUE)  \n\
  }  \n\
  \n\
  # Read cutoffs  \n\
  dat_cut = NULL  \n\
  if (cutoff != \"\") {  \n\
    dat_cut=read.table(cutoff, header = FALSE)  \n\
    dat_cut[,1]=as.factor(dat_cut[,1])  \n\
    dat_cut[,1]=factor(dat_cut[,1], levels=unique(dat_cut[,1]), ordered=TRUE)  \n\
  }  \n\
  \n\
  # x and y max  \n\
  y_max_given=TRUE;  \n\
  if (y_max == 0) {  \n\
    y_max=max(dat[dat[,1]!=\"read-total\" & dat[,1]!=\"read-only\" & dat[,2] > 3,]$Count)  \n\
    y_max_given=FALSE;  \n\
  }  \n\
  if (x_max == 0) {  \n\
    x_max=dat[dat[,3]==y_max,]$kmer_multiplicity  \n\
    x_max=x_max*2.5  \n\
  }  \n\
  if (! y_max_given) {  \n\
    y_max=y_max*1.1  \n\
  }  \n\
  if (zero != \"\") {  \n\
    y_max=max(y_max, sum(dat_0[,3]*1.1))	# Check once more when dat_0 is available  \n\
  }  \n\
  \n\
  outformat=\"png\"  \n\
  if (pdf) {  \n\
    outformat=\"pdf\"  \n\
  }  \n\
    \n\
  if (type == \"line\") {  \n\
    plot_line(dat, name, x_max, y_max, zero = dat_0, cutoff = dat_cut)  \n\
    save_plot(name=name, type=\"ln\", outformat, h=h, w=w)  \n\
  }  \n\
    \n\
  else if (type == \"fill\") {  \n\
    plot_fill(dat, name, x_max, y_max, zero = dat_0, cutoff = dat_cut)  \n\
    save_plot(name=name, type=\"fl\", outformat, h=h, w=w)  \n\
  }  \n\
    \n\
  else {  # type == \"stack\")  \n\
    plot_stack(dat, name, x_max, y_max, zero = dat_0, cutoff = dat_cut)  \n\
    save_plot(name=name, type=\"st\", outformat, h=h, w=w)  \n\
  }  \n\
}  \n\
  \n\
spectra_cn_plot(hist = args$file, name = args$output, zero = args$zero, cutoff = args$cutoff, h = args$ydim, w = args$xdim, x_max = args$xmax, y_max = args$ymax, type = args$type, pdf = args$pdf)  \n\
";
