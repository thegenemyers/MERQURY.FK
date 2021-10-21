static char *hap_plot_script =
"#!/usr/bin/env Rscript   \n\
   \n\
require(\"argparse\")   \n\
require(\"ggplot2\")   \n\
require(\"scales\")   \n\
   \n\
parser <- ArgumentParser(description = \"Make blob plots.\")   \n\
parser$add_argument(\"-f\", \"--file\", type=\"character\", help=\".count tdf; with headers; ie. <category> <seqId> <hap1Count> <hap2Count> <seqSize> (required)\", default=NULL)   \n\
parser$add_argument(\"-o\", \"--output\", type=\"character\", default=\"hapmers.blob\", help=\"output prefix [default %(default)s]\")   \n\
parser$add_argument(\"-x\", \"--xdim\", type=\"double\", default=6.5, help=\"width of output plot [default %(default)s]\")   \n\
parser$add_argument(\"-y\", \"--ydim\", type=\"double\", default=6, help=\"height of output plot [default %(default)s]\")   \n\
parser$add_argument(\"-p\", \"--pdf\", dest='pdf', default=FALSE, action='store_true', help=\"set to get output in .pdf. [default .png]\")   \n\
args <- parser$parse_args()   \n\
   \n\
blob_plot <- function(dat=NULL, out, w=6.5, h=6, pdf=FALSE) {   \n\
   \n\
  dat=read.table(dat, header=TRUE, comment.char=\"\")   \n\
   \n\
  max_total=max(max(dat[,3]), max(dat[,4])) * 1.01   \n\
  col_lab=names(dat)[1]   \n\
  x_lab=names(dat)[3]   \n\
  y_lab=names(dat)[4]   \n\
   \n\
  ggplot(dat, aes(x=dat[,3], y=dat[,4], fill=dat[,1], size=dat[,5])) +   \n\
    geom_point(shape=21, alpha=0.3) + theme_bw() +   \n\
    scale_color_brewer(palette = \"Set1\") +   \n\
    scale_fill_brewer(palette = \"Set1\") +     # Set1 -> Red + Blue. Set2 -> Green + Orange.   \n\
    scale_x_continuous(labels=comma, limits = c(0, max_total)) +   \n\
    scale_y_continuous(labels=comma, limits = c(0, max_total)) +   \n\
    scale_size_continuous(labels=comma, range = c(1, 10), name = \"Total k-mers (size)\") +   \n\
    theme(legend.text = element_text(size=11),   \n\
          legend.position = c(0.95,0.95),  # Modify this if the legend is covering your favorite circle   \n\
          legend.background = element_rect(size=0.5, linetype=\"solid\", colour =\"black\"),   \n\
          legend.box.just = \"right\",   \n\
          legend.justification = c(\"right\", \"top\"),   \n\
          axis.title=element_text(size=14,face=\"bold\"),   \n\
          axis.text=element_text(size=12)) +   \n\
    guides( size = guide_legend(order = 1),   \n\
            fill = guide_legend(override.aes = list(size=5, alpha=1), order = 2, title = col_lab)) +   \n\
    xlab(x_lab) + ylab(y_lab)   \n\
   \n\
  outformat=\"png\"   \n\
  if (pdf) {   \n\
    outformat=\"pdf\"   \n\
  }     \n\
  ggsave(file = paste(out, outformat, sep=\".\"), height = h, width = w)   \n\
   \n\
}   \n\
   \n\
blob_plot(dat = args$file, out = args$output, w = args$xdim, h = args$ydim, pdf = args$pdf)   \n\
";
