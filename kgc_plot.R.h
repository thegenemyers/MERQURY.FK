static char *kgc_plot =
"#!/usr/bin/env Rscript  \n\
  \n\
require(\"argparse\")  \n\
require(\"ggplot2\")  \n\
require(\"scales\")  \n\
require(\"viridis\")  \n\
  \n\
parser <- ArgumentParser(description = \"Make spectra-cn plots. Line, filled, and stacked spectra-cn plots will be generated.\")  \n\
parser$add_argument(\"-f\", \"--file\", type=\"character\", help=\".spectra-cn.hist file (required)\", default=NULL)  \n\
parser$add_argument(\"-o\", \"--output\", type=\"character\", help=\"output prefix (required)\")  \n\
parser$add_argument(\"-x\", \"--xdim\", type=\"double\", default=6, help=\"width of plot [default %(default)s]\")  \n\
parser$add_argument(\"-y\", \"--ydim\", type=\"double\", default=5, help=\"height of plot [default %(default)s]\")  \n\
parser$add_argument(\"-t\", \"--type\", type=\"character\", default=\"all\", help=\"available types: line, fill, stack, or all. [default %(default)s]\")  \n\
parser$add_argument(\"-p\", \"--pdf\", dest=\'pdf\', default=FALSE, action=\'store_true\', help=\"get output in .pdf. [default .png]\")  \n\
parser$add_argument(\"-s\", \"--source\", type=\"character\", help=\"source .ktab file (required)\", default=NULL)  \n\
args <- parser$parse_args()  \n\
  \n\
fancy_scientific <- function(d) {  \n\
  if (d[2] > 10000000) {  \n\
    for (i in 1:length(d)) {  \n\
      if (is.na(d[i])) {  \n\
        next  \n\
      }  \n\
      d[i] <- paste( as.character(as.integer(d[i])/1000000), \"M\", sep=\"\")  \n\
    }  \n\
  } else if (d[2] > 10000) {  \n\
    for (i in 1:length(d)) {  \n\
      if (is.na(d[i])) {  \n\
        next  \n\
      }  \n\
      d[i] <- paste( as.character(as.integer(d[i])/1000), \"K\", sep=\"\")  \n\
    }  \n\
  } else {  \n\
    for (i in 1:length(d)) {  \n\
      if (is.na(d[i])) {  \n\
        next  \n\
      }  \n\
      d[i] <- as.character(as.integer(d[i]))  \n\
    }  \n\
  }  \n\
  d  \n\
}  \n\
  \n\
percent <- function(d) {  \n\
  kmer = d[length(d)]  \n\
  for (i in 1:length(d)) {  \n\
    if (is.na(d[i])) {  \n\
      next  \n\
    }  \n\
    d[i] <- as.character(as.integer((as.numeric(d[i])*100)/kmer))  \n\
  }  \n\
  d  \n\
}  \n\
  \n\
format_theme <- function() {  \n\
    theme(legend.text = element_text(size=11),  \n\
          legend.position = \"right\",  \n\
          legend.title = element_text(angle = 90,hjust=.5),  \n\
          legend.key = element_blank(),  \n\
          plot.title = element_text(face=\"bold\",size=14,hjust=.5,vjust=.5,margin=margin(t=12,b=18,unit=\"pt\")),  \n\
          axis.title=element_text(size=12),  \n\
          axis.text=element_text(size=11),  \n\
          plot.background = element_blank(),  \n\
          panel.background = element_blank(),  \n\
          panel.grid.major = element_blank(),  \n\
          panel.grid.minor = element_blank(),  \n\
	  panel.border = element_rect(colour=\"black\", fill=NA, size=2))  \n\
}  \n\
  \n\
layer_on_contour <- function(with) {  \n\
  if (with)  \n\
    geom_contour(aes(z=Count), colour=\"white\", bins=6, show.legend=FALSE)  \n\
}  \n\
  \n\
plot_heat <- function(dat, source, h, with) {  \n\
  \n\
  y_max <- max(dat[,1]) + .5;  \n\
  x_max <- max(dat[,2]) + .5;  \n\
  \n\
  ggplot(data=dat, aes(x=KF,y=GCP,z=Count)) +  \n\
    geom_tile(aes(fill=Count)) +  \n\
    layer_on_contour(with) +  \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max), expand=FALSE) +  \n\
    scale_fill_viridis(name=\"Distinct K-mers per bin\", labels=fancy_scientific) +  \n\
    scale_y_continuous(breaks=seq(0,y_max,y_max/4), labels=percent) +  \n\
    guides(fill = guide_colorbar(title.position = \"left\", ticks = FALSE, draw.llim = TRUE, barheight = unit(h-1.15,\"in\"))) +  \n\
    xlab(paste(y_max,\"-mer frequency\", sep=\"\")) +  \n\
    ylab(\"GC percent\") +  \n\
    ggtitle(paste(\"KatGC of\", source, sep=\" \")) +  \n\
    format_theme()  \n\
}  \n\
  \n\
plot_contour <- function(dat, source) {  \n\
  \n\
  y_max <- max(dat[,1]) + .5;  \n\
  x_max <- max(dat[,2]) + .5;  \n\
  \n\
  ggplot(data=dat, aes(x=KF,y=GCP,z=Count)) +  \n\
    geom_contour(aes(z=Count, colour=..level..), bins=10, show.legend=TRUE) +  \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max), expand=FALSE) +  \n\
    scale_color_viridis(option=\"turbo\", name=\"Distinct K-mers per bin\", labels=fancy_scientific, breaks = 6) +  \n\
    scale_y_continuous(breaks=seq(0,y_max,y_max/4), labels=percent) +  \n\
    guides(color = guide_legend(title.position = \"left\", ticks = FALSE, draw.llim = TRUE, reverse=TRUE, keyheight=unit(24,\"pt\"))) +  \n\
    xlab(paste(y_max,\"-mer frequency\", sep=\"\")) +  \n\
    ylab(\"GC percent\") +  \n\
    ggtitle(paste(\"KatGC of\", source, sep=\" \")) +  \n\
    format_theme()  \n\
}  \n\
  \n\
save_plot <- function(name, type, outformat, h, w) {  \n\
  ggsave(file = paste(name, type, outformat, sep = \".\"), height = h, width = w)  \n\
}  \n\
  \n\
my_plot  <-  function(hist, name, h=4.5, w=6.0, type, pdf=FALSE, source) {  \n\
  \n\
  dat=read.table(hist, header=TRUE)  \n\
  \n\
  outformat=\"png\"  \n\
  if (pdf) {  \n\
    outformat=\"pdf\"  \n\
  }  \n\
  \n\
  if (type == \"combo\") {  \n\
    plot_heat(dat, source, h, TRUE)  \n\
    save_plot(name, \"st\", outformat, h, w)  \n\
  }  \n\
  else if (type == \"heat\") {  \n\
    plot_heat(dat, source, h, FALSE)  \n\
    save_plot(name, \"fi\", outformat, h, w)  \n\
  }  \n\
  else { # type == \"contour\" \n\
    plot_contour(dat, source)  \n\
    save_plot(name, \"ln\", outformat, h, w)  \n\
  }  \n\
  \n\
}  \n\
  \n\
my_plot(hist = args$file, name = args$output, h = args$ydim, w = args$xdim, type = args$type, pdf = args$pdf, source = args$source)  \n\
";
