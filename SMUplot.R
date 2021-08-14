#!/usr/bin/env Rscript

require("argparse")
require("ggplot2")
require("scales")
require("viridis")
require("cowplot")

parser <- ArgumentParser(description = "Make spectra-cn plots. Line, filled, and stacked spectra-cn plots will be generated.")
parser$add_argument("-f", "--file", type="character", help=".spectra-cn.hist file (required)", default=NULL)
parser$add_argument("-o", "--output", type="character", help="output prefix (required)")
parser$add_argument("-x", "--xdim", type="double", default=5, help="width of plot [default %(default)s]")
parser$add_argument("-y", "--ydim", type="double", default=5, help="height of plot [default %(default)s]")
parser$add_argument("-t", "--type", type="character", default="all", help="available types: line, fill, stack, or all. [default %(default)s]")
parser$add_argument("-p", "--pdf", dest='pdf', default=FALSE, action='store_true', help="get output in .pdf. [default .png]")
parser$add_argument("-s1", "--source1", type="character", help="source .ktab file (required)", default=NULL)
args <- parser$parse_args()

fancy_scientific <- function(d) {
  if (d[2] > 1000000) {
    for (i in 1:length(d)) {
      if (is.na(d[i])) {
        next
      }
      d[i] <- paste( as.character(as.integer(d[i])/1000000), "M", sep="")
    }
  } else if (d[2] > 1000) {
    for (i in 1:length(d)) {
      if (is.na(d[i])) {
        next
      }
      d[i] <- paste( as.character(as.integer(d[i])/1000), "K", sep="")
    }
  } else {
    for (i in 1:length(d)) {
      if (is.na(d[i])) {
        next
      }
      d[i] <- as.character(as.integer(d[i]))
    }
  }
  d
}

format_theme <- function() {
    theme(legend.text = element_text(size=8),
          legend.title = element_text(hjust=.5, size=8),
          legend.margin = margin(t=3, b=3, l=3, r=3, unit='pt'),
          plot.title = element_text(face="bold",size=14,hjust=.5,vjust=.5,
                                    margin=margin(t=12,b=18,unit="pt")),
          plot.margin =unit(c(0,0,6,6),"pt"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=11),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
	  panel.border = element_rect(colour="black", fill=NA, size=2))
}

layer_on_contour <- function(with) {
  if (with)
    geom_contour(aes(z=Count), colour="white", bins=6, show.legend=FALSE)
}

plot_heat <- function(dat, h, with) {

  y_max <- max(dat[,1]) + .5;
  x_max <- max(dat[,2]) + .5;

  M <- ggplot(data=dat, aes(x=KF2,y=KF1,z=Count)) +
    geom_tile(aes(fill=Count)) +
    layer_on_contour(with) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max), expand=FALSE) +
    scale_fill_viridis(name="kmer pairs", labels=fancy_scientific) +
    scale_x_continuous(breaks=c(12,17,25,33,50), labels=c("1/8","1/6","1/4","1/3","1/2")) +
    guides(fill = guide_colourbar(title.position = "top", frame.colour="black",
                                  ticks.colour="black", barheight = unit(.14*h,"in"))) +
    xlab(paste("Normalized minor kmer coverage: B/(A+B)", sep="")) +
    ylab(paste("Total coverage of the kmer pair: A+B", sep="")) +
    format_theme();

  M
}

plot_contour <- function(dat, h) {

  y_max <- max(dat[,1]) + .5;
  x_max <- max(dat[,2]) + .5;
  
  M <- ggplot(data=dat, aes(x=KF2,y=KF1,z=Count)) +
    geom_contour(aes(z=Count, colour=..level..), bins=12, show.legend=TRUE) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max), expand=FALSE) +
    scale_color_viridis(option="turbo", name="kmer pairs",
                        labels=fancy_scientific, n.breaks = 6) +
    guides(color = guide_legend(title.position = "top", frame.colour="black",
                                ticks.colour="black", reverse=TRUE,
                                keyheight = unit(2.4*h,"pt"))) +
    xlab(paste("Normalized minor kmer coverage: B/(A+B)", sep="")) +
    ylab(paste("Total coverage of the kmer pair: A+B", sep="")) +
    format_theme()

  M
}

save_plot <- function(name, type, outformat, h, w) {
  ggsave(file = paste(name, type, outformat, sep = "."), height = h, width = w)
}

my_plot  <-  function(hist, name, h=5, w=5, x_max, type, pdf=FALSE, s1) {

  dat=read.table(hist, header=TRUE)

  outformat="png"
  if (pdf) {
    outformat="pdf"
  }

  if (type == "combo")
    M <- plot_heat(dat, h, TRUE)
  else if (type == "heat")
    M <- plot_heat(dat, h, FALSE)
  else # type == "contour"
    M <- plot_contour(dat, h)

  hc = rgb(0.8352, 0.2431, 0.3098);

  E <- get_legend(M);

  M <- M + theme(legend.position='none');

  cvr <- aggregate(x=dat$Count,by=list(dat$KF1),FUN=max);
  C <- ggplot(cvr, aes(x=Group.1,y=x)) + geom_col(fill=hc, color="black", width=1.0) +
         coord_flip(xlim=c(0,50), expand=FALSE) +
         theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
         theme(plot.margin=unit(c(0,3,0,3),"pt")) +
         theme(panel.grid=element_blank(), panel.background=element_blank());

  hpr <- aggregate(x=dat$Count,by=list(dat$KF2),FUN=max)
  H <- ggplot(hpr, aes(x=Group.1,y=x)) + geom_col(fill=hc, color="black", width=1.0) +
         coord_cartesian(xlim=c(0,50), expand=FALSE) +
         theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
         theme(plot.margin=unit(c(3,0,3,0),"pt")) +
         theme(panel.grid=element_blank(), panel.background=element_blank());

  lmargin = (5*(.0865+.0166))/w;
  bmargin = (5*(.0750+.0166))/h;

  ggdraw() +
    draw_plot(M,x=0.0,y=0.0,width=0.8,height=0.8) +
    draw_plot(H,x=lmargin,y=0.8,width=0.80-lmargin,height=0.20) +
    draw_plot(C,x=0.80,y=bmargin,width=0.20,height=0.80-bmargin) +
    draw_plot(E,x=0.80,y=0.80,width=0.20,height=0.20) +
    draw_label(paste(s1), x=0.25/w, y=.95, hjust=0, fontface="bold.italic", size=40/.pt) +
    draw_text("estimated diploid", x=0.35/w, y=.90, hjust=0, size=32/.pt) +
    draw_text( c("AB     .60", "AAB    .10", "AABB   .05", "AAABB  .02", "AAABBB .01"),
               x = c(.83,.83,.83,.83,.83),
               y = .76 - c(.0,.03,.06,.09,.12) * (5/h),
               hjust=0, vjust=0, size=28/.pt, family="mono") +
    draw_text("1n = 158", x=.85, y=.31/h, hjust=0, vjust=0, size=32/.pt);

  if (type == "combo")
    save_plot(name=name, type="st", outformat, h=h, w=w)
  else if (type == "heat")
    save_plot(name=name, type="fi", outformat, h=h, w=w)
  else  # type == "contour"
    save_plot(name=name, type="ln", outformat, h=h, w=w)
}

my_plot(hist = args$file, name = args$output, h = args$ydim, w = args$xdim, type = args$type, pdf = args$pdf, s1 = args$source1)
