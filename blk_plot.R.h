static char *blk_plot_script =
"#!/usr/bin/env Rscript   \n\
   \n\
require(\"argparse\")   \n\
require(\"ggplot2\")   \n\
require(\"scales\")   \n\
   \n\
parser <- ArgumentParser(description = \"Make block N* or NG* plots. Applicable for scaffolds, contigs, and phased blocks.\")   \n\
parser$add_argument(\"-b\", \"--block\", type=\"character\", help=\"sorted .sizes file\", default=NULL)   \n\
parser$add_argument(\"-s\", \"--scaff\", type=\"character\", help=\"sorted .sizes file\", default=NULL)   \n\
parser$add_argument(\"-c\", \"--contig\", type=\"character\", help=\"sorted .sizes file\", default=NULL)   \n\
parser$add_argument(\"-o\", \"--out\", type=\"character\", help=\"output prefix (required) [default %(default)s]\", default=\"out\")   \n\
parser$add_argument(\"-g\", \"--gsize\", type=\"integer\", default=0, help=\"genome size for computing NG* (optional)\")   \n\
parser$add_argument(\"-x\", \"--xdim\", type=\"double\", default=6, help=\"width of plot [default %(default)s]\")   \n\
parser$add_argument(\"-y\", \"--ydim\", type=\"double\", default=5, help=\"height of plot [default %(default)s]\")   \n\
parser$add_argument(\"-p\", \"--pdf\", dest=\'pdf\', default=FALSE, action=\'store_true\', help=\"set to get output in .pdf. [default .png]\")   \n\
args <- parser$parse_args()   \n\
   \n\
fancy_scientific <- function(d) {   \n\
  # turn in to character string in scientific notation   \n\
  d <- format(d, scientific = TRUE)   \n\
  # quote the part before the exponent to keep all the digits and turn the \'e+\' into 10^ format   \n\
  d <- gsub(\"^(.*)e\\\\+\", \"'\\\\1'%*%10^\", d)  \n\
  # convert 0x10^00 to 0  \n\
  d <- gsub(\"\\\\'0[\\\\.0]*\\\\'(.*)\", \"'0'\", d)  \n\
  # return this as an expression   \n\
  parse(text=d)   \n\
}   \n\
   \n\
save_plot <- function(name, type, stats, outformat, h, w) {   \n\
  if (outformat == \"pdf\") {   \n\
    ggsave(file = paste(name, type, stats, outformat, sep = \".\"), height = h, width = w, units = \"in\")   \n\
  } else {   \n\
    ggsave(file = paste(name, type, stats, outformat, sep = \".\"), height = h, width = w, dpi=300)   \n\
  }   \n\
   \n\
}   \n\
   \n\
attach_n <- function(dat, gsize=0) {   \n\
  dat = read.table(dat, header = F)   \n\
  names(dat) = c(\"Type\", \"Group\", \"Size\")   \n\
  dat$Sum = cumsum(as.numeric(dat$Size))   \n\
     \n\
  if (gsize == 0) {   \n\
    # N*   \n\
    gsize = sum(dat$Size)   \n\
  }   \n\
  dat$N = 100*dat$Sum/gsize   \n\
  dat$N2 = 100*c(0, dat$Sum[-nrow(dat)]/gsize)   \n\
  return(dat)   \n\
}   \n\
   \n\
get_dummy <- function(dat=NULL, type) {   \n\
  x_max=max(dat$N)   \n\
  data.frame(Type = c(type), Group = c(\"dummy\"), Size = c(0), Sum = c(0), N = c(x_max), N2 = c(x_max))   \n\
}   \n\
   \n\
bind_blocks <- function(block, block_dummy, scaff, scaff_dummy, contig, contig_dummy) {   \n\
     \n\
  blocks = data.frame()   \n\
  if (!is.null(block)) {   \n\
    blocks=rbind(blocks, block, block_dummy)   \n\
  }   \n\
  if (!is.null(scaff)) {   \n\
    blocks=rbind(blocks, scaff, scaff_dummy)   \n\
  }   \n\
  if (!is.null(contig)) {   \n\
    blocks=rbind(blocks, contig, contig_dummy)   \n\
  }   \n\
  return(blocks)   \n\
}   \n\
   \n\
plot_block <- function(dat = NULL, stats) {   \n\
  # by phased block   \n\
  y_max=max(dat$Size)   \n\
  ggplot(data = dat, aes(x = dat[,5], y = dat[,3], fill = dat[,2], colour = dat[,2])) +   \n\
    geom_rect(xmin=dat[,6], xmax=dat[,5], ymin=0, ymax=dat[,3], alpha=0.7) +   \n\
    theme_bw() +   \n\
    theme(legend.text = element_text(size=11),   \n\
          legend.position = c(0.95,0.95),  # Modify this if the legend is covering your favorite circle   \n\
          legend.background = element_rect(size=0.1, linetype=\"solid\", colour =\"black\"),   \n\
          legend.box.just = \"right\",   \n\
          legend.justification = c(\"right\", \"top\"),   \n\
          axis.title=element_text(size=14,face=\"bold\"),   \n\
          axis.text=element_text(size=12)) +   \n\
    scale_fill_brewer(palette = \"Set1\", name= names(dat)[2]) +   \n\
    scale_colour_brewer(palette = \"Set1\", name= names(dat)[2]) +   \n\
    scale_x_continuous(limits = c(0, 100)) +   \n\
    scale_y_continuous(limits = c(0, y_max), labels = fancy_scientific) +   \n\
    xlab(stats) + ylab(\"Size (bp)\") +   \n\
    geom_vline(xintercept = 50, show.legend = FALSE, linetype=\"dashed\", color=\"black\")   \n\
}   \n\
   \n\
block_n <- function(block=NULL, scaff=NULL, contig=NULL, out, gsize = 0, w = 6, h = 5, pdf=FALSE) {   \n\
     \n\
  outformat=\"png\"   \n\
  if (pdf) {   \n\
    outformat=\"pdf\"   \n\
  }   \n\
     \n\
  # Read file   \n\
  if (!is.null(block)) {   \n\
    block = attach_n(dat = block, gsize = gsize)   \n\
    block_dummy = get_dummy(dat = block, type = \"block\")   \n\
  }   \n\
     \n\
  if (!is.null(scaff)) {   \n\
    scaff = attach_n(dat = scaff, gsize = gsize)   \n\
    scaff_dummy = get_dummy(dat = scaff, type = \"scaffold\")   \n\
  }   \n\
     \n\
  if (!is.null(contig)) {   \n\
    contig = attach_n(dat = contig, gsize = gsize)   \n\
    contig_dummy = get_dummy(dat = contig, type = \"contig\")   \n\
  }   \n\
     \n\
  stats = \"NG\"   \n\
  if (gsize == 0) {   \n\
    stats = \"N\"   \n\
  }   \n\
   \n\
  # Plot phase blocks filled by haplotypes   \n\
  plot_block(block, stats)   \n\
  save_plot(out, \"block\", stats, outformat, h = h, w = w)   \n\
     \n\
  dat = bind_blocks(block, block_dummy, scaff, scaff_dummy, contig, contig_dummy)   \n\
  y_max=max(dat$Size)   \n\
     \n\
  ggplot(data = dat, aes(x = N2, y = Size, colour = Type)) +   \n\
    geom_step() +   \n\
    theme_bw() +   \n\
    theme(legend.text = element_text(size=11),   \n\
          legend.position = c(0.95,0.95),  # Modify this if the legend is covering your favorite circle   \n\
          legend.background = element_rect(size=0.1, linetype=\"solid\", colour =\"black\"),   \n\
          legend.box.just = \"right\",   \n\
          legend.justification = c(\"right\", \"top\"),   \n\
          axis.title=element_text(size=14,face=\"bold\"),   \n\
          axis.text=element_text(size=12)) +   \n\
    scale_fill_brewer(palette = \"Set1\", name = \"Type\") +   \n\
    scale_colour_brewer(palette = \"Set1\", name= \"Type\") +   \n\
    scale_x_continuous(limits = c(0, 100)) +   \n\
    scale_y_continuous(limits = c(0, y_max), labels = fancy_scientific) +   \n\
    xlab(stats) + ylab(\"Size (bp)\") +   \n\
    geom_vline(xintercept = 50, show.legend = FALSE, linetype=\"dashed\", color=\"black\")   \n\
     \n\
  save_plot(out, \"continuity\", stats, outformat, h = h, w = w)   \n\
}   \n\
   \n\
block_n(block = args$block, scaff = args$scaff, contig = args$contig, out = args$out, gsize = args$gsize, w = args$xdim, h = args$ydim, pdf = args$pdf)   \n\
";
