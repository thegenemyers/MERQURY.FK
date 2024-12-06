DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = HAPmaker CNplot ASMplot HAPplot MerquryFK KatComp KatGC

all: $(ALL)

libfastk.c : gene_core.c libfastk.h
libfastk.h : gene_core.h

cn_plotter.c : cn_plot.R.h cn_plotter.h
asm_plotter.c : cn_plot.R.h asm_plotter.h
hap_plotter.c : hap_plot.R.h hap_plotter.h

HAPmaker: HAPmaker.c libfastk.c
	gcc $(CFLAGS) -o HAPmaker HAPmaker.c libfastk.c -lpthread -lm

CNplot: CNplot.c cn_plotter.c libfastk.c
	gcc $(CFLAGS) -o CNplot CNplot.c cn_plotter.c libfastk.c -lpthread -lm

ASMplot: ASMplot.c asm_plotter.c libfastk.c
	gcc $(CFLAGS) -o ASMplot ASMplot.c asm_plotter.c libfastk.c -lpthread -lm

HAPplot: HAPplot.c hap_plotter.c libfastk.c
	gcc $(CFLAGS) -o HAPplot HAPplot.c hap_plotter.c libfastk.c -lpthread -lm

MerquryFK: MerquryFK.c cn_plotter.c asm_plotter.c hap_plotter.c blk_plot.R.h libfastk.c
	gcc $(CFLAGS) -o MerquryFK MerquryFK.c cn_plotter.c asm_plotter.c hap_plotter.c libfastk.c -lpthread -lm

KatComp: KatComp.c libfastk.c kx_plot.R.h
	gcc $(CFLAGS) -o KatComp KatComp.c libfastk.c -lpthread -lm

KatGC: KatGC.c libfastk.c kgc_plot.R.h
	gcc $(CFLAGS) -o KatGC KatGC.c libfastk.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f MerquryFK.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf MerquryFK.tar.gz Makefile *.h *.c README.md
