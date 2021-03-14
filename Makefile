DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = CNplot ASMplot CNspectra

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

CNplot: CNplot.c cn_plotter.c cn_plotter.h libfastk.c libfastk.h
	gcc $(CFLAGS) -o CNplot CNplot.c cn_plotter.c libfastk.c -lpthread -lm

ASMplot: ASMplot.c asm_plotter.c asm_plotter.h libfastk.c libfastk.h
	gcc $(CFLAGS) -o ASMplot ASMplot.c asm_plotter.c libfastk.c -lpthread -lm

CNspectra: CNspectra.c cn_plotter.c cn_plotter.h asm_plotter.c asm_plotter.h libfastk.c libfastk.h
	gcc $(CFLAGS) -o CNspectra CNspectra.c cn_plotter.c asm_plotter.c libfastk.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f MerquryFK.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf MerquryFK.tar.gz Makefile *.h *.c README.md
