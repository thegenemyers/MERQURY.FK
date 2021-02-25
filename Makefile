DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = CNplot ASMplot

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

CNplot: CNplot.c cn.R.h libfastk.c libfastk.h
	gcc $(CFLAGS) -o CNplot CNplot.c libfastk.c -lpthread -lm

ASMplot: ASMplot.c cn.R.h libfastk.c libfastk.h
	gcc $(CFLAGS) -o ASMplot ASMplot.c libfastk.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f MerquryFK.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf MerquryFK.tar.gz Makefile *.h *.c README.md
