package = ssenss
version = 1.0.1
distdir = $(package)_$(version)

CC = mpicc
CFLAGS = -std=c99
LIBS = -lm -llapack -lblas -lhdf5
release: CFLAGS += -O3
debug: CFLAGS += -g -Wall -Wpedantic

export CC
export CFLAGS
export LIBS
export package

default: $(package)
all: $(package)
$(package): release

clean:
	cd src && $(MAKE) $@

debug:
	cd src && $(MAKE) -B $(package)
	cp src/$(package) .
	cd src && $(MAKE) clean

release:
	cd src && $(MAKE) $(package)
	cp src/$(package) .

dist: $(distdir).tar.gz

$(distdir).tar.gz: $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $@
	rm -rf $(distdir)

$(distdir):
	mkdir -p $(distdir)/src
	mkdir $(distdir)/doc
	mkdir $(distdir)/examples
	mkdir $(distdir)/misc
	cp Makefile $(distdir)
	cp README $(distdir)
	cp COPYING $(distdir)
	cp COPYRIGHT $(distdir)
	cp src/Makefile $(distdir)/src
	cp src/*.c $(distdir)/src
	cp src/*.h $(distdir)/src
	cp examples/setup $(distdir)/examples
	cp examples/*.h5 $(distdir)/examples

.PHONY: default all clean dist debug release
