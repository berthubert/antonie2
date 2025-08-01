-include sysdeps/$(shell uname).inc

VERSION=0.1
CXXFLAGS?=-Wall -O3 -fPIC -I/usr/include/python3.6m -ggdb -I. -Iext -Iext/libmba  -MMD -MP -pthread  $(CXX2014FLAGS) -Wno-strict-aliasing -Iext/nr_c304/code -Iext/fmt-8.0.1/include -fconcepts -Wno-reorder # -Wno-unused-local-typedefs 
CFLAGS=-Wall -I. -Iext/libmba -O3 -MMD -MP
LDFLAGS=$(CXX2014FLAGS) -pthread  # -Wl,-Bstatic -lstdc++ -lgcc -lz -Wl,-Bdynamic -static-libgcc -lm -lc
CHEAT_ARG := $(shell ./update-git-hash-if-necessary)

SHIPPROGRAMS=antonie 16ssearcher stitcher  fqgrep pfqgrep genex
PROGRAMS=$(SHIPPROGRAMS) digisplice gffedit gfflookup gfftool nwunsch fogsaa gtfreader cor2 correlo gcstats resscan hpcount skfit gcscan chagstats cpgstats genehisto\
	charsample chromopic exoexplore cluster simsearch fagrep

ifeq ($(CC),clang)
        CXXFLAGS+=-ftemplate-depth=1000
endif

all: $(PROGRAMS)

-include *.d

.PHONY:	antonie.exe codedocs/html/index.html check

MBA_OBJECTS = ext/libmba/allocator.o ext/libmba/diff.o ext/libmba/msgno.o ext/libmba/suba.o ext/libmba/varray.o 
ANTONIE_OBJECTS = antonie.o refgenome.o hash.o geneannotated.o misc.o fastq.o saminfra.o dnamisc.o githash.o phi-x174.o zstuff.o genbankparser.o $(MBA_OBJECTS)

dino: dino.o 
	$(CXX) $^ -o $@

strdiff: strdiff.o $(MBA_OBJECTS)
	$(CC) strdiff.o $(MBA_OBJECTS) -o $@

antonie: $(ANTONIE_OBJECTS)
	$(CXX) $(ANTONIE_OBJECTS) $(LDFLAGS) $(STATICFLAGS) -lz -o $@

SEARCHER_OBJECTS=16ssearcher.o hash.o misc.o fastq.o zstuff.o githash.o fastqindex.o stitchalg.o

16ssearcher: $(SEARCHER_OBJECTS)
	$(CXX)  $(SEARCHER_OBJECTS) -lz  $(LDFLAGS) $(STATICFLAGS) -o $@

digisplice: digisplice.o refgenome.o misc.o fastq.o hash.o zstuff.o dnamisc.o geneannotated.o genbankparser.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

stitcher: stitcher.o refgenome.o misc.o fastq.o hash.o zstuff.o dnamisc.o geneannotated.o genbankparser.o fastqindex.o stitchalg.o
	$(CXX) $(LDFLAGS) $^ -lz -pthread $(STATICFLAGS) -o $@

#renovo: renovo.o refgenome.o misc.o fastq.o hash.o zstuff.o dnamisc.o geneannotated.o genbankparser.o fastqindex.o stitchalg.o
#	$(CXX) $(LDFLAGS) $^ -lz -pthread $(STATICFLAGS) -o $@


libbridge.so: bridge.o
	g++ -shared -Wl,-soname,"libhello.so" bridge.o -lboost_python3 -fpic -o libbridge.so


invert: invert.o misc.o
	$(CXX) $(LDFLAGS) $(STATICFLAGS) $^ -lz -o $@

fqgrep: fqgrep.o misc.o fastq.o dnamisc.o zstuff.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

pfqgrep: pfqgrep.o misc.o fastq.o dnamisc.o zstuff.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

genex: genex.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -o $@

correlo: correlo.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

gcstats: gcstats.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

resscan: resscan.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@


hpcount: hpcount.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

chagstats: chagstats.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

cpgstats: cpgstats.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

simsearch: simsearch.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

fagrep: fagrep.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@


genehisto: genehisto.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@


chromopic: chromopic.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@


charsample: charsample.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o taxoreader.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@


gcscan: gcscan.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@



skfit: skfit.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o geneannotated.o ext/fmt-8.0.1/src/format.o ext/fmt-8.0.1/src/os.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@

cor2: cor2.o dnamisc.o zstuff.o misc.o hash.o nucstore.o refgenome2.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -pthread -lbz2 -o $@


gffedit: gffedit.o refgenome.o fastq.o dnamisc.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

exoexplore: exoexplore.o refgenome.o fastq.o dnamisc.o refgenome2.o genbankparser.o geneannotated.o nucstore.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

cluster: cluster.o refgenome.o fastq.o dnamisc.o refgenome2.o genbankparser.o geneannotated.o nucstore.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@


gfflookup: gfflookup.o geneannotated.o genbankparser.o refgenome2.o nucstore.o fastq.o dnamisc.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

gfftool: gfftool.o geneannotated.o genbankparser.o refgenome2.o nucstore.o fastq.o dnamisc.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@


gtfreader: gtfreader.o geneannotated.o genbankparser.o refgenome2.o nucstore.o fastq.o dnamisc.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@


gendump: gendump.o geneannotated.o genbankparser.o refgenome2.o nucstore.o fastq.o dnamisc.o zstuff.o misc.o hash.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@


nwunsch: nwunsch.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@

fogsaa: fogsaaimp.o
	$(CXX) $(LDFLAGS) $^ -lz $(STATICFLAGS) -o $@


install: antonie
	mkdir -p $(DESTDIR)/usr/bin/
	mkdir -p $(DESTDIR)/usr/share/doc/antonie/
	mkdir -p $(DESTDIR)/usr/share/doc/antonie/ext
	install -s $(SHIPPROGRAMS) $(DESTDIR)/usr/bin/
	cp report.html $(DESTDIR)/usr/share/doc/antonie
	cp -r ext/html $(DESTDIR)/usr/share/doc/antonie/ext

clean:
	rm -f *~ *.o $(MBA_OBJECTS) *.d $(PROGRAMS) githash.h 

package: all
	rm -rf dist
	DESTDIR=dist make install
	fpm -s dir -f -t rpm -n antonie -v 1.g$(shell cat githash) -C dist .
	fpm -s dir -f -t deb -n antonie -v 1.g$(shell cat githash) -C dist .	
	rm -rf dist

codedocs: codedocs/html/index.html

codedocs/html/index.html: 	
	doxygen

antonie.exe: 
	make clean
	STATICFLAGS="-static -static-libgcc -static-libstdc++" CXX=i686-w64-mingw32-g++ CC=i686-w64-mingw32-gcc make antonie
	mv antonie antonie.exe

16ssearcher.exe: 
	make clean
	CXXFLAGS="-Wall -O3 -I. -Iext/libmba -MMD -MP  $(CXX2014FLAGS)" STATICFLAGS="-static -static-libgcc -static-libstdc++" CXX=i686-w64-mingw32-g++ CC=i686-w64-mingw32-gcc make 16ssearcher
	mv 16ssearcher 16ssearcher.exe

check: testrunner
	./testrunner

testrunner: test-misc_hh.o test-nucstore_cc.o test-dnamisc_cc.o test-saminfra_cc.o testrunner.o misc.o dnamisc.o saminfra.o zstuff.o fastq.o hash.o nucstore.o taxoreader.o test-taxoreader_cc.o
	$(CXX) $^ -lboost_unit_test_framework -lz -o $@ 
