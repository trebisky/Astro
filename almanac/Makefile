# Makefile for the MMT almanac program
# Tom Trebisky  10-27-2005, 7-22-2009, 9-6-2011
# Tom Trebisky  9-18-2015 9-1-2016
# Tom Trebisky  11-27-2018


# This seems to be the biggest font that does not cause
# enscript to wrap lines.
FONT = Courier9

all:	almn19.ps almn19.pdf

almn19.f:	almn.proto mksrc mkdata
	./mksrc almn.proto >almn19.f

almn19:	almn19.f
	f77 -o almn19 almn19.f

almanac:	almn19 header
	rm -f almanac
	./almn19

almn19.txt:	almanac
	./mktxt almanac > almn19.txt

almn19.ps:	almn19.txt
	enscript --landscape --no-header --font $(FONT) --output=almn19.ps almn19.txt

almn19.pdf:	almn19.ps
	ps2pdf almn19.ps almn19.pdf

# We now do this on the nfs mounted /mmt from hacksaw, so don't need this.
#backup:
#	rsync -auv -e ssh /mmt/almanac hacksaw:/mmt

clean:
	rm -f almn19 almn19.f almn19.txt almanac almn19.ps almn19.pdf
