# ephem - astronomical calculations and ephemeris for the MMT
#
# Tom Trebisky  9-25-2017
# Tom Trebisky  11-27-2018

# This seems to be the biggest font that does not cause
# enscript to wrap lines.
FONT = Courier9

# Change this name for a new year
BASE=almn20

all:	$(BASE).ps $(BASE).pdf

$(BASE).txt:	ephem
	./ephem > $(BASE).txt

$(BASE).ps:	$(BASE).txt
	enscript --landscape --no-header --font $(FONT) --output=$(BASE).ps $(BASE).txt

$(BASE).pdf:	$(BASE).ps
	ps2pdf $(BASE).ps $(BASE).pdf

ephem:	ephem.c
	cc -o ephem ephem.c -lm

clean:
	rm -f *.o ephem
	rm -f *.ps *.txt *.pdf
