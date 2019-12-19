# ephem - astronomical calculations and ephemeris for the MMT
#
# Tom Trebisky  9-25-2017
# Tom Trebisky  11-27-2018

# This seems to be the biggest font that does not cause
# enscript to wrap lines.
FONT = Courier9

# Set this to last 2 digits of the desired year.
Y = 17

# Change this name for a new year
BASE=almn$(Y)
YEAR = 20$(Y)

all:	$(BASE).ps $(BASE).pdf

$(BASE).txt:	ephem
	./ephem -a$(YEAR)> $(BASE).txt

$(BASE).ps:	$(BASE).txt
	enscript --landscape --no-header --font $(FONT) --output=$(BASE).ps $(BASE).txt

$(BASE).pdf:	$(BASE).ps
	ps2pdf $(BASE).ps $(BASE).pdf

ephem:	ephem.c
	cc -o ephem ephem.c -lm

install: ephem
	cp ephem /home/tom/bin

clean:
	rm -f *.o ephem
	rm -f *.ps *.txt *.pdf
