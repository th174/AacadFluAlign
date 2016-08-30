

----------------USAGE HELP----------------

EXAMPLE USAGE:		perl AacadFluAlign.pl -in PB1remainder-all.fasta -seg PB1 -clean


Command Line Flags:

        -in <file.fasta>				Name of input file															(Mandatory)
        -seg <gene segment>				Gene segment to input														(Mandatory)
        -out <file.tsv>					Name of output file															(Default: <same-as-input.tsv>)
        -maxSeqs <integer>				Maximum sequences per MUSCLE alignment										(Default: 1)
        -maxThreads <integer>			Maximum number of parallel threads used in analysis							(Default: 32)
        -clean							Cleanup temporary alignment files after completing 							(Default: FALSE)
        -verbose						Display detailed runtime processes											(Default: FALSE)
        -ref <reference.fasta>			Name of reference sequence													(Default: Automatic)
        -showAll						Output all sequences not only flagged sequences								(Default: FALSE)
        -IRD							Use IRD error codes															(Default: FALSE)
        -hoffMatchMin <integer>			Minimum number of base pair matches to be classified as Hoffman Primer		(Default: 5)
        -blastMatchMin <integer>		Minimum size of extraneous sequence to be classified as cloning plasmid		(Default: 12)
        -printEndLength <integer>		Number of base pairs to print on 3' and 5' ends								(Default: 50)


To install:

First, make sure to have an updated version of Perl 5 installed. Visit www.perl.org/get.html to obtain the correct perl distribution for your system.

Then, download and unzip all the files in this repository into a local directory.

Navigate to the directory through the command line, and run this program by calling "perl AacadFluAlign.pl"

See the Usage Help above for more details.

---------------------------------------------------------------------------------------

This program was tested on three platforms, Windows 10, Linux, and Mac OS X.
