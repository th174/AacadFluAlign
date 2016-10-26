

----------------USAGE HELP----------------

EXAMPLE USAGE:		perl AacadFluAlign.pl -in PB1remainder-all.fasta -seg PB1 -clean


Command Line Flags:

        -in <file.fasta>				Name of input file															(Mandatory)
        -seg <gene segment>				Gene segment to input														(Mandatory, only internal genes supported)
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

This program operates by doing a systematic analysis of MUSCLE MSA results. It compares an aligned sample sequence dataset with a reference sequence, and attempts to identify and classify possible sequencing errors using procedures established in the Athens Academy Influenza Genome Sequencing Project. AacadFluAlign then generates a TSV (tab-separated-values) spreadsheet that presents and formats errors for humans to curate. 

AacadFluAlign first separates the sample FASTA dataset into smaller, more manageable files (Default size: 1 sequence per file). This facilitates multithreaded data processing through the usage of perlfork, optimizing CPU usage and significantly cutting down on computation time. Each separate alignment is firsed parsed then analyzed using an adaptation of Athens Academy Influenza Project's methodology. 

First, the program identifies and classifies extraneous sequences. This program recognizes any sample sequence aligned before the first base in the reference sequence as an extraneous sequence on the 5' end, and any sample sequence aligned after the last base in the reference sequences is recognized as an extraneous sequence on the 3' end. The program then attempts to classify the extraneous sequence based on its contents. 

        -- An extraneous lone 5' T or a lone 3' A is recognized as a Taq Polymerase error and indicated by Taq on the spreadsheet
        -- An extraneous sequence consistent with published influenza primers by Hoffman et al. is classified as a Primer error and indicated by Pri on the spreadsheet
        -- An extraneous sequence with length more than 12 bp is suspected to be cloning plasmid sequence, indicated with "clo" which the spreadsheet prompts the user to verify with BLAST
        -- All other extraneous sequences are classified as unknown, indicated as "unk".
        
Next, the program identifies potential errors in the conserved terminal regions, or the first 12 bp and the last 13 bp, of the chosen gene segment. Any base pairs that disagree between the reference sequence and the sample sequence within the first 12 bp or last 13 bp of the reference sequence are classified as a CTS error, and indicated as such on the spreadsheet.

Then, program then proceeds to identify and classify potential errors in the internal portion of the sample gene segment using the following algorithm. 

        -- A dash(-) in the reference sequence coinciding with a base pair in the sample sequence is classified as an insertion and indicated by "ins" on the spreadsheet.
        -- A base pair in the reference sequence coinciding with a dash(-) in the sample sequence is classified as a deletion and indicated as "del" on the spreadsheet.
        -- Adenine insertions or deletions in the sample sequence in positions coinciding with the Poly A tail of the reference sequence are classified as "Poly A", and are noted as such in the spreadsheet.
        -- Consecutive insertions and deletions of length equal to are multiple of 3 are classified as triplets, and are noted "3x" in the spreadsheet.
The program then further classifies these errors by the location in which they occur within the reference sequence. Errors ocurring before the start codon and after the stop codon in the reference sequence are classified as NCR (Non-coding region) erorrs, while errors ocurring between the start and stop codon are classified as CDS (Coding sequence) errors.

Next, the program records the number of mixed bases that occur in the sample sequence. A mixed base is defined as the bases R,Y,S,W,K,M,B,D,H,V, and N in the sample sequence. If the total amount of mixed bases is found to exceed 0.5% of length of the sample sequence, then they are indicated by "mix" as well as the region in which they occur, on the spreadsheet.

Over the course of this project, the above algorithm frequently and erroneously classified MUSCLE alignment mistakes as sequence errors. These alignment mistakes were generally characterized in the sample sequence by a short small number of base pairs, or orphan sequence, separated from the rest of the gene segment by large numbers of consecutive dashes. This program attempts to reduce the rate of false positives by identifying orphan sequences and removing them from the spreadsheet. This program classifies orphan sequences with the following rules

        -- An orphan sequence must either begin and end within the first 2% of the gene segment or the last 2% of the gene segment.
        -- An orphan sequence must be separated from the rest of the gene segment by a consecutive sequence of dashes with length 20 or more.
Sequences that do not meet the two rules above but are close to meeting them are identified as possible orphans on the spreadsheet, and prompts the user to manually check if deemed necessary.

Finally, this program outputs its data using a TSV (tab-separated values) formatted spreadsheet. TSV spreadsheets may be opened using Microsoft Excel or any similar spreadsheet editor. Do note that when importing a TSV into Excel, it is necessary to choose "none" for the text delimiter field, and explicitly specify that each column is text. Otherwise, some versions of Excel will assume that a dash(-) is a subtraction sign in a mathematical equation, and throw an error as a result.

-------------------------------------------------------------------------------------------------------------------

Credit for MUSCLE v3.0 goes to Bob Edgar http://www.drive5.com/muscle/manual/papers.html

Credit for this program goes to Timmy Huang, Nikki Chester, David Suarez, and the Athens Academy Influenza Genome Project
