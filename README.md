Created for Athens Academy's Influenza Sequencing Project.

----------------USAGE HELP----------------
    EXAMPLE USAGE:\tperl AacadFluAlign.pl -in PB1remainder-all.fasta -seg PB1 -clean
    Command Line Flags:
        -in <file.fasta>\t\tName of input file\t\t\t\t\t\t\t(Mandatory)
            -seg <gene segment>\t\tGene segment to input\t\t\t\t\t\t\t(Mandatory)
            -out <file.tsv>\t\tName of output file\t\t\t\t\t\t\t(Default: <same-as-input.tsv>)
            -maxSeqs <integer>\t\tMaximum sequences per MUSCLE alignment\t\t\t\t\t(Default: 1)
            -maxThreads <integer>\tMaximum number of parallel threads used in analysis\t\t\t(Default: 32)
            -clean\t\t\tCleanup temporary alignment files after completing \t\t\t(Default: FALSE)
            -verbose\t\t\tDisplay detailed runtime processes\t\t\t\t\t(Default: FALSE)
            -ref <reference.fasta>\tName of reference sequence\t\t\t\t\t\t(Default: Automatic)
            -showAll\t\t\tOutput all sequences not only flagged sequences\t\t\t(Default: FALSE)
            -IRD\t\t\tUse IRD error codes\t\t\t\t\t\t\t(Default: FALSE)
            -hoffMatchMin <integer>\tMinimum number of base pair matches to be classified as Hoffman Primer\t(Default: 5)
            -blastMatchMin <integer>\tMinimum size of extraneous sequence to be classified as cloning plasmid\t(Default: 12)
            -printEndLength <integer>\tNumber of base pairs to print on 3' and 5' ends\t\t\t\t(Default: 50)

			
	
