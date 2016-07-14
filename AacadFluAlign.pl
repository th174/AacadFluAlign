#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
# use warnings;


our $fileInput;                                                 
our $fileReference;
our $fileOutput;
our $boolShowCorrectSeq = 0;
our $boolIRDCodesOnly = 0;
our $boolCleanup = 0;
our $segment;
our $maxThreads = 32;
our $verbose = 0;

our $printEndLength = 50;

our $refName;

our $maxSequencesPerAlignment = 1;

our $lengthCTS5 = 12;                                                        
our $lengthCTS3 = 13;

our $NCR5;
our $NCR3;
our $regexPolyA;
our $hoff5 = "TATTCGTCTCAGGG";
our $hoff3= "AATACGAGACGATAT";
our $hoffThreshold = 5;
our $BLASTthreshold = 12;

GetOptions (
    'in=s'                      => \$fileInput,                                
    'ref=s'                     => \$fileReference,               
    'out=s'                     => \$fileOutput,     
    'maxThreads=i'              => \$maxThreads,
    'showAll'                   => \$boolShowCorrectSeq,                  
    'IRD'                       => \$boolIRDCodesOnly,                    
    'seg=s'                     => \$segment, 
    'verbose'                   => \$verbose,
    'hoffMatchMin=i'            => \$hoffThreshold,                       
    'blastMatchMin=i'           => \$BLASTthreshold,                      
    'printEndLength=i'          => \$printEndLength,                      
    'maxSeqs=i'                 => \$maxSequencesPerAlignment,            
    'clean'                     => \$boolCleanup                          
);

if (!$fileInput || !$segment){
    printHelp();
}

if (!$fileOutput) {
    $fileOutput = substr($fileInput,0,-6).".tsv";
}

if ($segment eq "PB1"){
    $NCR3 = "TAGTGAATTTAGCTTTGTCCTTCATGAAAAAATGCCTGTTTCTACT";
    $NCR5 = "AGCGAAAGCAGGCAAACCATTTGAATG";
    $regexPolyA = "TGAAA+TG";
} elsif($segment eq "PB2"){
    $NCR5 = "AGCGAAAGCAGGTCAAATATATTCAATATG";
    $NCR3 = "TAGTGTCGAATTGTTTAAAAACGACCTGTTTCTACT";
    $regexPolyA = "TTAAA+CG";
} elsif ($segment eq "PA") {
    $NCR5 = "AGCGAAAGCAGGTACTGATCCAAAATG";    #PA
    $NCR3 = "TAGTTGTGGCAATGCTACTATTTGCTATCCATACTGTCCAAAAAAGTACCTGTTTCTACT"; #PA
    $regexPolyA = "CCAAA+GT";
} elsif ($segment eq "MA") {
    $NCR5 = "AGCAAAAGCAGGTAGATATTGAAAGATG";   #MA
    $NCR3 = "TAGAGCTGGAGTAAAAAACTACCTTGTTTCTACT"; #MA
    $regexPolyA = "GTAAA+CT";;
} elsif ($segment eq "NS") {
    $NCR5 = "AGCAAAAGCAGGGTGACAAAAACATAATG";
    $NCR3 = "TGATAAAAAACACCCTTGTTTCTACT";
    $regexPolyA = "ATAAA+CA";
} elsif ($segment eq "NP") {
    $NCR5 = "AGCAAAAGCAGGGTAGATAATCACTCACCGTGTGACATCCACATCATG";
    $NCR3 = "TAAAGAAAAATACCCTTGTTTCTACT";
    $regexPolyA = "AGAAA+TA";
}
$NCR3 =~ /(.*)$regexPolyA/;
our $posPolyA = (length($1) > 1) ? length($1) : 1000;

our %fileRef = (
    PB1 => 'KR732492(PB1).fasta',
                PB2 => 'KR732460(PB2).fasta',
                PA  => 'EU742958(PA).fasta',
                NP  => 'DQ870889(NP).fasta',
                NS  => 'KT314336(NS).fasta',
                MA  => 'EU026075(MA).fasta'
);

#--Main Begin--
my $startTime = time;
readRef($segment);
my $prefix = substr($fileInput,0,-6)."-temp/".substr($fileInput,0,-6)."-part";
print "\nPreparing to begin...\n";
unlink <$prefix*>;
splitFasta($prefix);
my @pid = ();
my @files = <$prefix*>;
my @failedFiles = ();
open(fileOutput,">","$fileOutput") or die "Could not open $fileOutput";
print fileOutput "Sequence Accession\tOrganism\tStrain Name\tSegment\tSubtype\tHost\t5' End\t3'End\tSequence Length\tExpected Length\tTotal Mixed Bases\tError(s)\n";
close fileOutput;
my $countComplete = 0;
my $countProcesses = 0;
foreach my $f (@files) {
    $countProcesses++;
    while($countProcesses > $maxThreads){
        wait;
        $countComplete += ($? >> 8);
        $countProcesses--;
    }
    if (!$verbose) {
        print "\rCompleted $countComplete sequences\t        Using $countProcesses active threads    ";
    }
    my $pid = fork();
    if ($pid == 0){
        open(fileReference, "< $fileReference") or die "$fileReference failed to load";
        if ($verbose) {
            print "  Processing ./$f...        View progress at ./$f-muscle-log\n";
        }
        open(my $current, "<", $f) or die $!;
        open(my $combined, ">", "$f-combined") or die $!;
        while (<$current>){
            print $combined $_;
        }
        seek (fileReference,0,0);
		print $combined "\n";
        while (<fileReference>){
            print $combined $_;
        }
        close $current;
        close $combined;
        close fileReference;
        if ($^O =~ /linux/){
            system("./muscle-linux -in ./$f-combined -out ./$f-aligned -quiet -verbose -log ./$f-muscle-log");
        } elsif ($^O =~ /darwin/){
            system("./muscle-darwin -in ./$f-combined -out ./$f-aligned -quiet -verbose -log ./$f-muscle-log");
        } elsif ($^O =~ /Win32/){
            system("./muscle-win32.exe -in ./$f-combined -out ./$f-aligned -quiet -verbose -log ./$f-muscle-log");
        } elsif ($^O =~ /cygwin/){
            system("./muscle-cygwin.exe -in ./$f-combined -out ./$f-aligned -quiet -verbose -log ./$f-muscle-log");
        } else {
            die("Error: OS not recognized.");
        }
        my %currentSequenceSet = %{readInput("./$f-aligned")};
        chkExtra(\%currentSequenceSet);
        chkCTS(\%currentSequenceSet);
        chkInt(\%currentSequenceSet);
        printTSV(\%currentSequenceSet);
        if ($verbose) {
            print "  Completed processing ./$f\n";
        }
        my $completed = (keys %currentSequenceSet) - 1;
        exit ($completed);
    }
}

while (waitpid(-1,"WNOHANG") != -1){
    $countProcesses--;
    $countComplete += ($? >> 8);
    if (!$verbose) {
        print "\rCompleted $countComplete sequences\t        Using $countProcesses active threads      ";
    }
}

if ($boolCleanup) {
    print "\n\nCleaning up...\n";
    unlink <$prefix*>;
    rmdir (substr($fileInput,0,-6)."-temp");
}
my $endTime = time;
print "\n\nAll parts completed in ",$endTime - $startTime," seconds.\n\n";
print "Please review the results in $fileOutput.\n\n";
print "Created by Timmy Huang.\n\n";

#--Main End--



sub splitFasta {
    open(fileInput, "<","$fileInput") or die "$fileInput failed to load";
    my $countSequences = 0;
    while (my $line = <fileInput>){     
        if ($line =~ /^>/) {
            $countSequences++;
        }
    }
    print "\nTotal Sequences: $countSequences\n\n";
    my $digits = length($countSequences);
    mkdir substr($fileInput,0,-6)."-temp/";
    my $outputPrefix = $_[0];
    my $outputPart;
    seek (fileInput,0,0);
    my $countFiles = 0;
    $countSequences = 0;
    while (my $line = <fileInput>) {
        if ( $line =~ /^>/ ) {
            if ($countSequences % $maxSequencesPerAlignment == 0) {
                $countFiles++;
                if ($countSequences > 0) {
                    close($outputPart);
                }
                open($outputPart, sprintf("> %s%0${digits}d",$outputPrefix,$countFiles)) or die "Could not open:", sprintf("> %s%0${digits}d",$outputPrefix,$countFiles);
            }
            $countSequences++;
        }
        print $outputPart $line;
    }
    close($outputPart);
}

sub readInput {
    my $filePart = $_[0];
    open(filePart, "< $filePart") or die "$filePart failed to load";
    my $line = <filePart>;
    my $sequence = "";
    my %input;
    while ($line) {
        chomp $line;
        my $name;
        if ($line =~ /^[>\Z]/ ) {
            $name = $line;
            $sequence = "";
        }
        while (defined($line = <filePart>) && ($line !~ /^[>\Z]/)) {
            chomp $line;
            $sequence = $sequence.$line;
        }
        "" =~ /()/;
        $sequence =~ /^(\-+)/;
        my $offset5 = length($1);
        "" =~ /()/;
        $sequence =~ /(\-+)$/;
        my $offset3 = length($1);
        my $errorInt = {    
            insertion  =>  [], 
            deletion   =>  [],
            mixedBases =>  [],
        };
        $input{$name} = { 
            sequence    => $sequence, 
            offset5     => $offset5, 
            offset3     => $offset3, 
            bpLength    => length(my $temp = $sequence =~ s/-//gr),
            totalMixed  => 0,
            errorExtra  => [], 
            errorCTS    => [], 
            errorInt    => $errorInt
        };
    }
    close filePart;
    return \%input;
}

sub readRef {
	if (!$fileReference){
		$fileReference = $fileRef{$_[0]};
	}
    open(fileReference, "< $fileReference") or die "$fileReference failed to load";
    my $line = <fileReference>;
    chomp($line);
    $refName = $line;
    close fileReference;
}

sub chkExtra{
    my %currentSequenceSet = %{$_[0]};
    my $length5 = $currentSequenceSet{$refName}{offset5};
    my $length3 = $currentSequenceSet{$refName}{offset3};
    foreach (keys %currentSequenceSet) {
        my $extra5 = substr($currentSequenceSet{$_}{sequence}, 0, $length5);
        my $extra3 = substr($currentSequenceSet{$_}{sequence}, -$length3, $length3);
        my $extra5NoGaps = ($extra5 =~ s/\-//gr);
        my $extra3NoGaps = ($extra3 =~ s/\-//gr);
        if (length($extra5NoGaps) > 0 && $boolIRDCodesOnly){
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Ext5'");
        } elsif ($extra5NoGaps eq "T") {
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Taq5'");
        } elsif ($hoff5 =~ /$extra5NoGaps$/ && length($extra5NoGaps) >= $hoffThreshold ) {
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Pri5'");                 
        } elsif (length($extra5NoGaps) >= $BLASTthreshold){ 
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Clo5' (BLAST)");
        } elsif (length($extra5NoGaps) >= 1) { 
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Unk5'");
        }
        if (length($extra3NoGaps) > 0 && $boolIRDCodesOnly) {
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Ext3'");
        } elsif ($extra3NoGaps eq "A"){
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Taq3'");       
        } elsif ($hoff3 =~ /^$extra3NoGaps/ && length($extra3NoGaps) >= $hoffThreshold) {
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Pri3'");                 
        } elsif (length($extra3NoGaps) >= $BLASTthreshold){
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Clo3' (BLAST)");
        } elsif (length($extra3NoGaps) >= 1) {
            push (@{$currentSequenceSet{$_}{errorExtra}}, "NCR-Unk3'");
        }
    }
}

sub chkCTS {
    my %currentSequenceSet = %{$_[0]};
    my $refSeqNoDashes = substr($currentSequenceSet{$refName}{sequence},$currentSequenceSet{$refName}{offset5},length($currentSequenceSet{$refName}{sequence})-$currentSequenceSet{$refName}{offset5}-$currentSequenceSet{$refName}{offset3}) =~ s/\-//gr;
    my $seqNoGaps0 = substr($currentSequenceSet{$refName}{sequence},0,$currentSequenceSet{$refName}{offset5}).$refSeqNoDashes.(substr($currentSequenceSet{$refName}{sequence},-$currentSequenceSet{$refName}{offset3},$currentSequenceSet{$refName}{offset3}));
    foreach (keys %currentSequenceSet) {
        my $begin5 = ($currentSequenceSet{$refName}{offset5} >= $currentSequenceSet{$_}{offset5} ? $currentSequenceSet{$refName}{offset5} : $currentSequenceSet{$_}{offset5});
        my $end3 = ($currentSequenceSet{$refName}{offset3} >= $currentSequenceSet{$_}{offset3} ? $currentSequenceSet{$refName}{offset3} : $currentSequenceSet{$_}{offset3});
        my $seqNoGaps = substr($currentSequenceSet{$_}{sequence},0,$currentSequenceSet{$_}{offset5}).(substr($currentSequenceSet{$_}{sequence},$currentSequenceSet{$_}{offset5},length($currentSequenceSet{$_}{sequence})-$currentSequenceSet{$_}{offset5}-$currentSequenceSet{$_}{offset3}) =~ s/\-//gr).substr($currentSequenceSet{$_}{sequence},-$currentSequenceSet{$_}{offset3}, $currentSequenceSet{$_}{offset3});        
        if (($currentSequenceSet{$_}{offset5} < $lengthCTS5 + $currentSequenceSet{$refName}{offset5}) && (substr($seqNoGaps, $begin5, ($currentSequenceSet{$refName}{offset5} + $lengthCTS5)-$begin5) ne substr($seqNoGaps0, $begin5, ($currentSequenceSet{$refName}{offset5} + $lengthCTS5)-$begin5))){
            push (@{$currentSequenceSet{$_}{errorCTS}}, "NCR-CTS5'");
        }
        if (($currentSequenceSet{$_}{offset3} < $lengthCTS3 + $currentSequenceSet{$refName}{offset3}) && (substr($seqNoGaps, -($currentSequenceSet{$refName}{offset3}+$lengthCTS3), $currentSequenceSet{$refName}{offset3}+$lengthCTS3 - $end3)) ne substr($seqNoGaps0, -($currentSequenceSet{$refName}{offset3}+$lengthCTS3), $currentSequenceSet{$refName}{offset3}+$lengthCTS3 - $end3)){
            push (@{$currentSequenceSet{$_}{errorCTS}}, "NCR-CTS3'"); 
        }        
    }
}

sub chkInt {
    my %currentSequenceSet = %{$_[0]};
    foreach (keys %currentSequenceSet) {
        my $polyABegin = length($currentSequenceSet{$refName}{sequence}) - $currentSequenceSet{$refName}{offset3} - length($NCR3) + $posPolyA;
        my $polyAEnd = length($currentSequenceSet{$refName}{sequence}) - $currentSequenceSet{$refName}{offset3} - $lengthCTS3;
        my $refSeqLength = length($currentSequenceSet{$refName}{sequence}) - $currentSequenceSet{$refName}{offset5} - $currentSequenceSet{$refName}{offset3};
        for (my $i = $currentSequenceSet{$refName}{offset5} >= $currentSequenceSet{$_}{offset5} ? $currentSequenceSet{$refName}{offset5} : $currentSequenceSet{$_}{offset5}; $i<length($currentSequenceSet{$_}{sequence}) - ($currentSequenceSet{$refName}{offset3} >= $currentSequenceSet{$_}{offset3} ? $currentSequenceSet{$refName}{offset3} : $currentSequenceSet{$_}{offset3}); $i++){
            my $compare1 = substr($currentSequenceSet{$refName}{sequence}, $i, 1);
            my $compare2 = substr($currentSequenceSet{$_}{sequence},$i,1);
            if ($compare1 eq "-" && $compare1 ne $compare2){
                my $temp = "$i";
                my $start = $i;
                for (my $consecutive = 1; substr($currentSequenceSet{$refName}{sequence}, $i+1, 1) eq "-" && substr($currentSequenceSet{$refName}{sequence}, $i+1, 1) ne substr($currentSequenceSet{$_}{sequence}, $i+1,1); $consecutive++){
                    $i++;
                }
                if ($i - $start > 0) {
                    $temp .= "-".$i;
                }
                if (($i - $start) % 3 == 2){
                    $temp .= "(3x)";
                }
                if ($start > $polyABegin && $i < $polyAEnd && $compare2 eq "A" && $i - $start < 2){
                    $temp .= "(PolyA?)";
                }
                push(@{$currentSequenceSet{$_}{errorInt}{insertion}},$temp);
            } elsif ($compare2 eq "-" && $compare1 ne $compare2){
                my $temp = "$i";
                my $start = $i;
                for (my $consecutive = 1; substr($currentSequenceSet{$_}{sequence},$i+1,1) eq "-" && substr($currentSequenceSet{$refName}{sequence}, $i+1, 1) ne substr($currentSequenceSet{$_}{sequence}, $i+1,1); $consecutive++){
                    $i++;
                }
                if ($i - $start > 0) {
                    $temp .= "-".$i;
                }
                if (($i - $start) % 3 == 2){
                    $temp .= "(3x)";
                }
                if ($start > $polyABegin && $i < $polyAEnd && $compare1 eq "A" && $i - $start < 2){
                    $temp .= "(PolyA?)";
                }
                if (($start - $currentSequenceSet{$refName}{offset5}) / $refSeqLength > 0.02 && ($i - $currentSequenceSet{$refName}{offset5}) / $refSeqLength < 0.98){
                    if ($i - $start >= 20) {
                        $temp .= "(Orphan?)";
                    }
                } else {
                    if ($i - $start > 20){
                        $temp = "delete-orphan$temp";
                    } elsif ($i - $start >= 10) {
                        $temp .= "(Orphan?)";
                    }
                }
                push(@{$currentSequenceSet{$_}{errorInt}{deletion}},$temp);
            } elsif ($compare2 =~ /^[RYSWKMBDHVN]$/){
                my $temp = "$i";
                my $start = $i;
                my $mix = $compare2;
                for (my $consecutive = 1; substr($currentSequenceSet{$_}{sequence},$i+1,1) eq "$mix"; $consecutive++){
                    $currentSequenceSet{$_}{totalMixed}++;
                    $i++;
                }
                if ($i - $start > 0) {
                    $temp .= "-".$i;
                }
                $temp .= $mix;
                push(@{$currentSequenceSet{$_}{errorInt}{mixedBases}}, $temp);
                $currentSequenceSet{$_}{totalMixed}++;
            }
        }
    }
}

sub formatError{
    my @orphanStarts = ();
    my @orphanEnds = ();
    my $errorNCR5del = "";
    my $errorNCR5ins = "";
    my $errorCDSdel = "";
    my $errorCDSins = "";
    my $errorNCR3del = "";
    my $errorNCR3ins = "";
    my $errorNCR3delA = "";
    my $errorNCR3insA = "";
    my $errorNCR5mix = "";
    my $errorNCR3mix = "";
    my $errorCDSmix = "";
    my %current = %{$_[0]};
    my %currentSequenceSet = %{$_[1]};
    my @insertions = @{$current{errorInt}{insertion}};
    for (my $i = 0; $i < @insertions; $i++) {
        if ($insertions[$i] =~ /^(\d+)-*(\d*)/){
            if ($1 < length($NCR5) + $currentSequenceSet{$refName}{offset5}) {
                $errorNCR5ins .= ",".$insertions[$i];
            } elsif ($1 > length($currentSequenceSet{$refName}{sequence}) - length($NCR3) - $currentSequenceSet{$refName}{offset3} || $2 > length($currentSequenceSet{$refName}{sequence}) - length($NCR3) - $currentSequenceSet{$refName}{offset3}) {
                $errorNCR3ins .= ",".$insertions[$i];
            } else {
                $errorCDSins .= ",".$insertions[$i];
            }
        }
    }
    if (length($errorNCR5ins) > 0){
        $errorNCR5ins = substr($errorNCR5ins, 1);                       
        $errorNCR5ins = "\tNCR-ins5'($errorNCR5ins)";
    }
    if (length($errorNCR3ins) > 0){
        $errorNCR3ins = substr($errorNCR3ins, 1);                       
        $errorNCR3ins = "\tNCR-ins3'($errorNCR3ins)";
    }
    if (length($errorCDSins) > 0){
        $errorCDSins = substr($errorCDSins, 1);  
        $errorCDSins = "\tCDS-ins($errorCDSins)";
    }
    my @deletions = @{$current{errorInt}{deletion}};
    for (my $i = 0; $i < @deletions; $i++) {
        if ($deletions[$i] =~ /^(\d+)-*(\d*)/){
            if ($1 < length($NCR5) + $currentSequenceSet{$refName}{offset5}) {
                $errorNCR5del .= ",".$deletions[$i];
            } elsif ($1 > length($currentSequenceSet{$refName}{sequence}) - length($NCR3) - $currentSequenceSet{$refName}{offset3} || $2 > length($currentSequenceSet{$refName}{sequence}) - length($NCR3) - $currentSequenceSet{$refName}{offset3}) {
                $errorNCR3del .= ",".$deletions[$i];
            } else {
                $errorCDSdel .= ",".$deletions[$i];
            }
        }
        if ($deletions[$i] =~ /orphan/i && $deletions[$i] =~ /(\d+)-*(\d*)/){
            push(@orphanStarts, $1);
            push(@orphanEnds, $2);
        } 
    }
    if (length($errorNCR5del) > 0){
        $errorNCR5del = substr($errorNCR5del, 1);                    
        $errorNCR5del = "\tNCR-del5'($errorNCR5del)";
    }
    if (length($errorNCR3del) > 0){
        $errorNCR3del = substr($errorNCR3del, 1);                       
        $errorNCR3del = "\tNCR-del3'($errorNCR3del)";
    }
    if (length($errorCDSdel) > 0){
        $errorCDSdel = substr($errorCDSdel, 1);                        
        $errorCDSdel = "\tCDS-del($errorCDSdel)";
    }
    my @mixedBases = @{$current{errorInt}{mixedBases}};
    for (my $i = 0; $i < @mixedBases; $i++) {
        if ($mixedBases[$i] =~ /^(\d+)-*(\d*)/){
            if ($1 < length($NCR5) + $currentSequenceSet{$refName}{offset5}) {
                $errorNCR5mix .= ",".$mixedBases[$i];
            } elsif ($1 > length($currentSequenceSet{$refName}{sequence}) - length($NCR3) - $currentSequenceSet{$refName}{offset3} || $2 > length($currentSequenceSet{$refName}{sequence}) - length($NCR3) - $currentSequenceSet{$refName}{offset3}) {
                $errorNCR3mix .= ",".$mixedBases[$i];
            } else {
                $errorCDSmix .= ",".$mixedBases[$i];
            }
        }   
    }
    if (length($errorNCR5mix) > 0){
        $errorNCR5mix = substr($errorNCR5mix, 1);                    
        $errorNCR5mix = "\tNCR-mix5'($errorNCR5mix)";
    }
    if (length($errorNCR3mix) > 0){
        $errorNCR3mix = substr($errorNCR3mix, 1);                    
        $errorNCR3mix = "\tNCR-mix3'($errorNCR3mix)";
    }
    if (length($errorCDSmix) > 0){
        $errorCDSmix = substr($errorCDSmix, 1);                      
        $errorCDSmix = "\tCDS-mix($errorCDSmix)";
    }
    my $errorExtra = "";
    foreach (@{$current{errorExtra}}) {
        $errorExtra .= "\t$_";
    }
    my $errorCTS = "";
    foreach (@{$current{errorCTS}}) {
        $errorCTS .= "\t$_";
        if ($_ =~ /CTS5/) {
            foreach my $i (@orphanStarts) {
                if ($i < $lengthCTS5 + $current{offset5}){
                    $errorCTS .= "(Orphan?)";
                    last;
                }
            }
        }
        if ($_ =~ /CTS3/) {
            foreach my $i (@orphanEnds) {
                if ($i > length($current{sequence}) - $lengthCTS3 - $current{offset3}){
                    $errorCTS .= "(Orphan?)";
                    last;
                }
            }
        }
    }
    my $errorInt = $errorNCR5ins.$errorNCR5del.$errorCDSins.$errorCDSdel.$errorNCR3ins.$errorNCR3del;
    if ($current{totalMixed} > $currentSequenceSet{$refName}{bpLength} / 200) {
        $errorInt .= $errorNCR5mix;
        $errorInt .= $errorCDSmix;
        $errorInt .= $errorNCR3mix;
    }
    return "$errorExtra$errorCTS$errorInt";
}

sub printTSV { 
    my %currentSequenceSet = %{$_[0]};
    open(fileOutput,">>","$fileOutput") or die "Could not open $fileOutput";
    foreach (sort keys %currentSequenceSet){
        my $error = formatError($currentSequenceSet{$_},\%currentSequenceSet);
        if ($_ ne $refName && (length($error) > 0 || $boolShowCorrectSeq)){
			my $line = "";
            my @header = split(/[\|\:\>]/, $_);
            for (my $i = 2; $i < @header; $i += 2) {
                $line .= "$header[$i]\t";
            }
            $line .= substr($currentSequenceSet{$_}{sequence},0, $printEndLength)."\t";
            $line .= substr($currentSequenceSet{$_}{sequence}, -$printEndLength, $printEndLength)."\t";
            $line .= $currentSequenceSet{$_}{bpLength}."\t";
            $line .= $currentSequenceSet{$refName}{bpLength}."\t";
            $line .= $currentSequenceSet{$_}{totalMixed};
            $line .= $error;
            $line .= "\n";
            print fileOutput $line;
        }
    }
    close fileOutput;
} 

sub printHelp {
    die( 
    "\n",
    "----------------USAGE HELP----------------\n\n",
    "EXAMPLE USAGE:\tperl AacadFluAlign.pl -in PB1remainder-all.fasta -seg PB1 -clean\n\n",
    "Command Line Flags:\n\n",
    "    -in <file.fasta>\t\tName of input file\t\t\t\t\t\t\t(Mandatory)\n",
        "    -seg <gene segment>\t\tGene segment to input\t\t\t\t\t\t\t(Mandatory)\n",
        "    -out <file.tsv>\t\tName of output file\t\t\t\t\t\t\t(Default: <same-as-input.tsv>)\n",
        "    -maxSeqs <integer>\t\tMaximum sequences per MUSCLE alignment\t\t\t\t\t(Default: 1)\n",
        "    -maxThreads <integer>\tMaximum number of parallel threads used in analysis\t\t\t(Default: 32)\n",
        "    -clean\t\t\tCleanup temporary alignment files after completing \t\t\t(Default: FALSE)\n",
        "    -verbose\t\t\tDisplay detailed runtime processes\t\t\t\t\t(Default: FALSE)\n",
        "    -ref <reference.fasta>\tName of reference sequence\t\t\t\t\t\t(Default: Automatic)\n",
        "    -showAll\t\t\tOutput all sequences, not only flagged sequences\t\t\t(Default: FALSE)\n",
        "    -IRD\t\t\tUse IRD error codes\t\t\t\t\t\t\t(Default: FALSE)\n",
        "    -hoffMatchMin <integer>\tMinimum number of base pair matches to be classified as Hoffman Primer\t(Default: 5)\n",
        "    -blastMatchMin <integer>\tMinimum size of extraneous sequence to be classified as cloning plasmid\t(Default: 12)\n".
        "    -printEndLength <integer>\tNumber of base pairs to print on 3' and 5' ends\t\t\t\t(Default: 50)\n",
        "-----------------------------------------\n\n",
         "Email ",'Timmy Huang (th174@duke.edu)'," with any further questions.\n\n");
}

