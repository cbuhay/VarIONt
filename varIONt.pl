#!/usr/bin/perl -w
#
# Author: Christian Buhay
#
# Description: This is run after the mpileup
#  is run.  At this point, you need to change
#  the adjustable features.
#############################################

# Indicated lines added by Ian Gibson ('IBG')

use warnings;    # IBG
use strict;      # IBG
use Getopt::Std; 

### GLOBAL VARS ###
###################
my %fre;
my %opt;

### OPTIONS ###
###############
getopts("i:o:m:u:f:v:t",\%opt);
my $infile = $opt{i} || &USAGE;
# second invocation argument should be one of 'vcf', 'both' or 'mpileup'  # IBG
my $output_choice = $opt{o} || 'VCF'; # IBG, default output to VCF only
my $min_total_coverage  = $opt{m} || 20;      #default 20
my $min_usable_coverage = $opt{u} || 30;      #default: 30
my $min_var_freq        = $opt{f} || 0.05;    #default: 0.05
my $variant_output      = $opt{v} || 0;
my $fraction_filter    = $opt{t} || 0.03;       #default: 0.03
# Parse input file to create output file.
my $input_file_no_extension = $infile;   # IBG
$input_file_no_extension =~ s{\.[^.]+$}{};   # IBG, remove input file extension
my $outfile = "$input_file_no_extension.variont.mpileup";   # IBG, modified
#my $outfile2 = "$input_file_no_extension.ref.mpileup";

### MAIN CODE ###
#################
open( FIN, "$infile" )   || die "Can't open $infile: $!\n";
	
my @lines;   # IBG

while(<FIN>)
{
    chomp; 
    my $line = $_; 
    my ($chr, $coord, $ref, $total_coverage, $pileup_string) = split(/\s+/, $line);

    #Filter out indels and/or zero consensus quality
    if ($ref eq "*") {next;}

    #Filter out reference bases that are "N"
    if ($ref eq "N") {next;}

    #Filter out zero consensus quality & zero SNP quality
    #Hopefully, this lets in 100% homozygous vars
    #if ($consensus_qual eq "0" && $snp_qual eq "0") {next;}

    #Filter out total-cov below threshold
    #This includes ref, var allele, and noise
    if ($total_coverage < $min_total_coverage) {next;}

    my $pileup_analysis_results = &PARSE_PILEUP($ref, $pileup_string);
    my ($var_allele, $usable_cov, $var_allele_freq, $Acall, $Tcall, $Ccall, $Gcall, $VCFaddition) = split(/\t/, $pileup_analysis_results);
    my $base_output = "$var_allele_freq\t$total_coverage\t$pileup_string";
    my $ancillary_output = "$Acall\t$Tcall\t$Ccall\t$Gcall\t$VCFaddition";

    # IBG - changed all 'print OUT' lines to 'push (@lines, ...)':
    if ( $var_allele eq $ref )
    {
        $var_allele = $ref;
        push (@lines, "$chr\t$coord\t$ref\t$var_allele\t$base_output\t$ancillary_output\tref_match\n");
    } 
    elsif ( $usable_cov < $min_usable_coverage && $var_allele_freq >= $min_var_freq )
    {
        push (@lines, "$chr\t$coord\t$ref\t$var_allele\t$base_output\t$ancillary_output\tambiguous_var_lowcov\n");
    } 
    elsif ( $usable_cov >= $min_usable_coverage && $var_allele_freq >= $min_var_freq )
    {
        push (@lines, "$chr\t$coord\t$ref\t$var_allele\t$base_output\t$ancillary_output\tvariant\n");
    }
    elsif ( $usable_cov >= $min_usable_coverage && $var_allele_freq < $min_var_freq && $var_allele_freq >= $fraction_filter )
    {
    push (@lines, "$chr\t$coord\t$ref\t$var_allele\t$base_output\t$ancillary_output\tambiguous_var_lowfreq\n");
    }
    elsif ( $var_allele_freq < $fraction_filter )
    {
    push (@lines, "$chr\t$coord\t$ref\t$var_allele\t$base_output\t$ancillary_output\tref_match\n");
    }
}

# <IBG start>
if ( $output_choice eq 'both' ) {
    open( OUT, ">$outfile" ) || die "Can't open $outfile: $!\n";
    #open( OUT2, ">$outfile2" ) || die "Can't open $outfile: $!\n";
    if ( $variant_output ) {
	foreach (@lines) {
	    if ( $_ =~ /var/ ) {print OUT "$_";}
            #if ( $_ =~ /ref_match/ ) {print OUT2 "$_";}
	    }
    } else {
	foreach (@lines) {print OUT "$_";}
    }
    vcf_output();
}
elsif ( $output_choice eq 'mpileup' ) {
    open( OUT, ">$outfile" ) || die "Can't open $outfile: $!\n";
    #open( OUT2, ">$outfile2" ) || die "Can't open $outfile: $!\n";
    if ( $variant_output ) {
	    foreach (@lines) {
		if ( $_ =~ /var/ ) {print OUT "$_";}
                #if ( $_ =~ /ref_match/ ) {print OUT2 "$_";}
	    }
    } else {
	foreach (@lines) {print OUT "$_";}
    }
}
else { 
    print "Producing VCF-only file output. Re-run with second invocation argument ",
    "'mpileup' to output only the mpileup filtered results, or with 'both' to output ",
    "both formats in separate files.\n";
    vcf_output(); 
}

close(FIN);
close(OUT);
#close(OUT2);

# Read in the annotated output produced earlier, then output in VCF format.
sub vcf_output {
    # Open filhandles:
    open my $vcf_outfile, '>', "$input_file_no_extension-mpileup.vcf"
      or die "Cannot open VCF output file: $!\n";

    # Declare/initialize variables:
    my (
        $chr, $pos, $ref, $alt, $maf,
        $dp,  $pu, $geno, $ac, $tc, $cc, $gc, $st
    );
    my $id     = '.';
    my $qual =  '.';
    my $mmq =  '.';
    my $filter;
    my $format = "";

    # VCF-complying format, which doesn't include the pile-up string, since
    # the pile-up string contains illegal metacharacters:
    print $vcf_outfile '##fileformat=VCF4.1', "\n";
    print $vcf_outfile '##fileDate=', eval { `date +"%Y%m%d"` };   # not portable to non-Unix
    print $vcf_outfile '##generatedBy:', "$0", "\n";   # output program name to file
    print $vcf_outfile '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Fraction">', "\n";
    print $vcf_outfile '##INFO=<ID=ST,Number=.,Type=String,Description="Status">', "\n";
    print $vcf_outfile '##INFO=<ID=MMQ,Number=1,Type=Integer,Description="Maximum Mapping Quality">', "\n";
    print $vcf_outfile '##INFO=<ID=DP,Number=1,Type=Integer,Description="Reads Covering Site">', "\n";
    print $vcf_outfile '##INFO=<ID=A,Number=1,Type=Integer,Description="A-coverage">', "\n";
    print $vcf_outfile '##INFO=<ID=T,Number=1,Type=Integer,Description="T-coverage">', "\n";
    print $vcf_outfile '##INFO=<ID=C,Number=1,Type=Integer,Description="C-coverage">', "\n";
    print $vcf_outfile '##INFO=<ID=G,Number=1,Type=Integer,Description="G-coverage">', "\n";
    print $vcf_outfile '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "\n";
    print $vcf_outfile '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">', "\n";
    print $vcf_outfile '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', "\n";
    print $vcf_outfile '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference Read Depth">', "\n";
    print $vcf_outfile '##FORMAT=<ID=VR,Number=1,Type=Integer,Description="Major Variant Read Depth">', "\n";
    print $vcf_outfile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$input_file_no_extension\n";

    # Print columns in correct format:
    foreach (@lines) {
        (
            $chr, $pos, $ref, $alt, $maf,
            $dp,  $pu, $ac, $tc, $cc, $gc, $geno, $st
        ) = split(/\s+/);

        $alt = '.' if $ref eq $alt;    # ALT is '.' if identical to REF

	#CB: populate the "FILTER" field based on status
	if ($st eq "variant") {$filter = "PASS";
	} elsif ($st eq "ambiguous_lowcov") {$filter = "low_coverage";
	} elsif ($st eq "ambiguous_lowfreq") {$filter = "low_freq";
	} elsif ($st eq "ref_match") {$filter = ".";}

	#CB: provide a mechanism for only-variants output
	if ($variant_output) {
	    if ($_ =~ /ref_match/) {next;}
	}

        # VCF-complying output format:
        print $vcf_outfile "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t";
        print $vcf_outfile "MAF=$maf;ST=$st;MMQ=$mmq;DP=$dp;A=$ac;T=$tc;C=$cc;G=$gc\t";
	print $vcf_outfile "GT:VR:RR:DP:GQ\t$geno\n";
    }

    close $vcf_outfile;
}
# <IBG end>


### SUBROUTINES ###
###################
sub USAGE {
    print "\nUSAGE: $0 -i <pileup-input> -[omuf]\n";
    print "  -i: input file in pileup format. REQUIRED\n";
    print "  -o: output choices- vcf, pileup, both.  Default: VCF.\n";
    print "  -v: only outputs variants.\n";
    print "  -m: minimum coverage. Default: 20.\n";
    print "  -u: usable coverage (reference plus primary alternate allele).  Default: 30.\n";
    print "  -f: minor allele fraction. Default 0.05 (or 5%).\n\n";
    print "  -t: fraction filter. Default 0.03 (or 3%).\n\n";
    exit;
}

sub PARSE_PILEUP {
    my $ref_allele = shift;
    my $string = shift;
    $string =~ s/.\$//g; # remove variant calls at end of read
    #$string =~ s/(\-|\+){1}\w+//g; # remove indels in pileup
    $string =~ s/\^.{1}//g; # remove ascii character after ^ which is mapping qual
    $string =~ s/\*//g; # removes any pads made from deletions
    $string =~ tr/a-z/A-Z/; # make all letters uppercase
    my $foo = $string;
    $foo =~ s/[^0-9]/ /g;
    my @arrayA = split(/\s+/, $foo);
    my $arraylast = scalar(@arrayA);
    #print "$arraylast\n";
    #print "@arrayA\n";
    for (my $j=0; $j<$arraylast; $j++) {
        $string =~ s/(\+|\-)($arrayA[$j])[ACGTNacgtn]{$arrayA[$j]}//g; #removes indels
    }
    #print "$foo\n";
    #print "$string\n";
    my @array = split(//, $string);
    my $tot_cov = @array; #turns the array into a scalar
    
    my %histogram;
    $histogram{$_}++ for @array;
    
    if ( ! $histogram{'.'} ) {
        $histogram{'.'} = '0';
    } 
    if ( ! $histogram{','} ) {
        $histogram{','} = '0';  
    }
    my $ref_sum = $histogram{','} + $histogram{'.'};
    delete $histogram{','};
    delete $histogram{'.'};

    my @variant_allele_array = sort {$histogram{$b} <=> $histogram{$a}} keys %histogram; #sorts the hash largest value to smallest
    my $variant_allele = $variant_allele_array[0];
    if ( ! $variant_allele || $variant_allele !~ /(A|T|C|G)/ ) {
        $variant_allele = $ref_allele;
    }
    #the problem here is that if there are two variants with same size, it picks the alphabetically first.
    
    # Add in the reference; give undefined bases 0s.
    $histogram{$ref_allele} = $ref_sum;
    my @allele_bases = ("A","C","T","G");
    for my $base (@allele_bases) {
        if ( ! $histogram{$base} ) {
            $histogram{$base} = '0'; 
        }
    }

    my $usable_sum = $ref_sum;
    my $variant_allele_freq = 0;
    if ( $ref_allele ne $variant_allele ) {
        $usable_sum = $histogram{$ref_allele} + $histogram{$variant_allele};
        $variant_allele_freq = $histogram{$variant_allele} / $usable_sum;
    } 

    #CB: addendum for VCF genotype string
    #first, suss out the genotype
    my $genotype;
    my $geno_var_cov;
    if ( $variant_allele_freq < .01 ) {
	$genotype = "0\/0";
	$geno_var_cov = ".";
    } elsif ($variant_allele_freq > .90) {
	$genotype = "1\/1";
	$geno_var_cov = $histogram{$variant_allele};
    } else {
	$genotype = "0\/1";
	$geno_var_cov = $histogram{$variant_allele};
    }
    #then, recreate GT:VR:RR:DP:GQ; genotype:var_cov:ref_cov:read_depth:genotype_qual
    #I'm making a conscious decision to use total coverage vs. usable coverage
    my $VCF_string = "$genotype:$geno_var_cov:$ref_sum:$tot_cov:.";

    $variant_allele_freq = sprintf "%.3f", $variant_allele_freq; #round to 3 decimals.
    #print "$tot_cov\t$ref_allele\t$variant_allele\t$variant_allele_freq\t$histogram{A}\t$histogram{T}\t$histogram{C}\t$histogram{G}\n";
    my $pileup_parse = "$variant_allele\t$usable_sum\t$variant_allele_freq\t$histogram{A}\t$histogram{T}\t$histogram{C}\t$histogram{G}\t$VCF_string";
    return $pileup_parse;
}
