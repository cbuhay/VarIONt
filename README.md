## The varIONT Single Nucleotide Variant Caller ##

VarIONt is a SNV caller that is trained initially on Ion Torrent PGM data, but is versatile and extensible to be used in other next-gen sequencing platforms like Illumina.  It’s strength stems from a multitude of user options to be able to fine tune parameters.  Options like minimum coverage and minor allele fraction are adjustable.  Moreover, the concept of “minimum usable coverage” can also be adjusted.  Minimum usable coverage is the coverage of the reference allele and the primary variant allele.  In sequencing applications that employ very deep coverage (>2000X), the minimum usable coverage is advantageous because it ignores much of the noise stemming from random error.  

## Usage ##

	$ perl varIONt.pl -i <mpileup file> -o <output format> -m <minumum coverage> -u <usable coverage> -f <alelle frequency> -t <minumum allele fraction> -v

-i: input file in pileup format. REQUIRED
-o: output choices- vcf, pileup, both.  Default: VCF.
-v: only outputs variants; or else outputs all sites.
-m: minimum coverage. Default: 20.
-u: usable coverage (reference plus primary alternate allele).  Default: 30.
-f: minor allele fraction. Default 0.05 (or 5%).
-t: min recorded allele fraction. Default 0.03 (or 3%).

Example:
	$ perl varIONt.pl -i mypileupfile -o vcf -m 10 -u 20 -f 0.10 -t 0.05

The output will be in VCF format.
The miniumum read coverage surveyed will be 10X, all sites <10X will not be considered.
The minimum usable coverage is 20-- the coverage of major + most prevalent minor allele must be ≥20.
The minor allele fraction which to call a PASS variant is 10%
The alllele fraction which to consider ANY variant is 5%

## Output ##

Able to output a VCF file, pileup-style file, or both.
