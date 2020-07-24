#!/usr/bin/perl

##############################################################################
# Script takes list of chr/positions, annotates with rate and k-mer according
# to specified file
#
# Build under Perl v5.18.2
##############################################################################

use strict;
use warnings;
use POSIX;
use FindBin;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# my $parentdir="$FindBin::Bin/../";
my $fastadir="/usr/local/share/data";
my $assetdir="/var/www/jedidiahcarlson.com/assets";

my $help=0;
my $man=0;
# my $chr;
my $adj=3;
my $f_fasta = "$fastadir/human_g1k_v37.fasta";
my $f_positions = "$assetdir/test_sites.txt";
my $seqflag;
my $scale=0;
my $sciflag;
my $subseq1 = $adj*2+1;
my $rates = "$assetdir/ERV_${subseq1}bp_rates.txt";

GetOptions (
# 'chr=i'=> \$chr,
'adj=i' => \$adj,
'rates=s' => \$rates,
'ref=s' => \$f_fasta,
'in=s' => \$f_positions,
'seq' => \$seqflag,
'scale=s' => \$scale,
'sci' => \$sciflag,
'help|?'=> \$help,
'ref=s' => \$f_fasta,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

my $subseq = $adj*2+1;
$rates =~ s/$subseq1/$subseq/ee; #"$parentdir/assets/ERV_${subseq}bp_rates.txt";

# Target mutation rate
my $factor;
if($scale!=0){
  my $mu = $scale;
  my $nbases = 3e9;
  my $indmuts = $mu*$nbases;
  my $erv_count = 36087319;
  $factor = $indmuts/$erv_count;
  if($sciflag){
    $factor=$factor*1e-8;
  }
} else {
  $factor=1;
}

my %cathash = (
  AT_CG => 0,
  AT_GC => 1,
  AT_TA => 2,
  GC_AT => 3,
  GC_CG => 4,
  GC_TA => 5,
);

# Initialize and hash rate table
open my $rateFH, '<', $rates or die "can't open $rates: $!";
readline($rateFH); #<-throws out header

my %hash=();
while (<$rateFH>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[0];
	# print "$key\n";
	$hash{$key}=[@line[1 .. $#line]];
}
close $rates;

my @sites;
my $prevchr = '';
my $seq;
my $altseq;

# my $firstline = readline($positions);
# chomp($firstline);
# $firstline =~ tr/\x{d}\x{a}//d;
# my @line=split(/\t/, $linestr);
# my $sitechr=$line[0];
# if($sitechr =~ /^C/){
#   print "$firstline\tCATEGORY\tMU";
#   if($seq){
#     print "\tMOTIF\n";
#   }
#   next;
# };

open my $positions, '<', $f_positions or die "can't open $f_positions: $!";
while(<$positions>){
	# print "$_\n";
	chomp;
	my $linestr=$_;
  $linestr =~ tr/\x{d}\x{a}//d;
	my @line=split(/\t/, $linestr);
	my $sitechr=$line[0];
  if($.==1 && $sitechr =~ /^C/){
    if($seqflag){
      print "$linestr\tCATEGORY\tMU\tMOTIF\n";
    } else {
      print "$linestr\tCATEGORY\tMU\n";
    }
    next;
  } elsif($.==1 && $sitechr !~ /^C/){
    print "CHR\tPOS\tREF\tALT\t";
    if($#line>3){
      for my $i(4..($#line-1)){
        my $colnum=$i+1;
        print "V$colnum\t";
      }
      my $lastcol = $#line+1;
      print "V$lastcol\t";
    }
    print "CATEGORY\tMU";
    if($seqflag){
      print "\tMOTIF\n";
    } else {
      print "\n";
    }
  }

	my $pos=$line[1];
  my $ref=$line[2];
  my $alt=$line[3];
  my $categ=&getCat($ref, $alt);

	# my $categ=$line[5];
	my $catind=$cathash{$categ};
  # print "$categ\t$catind\n";

	if ($sitechr ne $prevchr) {
      $seq=&getRef($sitechr);
	}

	my $base=substr($seq, $pos, 1);

	my $localseq = substr($seq, $pos-$adj-1, $subseq);
	if($localseq!~/[MNSW]/){
		# my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);
    my $altlocalseq = reverse $localseq;
    $altlocalseq  =~ tr/ACGT/TGCA/;
		my $sequence;
		if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
			$sequence = $localseq . '(' . $altlocalseq . ')';
		} else {
			$sequence = $altlocalseq . '(' . $localseq . ')';
		}

      my $mu = 'NA';
      if (exists $hash{$sequence}) { # Ref sequence may contain additional codes such as 'R'. 
         $mu = $factor*$hash{$sequence}[$catind];
      }

		# print OUT "$chr\t$i\t$hash{$sequence}\n";
		if($seqflag){
			# print OUT "$linestr\t$sequence\t$hash{$sequence}[$catind]\t\n";
      print "$linestr\t$categ\t$mu\t$sequence\n";
		}	else {
      print "$linestr\t$categ\t$mu\n";
		}
	}

	$prevchr=$sitechr;
}

close $positions;


##############################################################################
# Subroutine for getting mutation type from REF/ALT fields
##############################################################################
sub getCat{
  my $ref=shift;
  my $alt=shift;

  my $cat="${ref}${alt}";
  my $categ;
  if($cat =~ /AC|TG/){
    $categ = "AT_CG";
  } elsif($cat =~ /AG|TC/){
    $categ = "AT_GC";
  } elsif($cat =~ /AT|TA/){
    $categ = "AT_TA";
  } elsif($cat =~ /GA|CT/){
    $categ = "GC_AT";
  } elsif($cat =~ /GC|CG/){
    $categ = "GC_CG";
  } elsif($cat =~ /GT|CA/){
    $categ = "GC_TA";
  }

  return $categ;
}

##############################################################################
# Subroutine for getting reference sequence
##############################################################################
sub getRef{
	my $chr=shift;
	my $nextchr;
	my $hasprefix=($chr =~ /^chr/);

   $chr =~ s/chr//g;
   if ($chr<22) {
      $nextchr=$chr+1;
   } elsif ($chr==22) {
      $nextchr="X";
   }

   if ($hasprefix) {
      $chr = "chr" . $chr;
      $nextchr = "chr" . $nextchr;
   }

	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";
	# print "Getting reference sequence for chromosome $chr...\n";

	my $seq;
	while (<$fasta>) {
		chomp;
		if (/>$chr /../>$nextchr /) {
			next if />$chr / || />$nextchr /;
			$seq .=$_;
		}
	}

	return $seq;
	close $fasta;
}

__END__
=head1 NAME

mr_eel.pl - Annotation Utility for [M]utation [R]ate [E]stimation using [E]xtremely rare variants and [L]ocal sequence context

=head1 SYNOPSIS

        anno_rate.pl [OPTIONS]
        Options:
		--help			program documentation
    --in        /path/to/inputfile.txt
    --seq       print motif
		--rates		/path/to/{K}bp_rates.txt.
		--ref			/path/to/human_g1k_v37.fasta


=head1 OPTIONS

=over 8

=item B<--help>

Display this documentation

=item B<--in>

File to be processed. A sample file is mr_eel/assets/test_sites.txt. Positions in file should be mapped to the GRCh37 reference genome, and the file should be formatted as described at http://www.jedidiahcarlson.com/mr-eel/.

=item B<--seq>

This flag prints an additional column displaying the K-mer motif of each position

=item B<--rates>

Path to rate table. Default is mr_eel/assets/ERV_7bp_rates.txt

=item B<--ref>

Path to GRCh37 reference genome.  Default is mr_eel/assets. You will either need to download one to the /assets/ folder using wget or the included script (mr_eel/bin/download_ref.pl), or specify a custom location. This utility was built using the human_g1k_v37 fasta file from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/. If using a different assembly for GRCh37, please ensure that the the heading lines start with the chromosome number, and not "chr" (e.g., ">1", not ">chr1" or ">chrom1"). This requirement will be removed in a future release.

=back

=head1 DESCRIPTION

B<mr_eel.pl> - Annotation Utility for [M]utation [R]ate [E]stimation using [E]xtremely rare variants and [L]ocal sequence context

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Bioinformatics E<10> University of Michigan

=cut
