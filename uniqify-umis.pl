#!/usr/bin/env perl

# Produces a bam file with just the unique reads to unique.bam (with
# multiplicity in tag mu:i:N). Duplicates (defined as those having
# identical chromosome, read1 position, read2 position, read1 umi AND
# read2 umi are) are written to file 'umi-duplicates.bam'; improper
# (partially or fully unmapped and discordantly mapped genes) to
# 'improper.bam', multimappers to multimappers.bam
# Written by Philip Lijnzaad <plijnzaad@gmail.com>

use strict;

my $script="uniquify-umis.pl";
my $version="v2";

my @argv_copy=@ARGV;

use FileHandle;
use Number::Format;
use Getopt::Std;

use lib "$ENV{HOME}/git/demultiplex"; # find it at https://bitbucket.org/princessmaximacenter/demultiplex/

use mismatch;

use vars qw($opt_h $opt_u $opt_m $opt_o $opt_c);

$opt_m=1;                               # allow one mismatch

my $usage=" $script [ -c chromosome ] -u barcodes [ -m 1 ] [ -o dir ] [ file.bam ]    > uniq.sam";

## -c filters by chromosome so we can run in parallel

if ( !getopts("u:c:m:ho:") || ! $opt_u ||  $opt_h ) {
    die $usage; 
}

my $barcodes_mixedcase = mismatch::readbarcodes($opt_u); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'

my $mismatch_REs=undef;

if ($opt_m) { 
  $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, 
                                                allowed_mismatches =>$opt_m);
}
$barcodes_mixedcase=undef;              # not used in remainder, delete to avoid confusion


my $nexact=0;
my $nmismatched=0;
my $nunknown=0;

$opt_o = '.' unless $opt_o;
die "directory $opt_o does not exist" unless -d $opt_o;

my $outnames={MULT => "$opt_o/multimappers.bam",
              IMPR => "$opt_o/improper.bam", ## one or both unmapped or discordant
              UNIQBAM => "$opt_o/unique.bam", ## first readname is arbitrarily 'the unique one'
              DUPL => "$opt_o/umi-duplicates.bam"## the second (and third etc.) arbitrarily called duplicate
}; 

my $outfiles={};
for my $key (keys %$outnames) { 
  $outfiles->{$key}=open_bamout($key, $outnames->{$key});
}

my $fmt=new Number::Format(-thousands_sep => ',');
sub commafy {   $fmt->format_number($_[0]); }

my $samview="samtools view -h";         # alternatively, specify full path, or e.g. sambamba

my $flags = ['PE',                      # 0
             'proper pair',              # 1
             'unmapped',                 # 2
             'mate unmapped',            # 3
             'reverse strand',           # 4
             'mate on reverse strand',   # 5
             'read1', # i.e. first of two halfpairs # 6
             'read2', # i.e. second of two halfpairs # 7
             'part of a secondary alignment', # 8
             'not passing QC',                # 9
             'PCR or optical duplicate',      # 10
             'chimeric/supplementary alignment']; # 11

## make code slightly more readable:
my($FPE, $Fproper, $Funmapped, $Fmateunmapped, $Frev, $Fmaterev, 
   $Fread1, $Fread2, $Fsecondary, $FnopassQC, $Fduplicate, $Fchimeric ) =
    map { 1<<$_ } 0..11;

my $stats= {};

my $fh;
if (@ARGV) { 
  die "Only single file allowed " unless @ARGV==1;
  my $f=$ARGV[0];
  if ($f =~ /\.(b|cr)am$/) {
    my $cmd="$samview $f |";
    $fh = FileHandle->new($cmd) or die "$cmd: $!";
  } else {
    $fh = FileHandle->new("< $f") or die "$f: $!";
  }
} else {                                # stdin
    $fh = FileHandle->new_from_fd(0, "<") or die "stdin: $!";
}

## -- main --
$_=<$fh>;

die "sam/bam file must be name sorted by name first, sorry" unless $_ =~ /^\@HD\s+VN:\d+\.\d+\s+SO:queryname/;

foreach my $key (keys %$outfiles) {   $outfiles->{$key}->{fh}->print($_);}

my $seenpg=0;
my $seenumis=0;
my $seenhisat2=0;

my $hash={};
my $uniqbam={};

LINE:
while(<$fh>) { 
  if  ( /^\@/  ) {
    if(/^\@PG/ && ! $seenpg++) { 
      my $line="\@PG\tID:$script\tPN:$script\tVN:$version\tCL:$script @argv_copy\n" ;
      foreach my $key (keys %$outfiles) { $outfiles->{$key}->{fh}->print($line);}
## NOTE: this prepends the line, should be appended after the last @PG line ... change this
    }
    foreach my $key (keys %$outfiles) { $outfiles->{$key}->{fh}->print($_);}
    $seenhisat2++ if /^\@PG.*PN:hisat2/;
    $seenumis++ if /^\@PG.*ID:addumis2bam/;
    next LINE;
  }

  die "Expected hisat2 input, didn't find it in header ... " unless $seenhisat2;
  die "Expected addumis2bam in header ... " unless $seenumis;


  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
     $seq, $qual, @optionals)=split("\t", $_);
  my($NH, $u1, $u2)=@optionals[(-3,-2,-1)];

  next LINE if $opt_c && ($rname ne $opt_c);

  $stats->{'total'}++;

## exclude things:
  if(!   ($flag & $FPE)
     || !($flag & $Fproper) 
     || ($flag & $Funmapped)
     || ($flag & $Fmateunmapped)
     || ( $rnext ne '=' && $rname ne $rnext) 
      ) { 
    $stats->{'improper'}++;
    $outfiles->{IMPR}->{fh}->print($_);
    next LINE;
  };

  ## following depends on it being Hisat2 output ...
  if ( $NH ne "NH:i:1" ) {  # also test ($flags & $Fsecondary) ?
    $stats->{'multimappers'}++;
    $outfiles->{MULT}->{fh}->print($_);
    next LINE;
  };

  my $line1 = $_;
  $stats->{'concordant'} ++ ;

  my $line2=<$fh>;  

  my($qname2,$flag2, $rname2, $pos2, $mapq2, $cigar2, $rnext2, $pnext2, $tlen2,
     $seq2, $qual2, @optionals2)=split("\t", $line2);
  my($NH2, $u1_2, $u2_2)=@optionals2[(-3,-2,-1)];

  next LINE if $opt_c && ($rname2 ne $opt_c);

  $stats->{'total'} ++;

  die "not proper mates around line $.\n$line1\nversus\n$line2\n" unless
      $qname eq $qname2 
      && $pos == $pnext2 
      && $pos2 == $pnext
      && $rname eq $rname2
      && $NH ==$NH2
      && $u1 ==$u1_2
      && $u2 ==$u2_2;
  $stats->{'concordant'} ++ ;
  
  $u1 =~ s/u[12]:Z://;
  $u2 =~ s/u[12]:Z://;
  $u2 =~ s/[\n\r]$//;

  my $nrescued=0;
  if (! exists $barcodes->{$u1} ) {
    if ( $opt_m ) { 
      $u1=mismatch::rescue($u1, $mismatch_REs);
      $nrescued += defined($u1);
    } else {
      $u1=undef;
    }
  }

  if (! exists $barcodes->{$u2} ) {
    if ( $opt_m ) { 
      $u2=mismatch::rescue($u2, $mismatch_REs);
      $nrescued += defined($u2);
    } else {
      $u2=undef;
    }
  }
  my $goodumis=(defined( $u1) + defined( $u2));
  if ($goodumis < 2) { 
    $stats->{'lost umi(s)'} ++;
    $stats->{'    lost one'} += ($goodumis==1);
    $stats->{'    lost both'} += ($goodumis==0);
    next LINE;
  }
  $stats->{'good umis'} ++ ;
  $stats->{'    rescued umis'} += $nrescued;

  ## Counting starts here.  We should take the 'leftmost' read as the
  ## first one (regardless of whether it's read1 or read2)
  my ($combpos, $combumi);
  if ( $pos < $pos2) {
    $combpos = $pos.$pos2;
    $combumi = $u1.$u2;
  }  else {
    $combpos = $pos2.$pos;
    $combumi = $u2.$u1;
  }
  $combpos .= $rname unless $opt_c;     # r1 & r2 are on same chromo since we filter for it

  $hash->{$combpos}->{$combumi}++;
  if( $hash->{$combpos}->{$combumi} == 1  ) { 
    $uniqbam->{$combpos}->{$combumi} = "$line1$line2"; # print later with the actual multiplicity
  } else {
    $stats->{'umi-duplicates'}++;
    $outfiles->{DUPL}->{fh}->print($line1);
    $outfiles->{DUPL}->{fh}->print($line2);
  }
}                                       # while

my $nunique=0;
for my $combpos ( keys  %$hash ) { 
  $nunique += int keys( %{$hash->{$combpos}} );
}

$stats->{'umi-uniques'} = $nunique;

$fh->close();
print_stats($stats);
print_uniqbam($hash, $uniqbam);

foreach my $key (keys %$outfiles) {   close_bamout($outfiles->{$key}); } 

exit 0;

## ------------------------------------------------------------------------

sub print_uniqbam {
  my($hash, $uniqbam)=@_;

  for my $combpos (keys %{$hash}) {
    for my $combumi (keys %{$hash->{$combpos}}) {
      my $multip = $hash->{$combpos}->{$combumi};
      my $bamlines=$uniqbam->{$combpos}->{$combumi};
      $bamlines =~ s/\n/\tmu:i:$multip\n/mg;
      $outfiles->{UNIQBAM}->{fh}->print($bamlines);
    }
  }
}                                       # print_uniqbam

sub print_stats { 
  my ($stats)=@_;

  my $order=[
    'total',
    'improper',
    'multimappers',
    'concordant',
    'good umis',
    '    rescued umis',
    'lost umi(s)',
    '    lost one',
    '    lost both',
    'umi-uniques',
    'umi-duplicates'
      ];
  
  my $ntot=$stats->{'total'};
  print STDERR "# $version. All numbers are per mate, not per readpair?\n";
  for my $key (@$order) { 
    print STDERR join("\t", ($key, commafy($stats->{$key}), 
                             sprintf("%.1f%%", 100*$stats->{$key}/$ntot)))."\n";
  }
}                                       # print_stats

sub open_bamout { 
  my($key, $file)=@_;
  warn "opening $file ...\n";
  my $fh=FileHandle->new(" | samtools view -h -b > $file") or die "Could not open $file";
  {key=>$key, file=>$file, fh=>$fh};
}

sub close_bamout { 
  my($h)=@_;
  $h->{fh}->close() || die "Could not write (or close) $h->{file}";
  undef;
}
