#!/usr/bin/env perl
#

=head1 NAME

calconaxis.pl -- call ecoordconv to get on-axis position and store it

=head1 SYNOPSIS

./calconaxis.pl eventfile.fits ccf.cif sum.sas  # stores ra,dec of on-axis position in onaxis.yaml
./obsid-sim.pl ...              # reads from onaxis.yaml

=head1 DESCRIPTION

See this helpdesk ticket: http://xmm2.esac.esa.int/xmmhelp/Calibration?id=16476;page=8;user=guest
where M. Guainazzi explains:

<<The expected location (DETX,DETY) of the proposal coordinates in EPIC
cameras for pn prime is

 MOS1: (310,-221) (169,-140 if RGS1 prime)
 MOS2: ( 38,  87) (-42, -55 if RGS1 prime)
 pn  : ( 91,-228) (172, -88 if RGS1 prime)

Major sources of scattering for these values when measured in actual SAS
processed data may be

(a) pointing error (typically 2-3 arcsec)

(b) accuracy of target coordinates (1-2arcsec)>>


Here we only consider EPIC as the prime instrument.


=cut


use Fcntl qw/:DEFAULT :flock/;
use PDL;

use FindBin;
use lib "$FindBin::Bin/lib";

use YAML::Tiny;
use XMMSAS::Tiny;

use strict;
use warnings;

die "Need to provide all three parameters: event file, ccf.cif, ...sum.sas.\n";

my ($reffile,$ccfcif,$sumsas) = @$ARGV;


die "Cannot find file $reffile" unless (-e $reffile);

my ($ra,$dec) = calconaxis( $reffile,$ccfcif,$sumsas );


my $conffile = 'onaxis.yaml';

sysopen(my $fh, $conffile, O_RDWR) or die "Cannot open $conffile\n";
flock( $fh, LOCK_EX ) or die "Cannot lock $conffile\n";

# does the file already exists?
my $yaml;
unless ($yaml = YAML::Tiny->read( $conffile )) {
    $yaml = YAML::Tiny->new;
}

$yaml->[0]->{$reffile} = {
			  RA => $ra,
			  DEC => $dec,
			 };

$yaml->write( $conffile );

# lock is released by exit

exit;


sub calconaxis {
    my ($reffile,$ccfcif,$sumsas) = @_;
    my $sas = XMMSAS::Tiny->new;

    my $hdr = rfitshdr($reffile);

    my ($x,$y);
    if ($hdr->{INSTRUME} eq 'EMOS1') {
	$x = 310;
	$y = -221;
    } elsif ($hdr->{INSTRUME} eq 'EMOS2') {
	$x = 38;
	$y = 87;
    } elsif ($hdr->{INSTRUME} eq 'EPN') {
	$x = 91;
	$y = -228;
    }

    $sas->ccf($ccfcif);
    $sas->odf($sumsas);
    my $cmd = "ecoordconv imageset=$reffile x=$x y=$y coordtype=det";

    my @out = $sas->call($cmd);

    my @radec = grep('RA: DEC:', @out);
    my $r = $radec[0];
    chomp($r);
    my @elems = split(' ',$r);
    return($elems[2],$elems[3]);
}
