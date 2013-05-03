#!/usr/bin/perl
#


# the XMM-CDFS simulator of astronomical X-ray observations
# Copyright (C) 2013 Piero Ranalli
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
# 
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.



=head1 NAME

obsid-sim - simulations for xmm-cdfs

=head1 SYNOPSIS


usage:  SAS_CCF=/path/to/ccf.cif ./obsid-sim.pl \
        counts-positions.dat spectrum.qdp real-event-file.fits \
        expmap.fits sim.fits [--oldpsf] [--pimin=2] [--pimax=10]

=head1 DESCRIPTION

The file counts-positions.dat has the following format
 COUNTS    RA    DEC

NB if the expmap is specified, then RATES are expected instead of counts;

spectrum.qdp is a model obtained in XSPEC with the following
commands:

 XSPEC> iplot model
 IPLOT> wd spectrum

and then editing the file to remove the commands (READ SERR etc.)
at its beginning;

and real-event-file.fits is a real-data event file, which is used
to extract the information about the boresight, by both this program
and the SAS task esky2det.

The expmap is optional (use 'none' instead of filename).

Output: a file named sim.fits .

As of SAS 8.0.0 there were no differences in the PSF models for
MOS1/MOS2/PN); this changed in SAS 9. The old behaviour can be
reproduced with the --oldpsf flag.

The effective area is different for the
different cameras, and the relevant camera is taken from the
real-event-file.


=head1 VERSION

Version 2: now multiplies by expmap after finishing the simulation, 
to reproduce also the CCD structure: gaps, dead columns, etc.

Version 2.1:
un'altra cosa sistemata sono le sorgenti che contribuiscono un numero
frazionario di eventi (molto deboli, servono a fare il fondo cosmico); la
soluzione e' arrotondare $counts prima del ciclo principale, conservandone
la somma.

Version 2.2:
fast resampling of PSFs with cumulants

Version 2.3:
added ftools pfiles tmp directories;
splitted mos/pn columns in counts-pos

Version 2.4 - 2011/12/1:
uses new PSFs by default  (retains old behaviour with --oldpsf flag);
the energy_cuts can be chosen from command line (or default to 2-10 keV)

Version 2.5 - 2011/12/2:
several bugs corrected - release used to produce the PSF for P34H_352

Version 2.51 - 2011/12/13:
old subroutines no longer used (samplePSF, metropolis) moved out of scope

Version 2.6 - 2012/4/6
added header fixes, the outfiles can now be read by SAS

Version 2.61 - 2012/5/16
move WCS dependance to WCS2
runs also on hobbes

=cut


use Astro::FITS::Header;
use PDL;
use PDL::NiceSlice;
use Piero::WCS2;
use Piero::Ftools::Pfiles;
use Getopt::Long;
use FixFitsHdr;

use strict;
use warnings;


$PDL::BIGPDL=1;

my ($oldpsf,$pimin,$pimax);
GetOptions('oldpsf' => \$oldpsf,
	   'pimin=f' => \$pimin,
	   'pimax=f' => \$pimax );

my $counts_posfile = $ARGV[0];
my $spectrumfile = $ARGV[1];
my $realeventfile= $ARGV[2];
my $expmapfile = $ARGV[3] eq 'none' ? 0 : $ARGV[3];
my $outfile = $ARGV[4];

my @energy_cuts =
    ($pimin and $pimax) ?
    ($pimin,$pimax) :
    (2,10); # keV

printf("Ref: %s\n",$realeventfile);
printf("Energy cuts: %s-%s keV\n",@energy_cuts);

our $pfiles = Pfiles->new; # a global object

# get camera, boresight and position angle
my $refhdr = rfits($realeventfile,{data=>0});
my ($RA_PNT,$DEC_PNT,$PA_PNT) = getboresight($refhdr);
my $camera = getcamera($refhdr);

my ($counts_mos,$counts_pn,$ra,$dec) = rcols($counts_posfile);
my $counts = $camera eq 'EPN' ?
    $counts_pn :
    $counts_mos;


# reads model spectrum, cut, multiply by effective area and return
# in the form of a (pdl) event list
my $evtlist = spec2evtlist($spectrumfile,$camera,@energy_cuts);

# read expmap
my ($expmap, $exposure);
if ($expmapfile) {
    $expmap = rfits($expmapfile);
    # normalise expmap
    my $exposure = $expmap->flat->maximum->sclr;
    $expmap /= $exposure;
    # and convert rates into counts
    $counts *= $exposure;
}


# calculate offaxis
my $offaxis = sqrt(($ra-$RA_PNT)**2 + ($dec-$DEC_PNT)**2); # in degrees
$offaxis *= 60;  # in arcmin

# eliminate out-of-FOV sources
my $idx = which($offaxis<=17); # the threshold here is a bit larger
                               # than the one is samplePSF, to allow
                               # the check of alignment
$offaxis = $offaxis->index($idx);
$ra = $ra->index($idx);
$dec = $dec->index($idx);
$counts = $counts->index($idx);

# order sources according to offaxis angle
$idx = qsorti($offaxis);
$offaxis = $offaxis->index($idx);
$ra = $ra->index($idx);
$dec = $dec->index($idx);
$counts = $counts->index($idx);

# counts may be fractionary: round them, conserving their sum
# (otherwise all the faint sources get lost: with 0.1 counts each,
# they are undetectable yet make the cosmic background)
$counts=round_with_sum($counts);

my $evtheader = rfits($realeventfile.'[1]',{data=>0});
my ($x,$y) = wcs_evt_transfinv($evtheader,$ra,$dec);
my ($x_PNT,$y_PNT) = wcs_evt_transfinv($evtheader,$RA_PNT,$DEC_PNT);
#print(join(' ',$x_PNT,$y_PNT,$evtheader->{REFXCRPX},$evtheader->{REFYCRPX},$evtheader->{REFXCRVL},$header->{RA_PNT}"\n"));
#print(join(' ',$evtheader->{TCRVL6},$evtheader->{RA_PNT},$evtheader->{TCRVL7},$evtheader->{DEC_PNT},"\n"));


my $angle = atan2($y-$y_PNT,$x-$x_PNT)*180/3.141592 +90; # - $PA_PNT +90;



my $PSF_Ebin_centre = 1.5*sequence(11);


my $src_ra  = pdl([]);
my $src_dec = pdl([]);
my $src_E  = pdl([]);


# main cycle, on sources
for my $i (0..($ra->dim(0) -1)) {

    next unless ($counts->($i)->sclr);

    my $thissrc_ra = my $thissrc_dec = pdl([]);
    # bootstrap on energies (sorting, so that the energies will
    # reflect the PSF)
    my $thissrc_E = qsort(bootstrap($evtlist,$counts->($i)->sclr));
    # bin on the available PSF energies
    my $src_Ehisto = histogram($thissrc_E,1.5,-.75,11);

    # cycle on energy bins
    for my $j (0..10) {

	# do nothing if there are no counts in the bin
	next unless ($src_Ehisto->($j)->sclr);

	my ($d_ra,$d_dec) = samplePSF_fast($offaxis->($i)->sclr,
				      $PSF_Ebin_centre->($j)->sclr,
				      $angle->($i)->sclr, 
				      $src_Ehisto->($j)->sclr,
				      $camera,
	                              $oldpsf );

	$thissrc_ra  =  $thissrc_ra->append($d_ra + $ra->($i)->sclr);
	$thissrc_dec = $thissrc_dec->append($d_dec+$dec->($i)->sclr);

    }


    $src_ra  =  $src_ra->append($thissrc_ra);
    $src_dec = $src_dec->append($thissrc_dec);
    $src_E   =   $src_E->append($thissrc_E);


}


# apply expmap
$idx = screen_by_expmap($src_ra,$src_dec,$expmap); # reuse $idx
if ($idx->dim(0)<$src_ra->dim(0)) {  # BUG corrected - was $ra
    #print "yeah\n";
    $src_ra=$src_ra->index($idx);
    $src_dec=$src_dec->index($idx);
    $src_E=$src_E->index($idx);
}







# get the X,Y columns
my ($src_x,$src_y) = wcs_evt_transfinv($evtheader,$src_ra,$src_dec);


# get the TIME column
# take times as the first n timestamps in the realeventfile
my $src_time = sampletime($realeventfile,$src_x->dim(0));


# prepare the fits file for output (still missing RAW and DET coordinates):
my $sim = {};
$sim->{tbl} = 'binary';
$sim->{hdr} = $evtheader;

adjustheader($sim);

# add columns
$sim->{RA} = $src_ra;
$sim->{DEC} = $src_dec;
$sim->{X} = $src_x;
$sim->{Y} = $src_y;
$sim->{TIME} = $src_time;
$sim->{PHA} = 1000*$src_E;
$sim->{PI} = 1000*$src_E;
$sim->{PATTERN} = zeroes($src_E);
$sim->{FLAG} = zeroes($src_E);

# output
wfits($sim,$outfile);



# fix the header
fixfitshdr($outfile, $evtheader, $realeventfile, $pfiles, {COPY_XY=>1});



# end



# columns that should be present in the xxsim.fits:
# TIME   can be random or fixed, pwxdetect should not care
# RAWX/Y  depends on RA,DEC + OBSID + CAMERA
# DETX/Y  depends on RA,DEC + OBSID + CAMERA
# X/Y     depends on RA,DEC + OBSID + CAMERA
# PHA,PI   can be fixed - or should reflect spectrum if one wants to
#          filter the file after?
# FLAG = 0
# PATTERN = 0
# CCDNR    depends on RA,DEC + OBSID + CAMERA



sub spec2evtlist {
    my $spectrumfile = shift;
    my $camera = shift;
    my $Emin = shift;
    my $Emax = shift;

    # reads model spectrum
    my ($spec_E,$spec_dE,$spec_F) = rcols($spectrumfile);

    # cut
    $spec_F *= ($spec_E >= $Emin);
    $spec_F *= ($spec_E <= $Emax);

    # read the effective area from the CCF
    my $effarea;
    if ($camera =~ m/EMOS1/) {
	$effarea = rfits($ENV{SAS_CCFPATH}.'/XRT1_XAREAEF_0008.CCF[1]');
    } elsif ($camera =~ m/EMOS2/) {
	$effarea = rfits($ENV{SAS_CCFPATH}.'/XRT2_XAREAEF_0009.CCF[1]');
    } elsif ($camera =~ m/EPN/) {
	$effarea = rfits($ENV{SAS_CCFPATH}.'/XRT3_XAREAEF_0011.CCF[1]');
    } else {
	die "Camera not recognised in $spectrumfile.\n";
    }
    # interpolate effective area on the spectrum
    my ($spec_Area,$err) = interpolate($spec_E,$effarea->{ENERGY}/1000,
				       $effarea->{AREA});

    # multiply spectrum by eff area
    $spec_F *= $spec_Area;


    # normalise: we want a "number of photons per bin", totalling 1000 photons on
    # the whole spectrum
    my $norm = $spec_F->sumover;
    $spec_F = rint( $spec_F * 1000 / $norm );
    # now convert to list of photons:
    my @lspec_E  = $spec_E->list;
    my @lspec_dE = $spec_dE->list;
    my @lspec_F  = $spec_F->list;

    my $evtlist = pdl([]);
    for my $i (0..$#lspec_E) {

	next unless $lspec_F[$i];  # nothing to do with 0 photons in the bin

	my $energies = $lspec_E[$i] + $lspec_dE[$i]*random($lspec_F[$i]);
	$evtlist = $evtlist->append($energies);
    }
    # the list of photons is ready to be bootstrapped later
    return($evtlist);
}


sub getboresight {
    my $hdr = shift;
    return($hdr->{RA_PNT},$hdr->{DEC_PNT},$hdr->{PA_PNT});

}

sub getcamera {
    my $hdr = shift;
    return($hdr->{INSTRUME});

}



sub samplePSF_fast {
# 	my ($d_ra,$d_dec) = samplePSF($offaxis->(($i)),
# 				      $PSF_Ebin_centre->(($j)),
# 				      $PA_PNT, 
# 				      $counts->(($i)) );
    use PDL::NiceSlice;
    use lib '..';
    use PSF;

    my ($offaxis,$energy,$PA,$counts,$camera,$oldpsf) = @_;

    my $psf = PSF->new($offaxis,$energy,$camera,$oldpsf);
    my $PSFtable = $psf->gettable;
    # - the pixel size of the PSF table is 1"x1"
    # - the table is 512x512 pixels
    # - the table values are normalised so that the sum is 625

    # sample using the cumulants method: 
    # (find_value should thread on both the indexes and the table dimension)
    # first sample the y
    my $j = 512*random($counts);
    my $cumu = $PSFtable->sumover->cumusumover;
    $cumu *= 512/625; # normalize to 512
    my $y = vsearch($j,$cumu);
    # now the x
    $j = 512*random($counts);
    $cumu = $PSFtable->dice('X',$y)->cumusumover;
    $cumu *= 512/$cumu(-1,:); # normalize each row to 512
    my $x = vsearch($j,$cumu);

    # convert to ra & dec
    $x -= 256;
    $y -= 256;

    $PA *= 3.141592 / 180; # deg->rad
    my $dra =  $x*cos($PA) + $y*sin($PA);
    my $ddec= -$x*sin($PA) + $y*cos($PA);

    $dra  /= 3600; # arcsec->deg
    $ddec /= 3600; # arcsec->deg

    return($dra,$ddec);

}




sub bootstrap {
    my $piddle = shift;
    my $n = shift;

    my $idx = rint( ($piddle->dim(0) -1) * random($n) );

    return($piddle->index($idx));
}




sub sampletime {
    my ($evt,$counts) = @_;
    my $pf = $main::pfiles->env;

    # how many rows in the file?
    my $n2 = `$pf fkeyprint $evt+1 NAXIS2 | awk '/NAXIS2  =/ { print \$3}'`;   
    my $samples = $n2 < $counts ? $n2 : $counts;

    open(TIME, "$pf fdump $evt+1 STDOUT TIME 1-".$samples.' prhead=no showcol=no showunit=no showrow=no page=no  | ');
    my $src_time = rcols(*TIME,{IGNORE=>'/^\s*$/'});
    close(TIME);
#     if ($src_time->dim(0) != $counts->(($i))) {
# 	die "Not enough photons in $realeventfile,\n from which I'd like to take timestamps.\n";
#     }

    # need to resample?
    if ($n2 < $counts) {
	$src_time = bootstrap($src_time,$counts);
    }

    return($src_time);
}

sub screen_by_expmap {
# my $idx = screen_by_expmap($thissrc_ra,$thissrc_dec,$expmap);

    my ($ra,$dec,$expmap) = @_;

    # get image coordinates for ra,dec pairs
    my ($img_x,$img_y) = wcstransfinv($expmap->hdr,$ra,$dec);
    $img_x = rint($img_x);
    $img_y = rint($img_y);

    # check (and correct) for absurd coordinates
    my $msk = $img_x >= 0;
    $msk *= $img_x <= $expmap->hdr->{NAXIS1};
    $msk *= $img_y >= 0;
    $msk *= $img_y <= $expmap->hdr->{NAXIS2};

    if (any(! $msk)) {
	$img_x = $img_x->where($msk);
	$img_y = $img_y->where($msk);
	print "absurd coordinates corrected\n";
    }

    #wfits($img_x,"deb-imgx.fits");
    #wfits($img_y,"deb-imgy.fits");

    # get probability values from expmap
    my $prob = $expmap->index2d($img_x,$img_y);
    #wfits($prob,"deb-prob.fits");


    # and return the indices of surviving photons
    my $rand = random($prob);
    my $idx = which($rand<$prob);

    #wfits($idx,"deb-idx.fits");
    return($idx);
}

sub mul_by_expmap {
# my $counts = mul_by_expmap($rates,$ra,$dec,$expmap);

    my ($rates,$ra,$dec,$expmap) = @_;

    # get image coordinates for ra,dec pairs
    my ($img_x,$img_y) = wcstransfinv($expmap->hdr,$ra,$dec);
    $img_x = rint($img_x);
    $img_y = rint($img_y);

    # check (and correct) for absurd coordinates
    my $msk = $img_x >= 0;
    $msk *= $img_x <= $expmap->hdr->{NAXIS1};
    $msk *= $img_y >= 0;
    $msk *= $img_y <= $expmap->hdr->{NAXIS2};

    if (any(! $msk)) {
	$img_x = $img_x->where($msk);
	$img_y = $img_y->where($msk);
	$rates = $rates->where($msk);
	print "Out-of-expmap sources eliminated.\n";
    }
    my $idx = which($msk);

    # get values from expmap
    my $exposures = $expmap->index2d($img_x,$img_y);

    # and return the number of counts
    return($exposures * $rates,$idx);
}




sub round_with_sum {
    # round a piddle conserving (in a statistical sense) the sum
    # the conservation of the sum is accurate

    # fotoni interi:  floor(cts)
    # fotoni frazionari: cf = cts-floor(cts)
    # somma delle frazioni: \sum cf
    # probabilita' di emettere un fotone:  cf_i

    my $x = shift;
    my $x0 = floor($x);   # integer photons
    my $cf = $x - $x0;    # fractionary counts
    my $r = random($x);   # probability of receiving a frac.count
    my $y = $r <= $cf;    # quantization
    return( $x0+$y );

}





__END__

