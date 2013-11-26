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



package PSF;
# taken from cdfs-extract; should always be kept up to date..

=head1 NAME

PSF -- manipulation of XMM PSFs

=head1 SYNOPSIS

my $psf = PSF->new($offaxis,$energy[,$camera][,$oldpsf]);
 # $offaxis in arcmin, $energy in keV

my $pdl2d = $psf->gettable;

# the following are used by cdfs-extract, see that code for documentation
$psf->rotatepsf($angle);
$psf->place_in_observation($expmap); # $expmap is an object
$psf->eef($expmap);

=head1 DESCRIPTION

* PSF->new($offaxis,$energy[,$camera][,$oldpsf]);
 reads the correct PSF from the latest CCF. If $oldpsf==1 or $camera==undef,
 then use old XRT1_XPSF007.CCF where PSF were the same for all cameras

=cut

use Carp;
use PDL;

sub new {
    my ($class, $offaxis,$energy,$camera,$oldpsf) = @_;
    my $self = {};

    $self->{OFFAXIS} = $offaxis;
    $self->{ENERGY} = $energy;

    # check for values
    if ($offaxis<0 or $energy<0 or $energy>15) {
	die "Something wrong with offaxis (=$offaxis) or energy (=$energy).\n";
    }
    if ($offaxis>16) {
	print("Large value of offaxis ($offaxis).\nContinuing, but please check the alignment of event and expmap files.\n");
	$offaxis=16;
    }

    # locate PSF library
    my $CCF;
    #printf("camera=%s\n",$camera);
    if ($oldpsf or not defined($camera)) { 
	$CCF = $ENV{SAS_CCFPATH}.'/XRT1_XPSF_0007.CCF';
    } else {
	if ($camera eq 'EMOS1') {
	    $CCF = $ENV{SAS_CCFPATH}.'/XRT1_XPSF_0013.CCF';
	} elsif ($camera eq 'EMOS2') {
	    $CCF = $ENV{SAS_CCFPATH}.'/XRT2_XPSF_0013.CCF';
	} elsif ($camera eq 'EPN') {
	    $CCF = $ENV{SAS_CCFPATH}.'/XRT3_XPSF_0013.CCF';
	} else {
	    croak("Cannot understand which instrument $camera is.\n");
	}
    }

    # which PSF table do we need?
    # in the CCF, the PSF images start from extension #4, and cycle like this:
    #  100 eV - 00 arcmin    ext  #4
    #  100 eV - 03 arcmin    ext  #5
    #  100 eV - 06 arcmin    ext  #6
    #  100 eV - 09 arcmin    ext  #7
    #  100 eV - 12 arcmin    ext  #8
    #  100 eV - 15 arcmin    ext  #9
    # 1500 eV - 00 arcmin    ext #10
    # 1500 eV - 03 arcmin    ext #11
    # 1500 eV - ...
    # 3000 eV - ...     from ext #16
    # 4500 eV - ...     from ext #22
    # 6000 eV - ...     from ext #28
    # 7500 eV - ...     from ext #34
    # 9000 eV - ...     from ext #40
    #10500 eV - ...     from ext #46
    #12000 eV - ...     from ext #52
    #13500 eV - ...     from ext #58
    #15000 eV - ...     from ext #64
    my $extension = 4 + 6*rint($energy/1.5) + rint($offaxis/3);

    $self->{PSFIMAGE} = rfits($CCF."[$extension]");
    # - the pixel size of the PSF table is 1"x1"
    # - the table is 512x512 pixels
    # - the table values are normalised so that the sum is 625

    bless( $self, $class );
    return($self);
}


sub gettable {
    my $self = shift;

    return $self->{PSFIMAGE};
}


sub rotatepsf {
    use PDL::Transform;

    my ($self,$angle) = @_;

    my $t = new PDL::Transform::Linear
	({rot=>$angle, pre=>[-256,-256],post=>[+256,+256]});

    # while method=>linear is default, here we specify it to avoid a
    # warning (Use of uninitialized value in pattern match (m//) at
    # /Library/Perl/5.8.6/darwin-thread-multi-2level/PDL/Transform.pm
    # line 902) which appears when 'use strict' is on;
    $self->{PSFIMAGE} = $t->map( $self->{PSFIMAGE},
				 {pix=>1,method=>'linear'} );
    return;
}



sub frac_in_radius {
    my ($self,$radius) = @_;	# radius in arcsec

    my $mask = rvals($self->{PSFIMAGE},{squared=>1}) < ($radius*$radius);

    return( $self->{PSFIMAGE}->where($mask)->sum / 625 );
}


sub place_in_observation {
    use PDL::Transform;

    my ($self,$expmap) = @_;

    my $scale = 1/abs($expmap->{IMG}->hdr->{CDELT1} * 3600);
    # (psfs are already binned at 1 arcsec/pixel)

    my ($x,$y,$r) = $expmap->first_src_area;

    my $transform = new PDL::Transform::Linear
	({pre=>[-256,-256], scale=>$scale, post=>[$x,$y]});

    my @naxis = ($expmap->{IMG}->hdr->{NAXIS1},$expmap->{IMG}->hdr->{NAXIS2});
    my $new = $transform->map($self->{PSFIMAGE},\@naxis,{pix=>1,method=>'linear'});
    $new  /= $new->sum;  # renorm to 1
    return($new);

}


sub eef { # encircled energy fraction
    my ($self,$expmap) = @_;

    my $placed_psf = $self->place_in_observation($expmap);
    wfits $placed_psf,'a.fits';
    wfits $expmap->{IMG},'b.fits';

    (my $mask, undef) = $expmap->exposed_masks;
    wfits $mask,'c.fits';

    return( $placed_psf->where($mask)->sum );

}




1;
