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



package FixFitsHdr;

# fix the headers of the fits files produced by the simulator,
# so that they can be read by the SAS

# this file should always come with /Volumes/data/SIM-CDFS/Data/cdfs-sim-fixhdr1.dat

# PR 2012/4/6

use Exporter;
@ISA = qw/Exporter/;

@EXPORT = qw/adjustheader fixfitshdr/;

use PDL;
use File::Temp 'mktemp';
use 5.010;
use strict;
use warnings;
use Carp;

our $EVTCHUNKMAXSIZE = 20000;
our $pf; # PFILES environment

sub adjustheader {
    my $fits = shift;

    # to save WCS keys for later use:
    # find which columns are X and Y
    my ($xcol,$ycol,$racol,$deccol);
    while (my ($key,$value) = each %{$fits->{hdr}}) {

	if ($key =~ m/TTYPE/ and $value eq 'X') {
	    ($xcol) = ($key =~ m/TTYPE(\d)/);
	}
	if ($key =~ m/TTYPE/ and $value eq 'Y') {
	    ($ycol) = ($key =~ m/TTYPE(\d)/);
	}
	if ($key =~ m/TTYPE/ and $value eq 'RA') {
	    ($racol) = ($key =~ m/TTYPE(\d)/);
	}
	if ($key =~ m/TTYPE/ and $value eq 'DEC') {
	    ($deccol) = ($key =~ m/TTYPE(\d)/);
	}
    }
    # and save WCS info relative to those columns
    my @keystosave = qw/TNULL TLMIN TDMIN TLMAX TDMAX TCRPX TCTYP TCRVL TCDLT TCUNI/;
    my %savedvalues;
    my @list = ($xcol,$ycol);
    push(@list,$racol) if defined($racol);
    push(@list,$deccol) if defined($deccol);
    for my $i (@list) {
	for my $key (@keystosave) {
	    $savedvalues{"$key$i"} = $fits->{hdr}->{"$key$i"};
	}
    }


    # delete all keys about columns
    for my $key (keys %{$fits->{hdr}}) {
	if ($key =~ m/T[LD]M(IN|AX)/ or
	    $key =~ m/TC?TYP/ or
	    $key =~ m/TFORM/ or
	    $key =~ m/TNULL/ or
	    $key =~ m/TCRPX/ or
	    $key =~ m/TCRVL/ or
	    $key =~ m/TCDLT/ or
	    $key =~ m/TCUNI/ or
	    $key =~ m/^DS/ or
	    $key =~ m/TUNIT/ ) {

	    delete($fits->{hdr}->{$key});
	}
    }

    # and now re-add them...
    my @type = qw/TIME X Y PHA PI PATTERN FLAG RA DEC/;
    my @form = qw/D J J I I B B D D/;
    my @unit = qw/s pixel pixel CHAN CHAN _ _ deg deg/;
#     my @type = qw/TIME X Y RA DEC PHA PI PATTERN/;
#     my @form = qw/D J J D D I I B/;
#     my @unit = qw/s pixel pixel deg deg CHAN CHAN _/;

    for my $i (1..$#type+1) {
	$fits->{hdr}->{"TTYPE$i"} = $type[$i-1];
	$fits->{hdr}->{"TFORM$i"} = $form[$i-1];
	$fits->{hdr}->{"TUNIT$i"} = $unit[$i-1];
    }

    # readd also WCS info
    # in the lists above, X and Y are 2 and 3 respectively
    for my $key (@keystosave) {
	$fits->{hdr}->{"${key}2"} = $savedvalues{"$key$xcol"};
	$fits->{hdr}->{"${key}3"} = $savedvalues{"$key$ycol"};
    }

}


sub fixfitshdr {
    my $outfile = shift;  # the file to be modified
    my $evtheader = shift; # the template header
    my $infile = shift;   # the template file
    my $pfiles = shift;   # a Piero::Ftools::Pfiles object
    my $opt = shift;

    # proton simulations have the problem that the WCS info for X,Y
    # is lost when the data type is changed from D to I (see below)
    # and the WCS is not recovered as in other cases
    # So we save it and re-add it later
    my $savehdr = rfitshdr($outfile.'[1]') if ($$opt{COPY_XY});

    unless (defined($pfiles)) {
	require Piero::Ftools::Pfiles;
	$pfiles = Piero::Ftools::Pfiles->new;
    }

    my $fixhdrdat = $$opt{FIXHDRDAT} // 'cdfs-sim-fixhdr1.dat';

    $pf = $pfiles->env;  # global

    # delete the continue keywords, set creator and exp_id
    system("$pf fthedit $outfile+0 @$fixhdrdat");
    system("$pf fthedit $outfile+1 @$fixhdrdat");

    # aggiungere all'estensione +0 anche le keyword INSTRUME,RA_PNT,DEC_PNT,PA_PNT,DATE-OBS
    my $tmp1 = mktemp('ftoolstmp-XXXXXX');
    open(my $fix,'>', $tmp1);
    for my $key (qw/INSTRUME RA_PNT DEC_PNT PA_PNT DATE-OBS DATAMODE/) {
	print $fix "$key = ".$evtheader->{$key}."\n";
    }
    close($fix);
    system("$pf fthedit $outfile+0 ".'@'.$tmp1);
    unlink($tmp1);

    # now we have to add the RAWX/Y, DETX/Y, CCDNR and FLAG
    # this is a task for esky2det, which can add the needed columns to
    # the file

    #### following code no longer needed after esky2det bug corrected in SAS 13
    #### (and anyway, this workaround was not working)
    # # check if the event file has NAXIS2>32000; if so, split and run esky2det
    # # on the single chunk
    # my $hdr = rfitshdr($outfile.'[1]');
    # if (0 and $hdr->{NAXIS2} > $EVTCHUNKMAXSIZE) {
    # 	my @pieces = splitfits($outfile,$hdr->{NAXIS2});
    # 	esky2det($_) for @pieces;
    # 	rejoin(@pieces,$outfile);
    # } else {
    # 	esky2det($outfile);
    # }

    esky2det($outfile);


    my $tmp2 = mktemp('ftoolstmp-XXXXXX');
    my $tmp3 = mktemp('ftoolstmp-XXXXXX');
    my $tmp4 = mktemp('ftoolstmp-XXXXXX');
    #my ($tmp1,$tmp2,$tmp3,$tmp4) = qw/tmp1 tmp2 tmp3 tmp4/;


    # the DETX,DETY column have type E. Convert to I:
    system("$pf ftcalc $outfile+1 $tmp1 DX '(int) DETX' tform='1I'");
    system("$pf ftcalc $tmp1+1 $tmp2 DY '(int) DETY' tform='1I'");
    system("$pf fdelcol $tmp2+1 DETX N Y");
    system("$pf fdelcol $tmp2+1 DETY N Y");
    system("$pf ftcalc $tmp2+1 $tmp3 DETX DX tform='1I'");
    system("$pf ftcalc $tmp3+1 $tmp4 DETY DY tform='1I'");
    system("$pf fdelcol $tmp4+1 DX N Y");
    system("$pf fdelcol $tmp4+1 DY N Y");
    system("mv -f $tmp4 $outfile");
    unlink($tmp1,$tmp2,$tmp3,$tmp4);


    # X,Y also are E. Convert to J:
    system("$pf ftcalc $outfile+1 $tmp1 DX '(int) X' tform='1J'");
    system("$pf ftcalc $tmp1+1 $tmp2 DY '(int) Y' tform='1J'");
    system("$pf fdelcol $tmp2+1 X N Y");
    system("$pf fdelcol $tmp2+1 Y N Y");
    system("$pf ftcalc $tmp2+1 $tmp3 X DX tform='1J'");
    system("$pf ftcalc $tmp3+1 $tmp4 Y DY tform='1J'");
    system("$pf fdelcol $tmp4+1 DX N Y");
    system("$pf fdelcol $tmp4+1 DY N Y");
    system("mv -f $tmp4 $outfile");
#    system("cp $tmp4 $outfile");
    unlink($tmp1,$tmp2,$tmp3,$tmp4);



    # copy WCS info for DETX/Y from template to simulation
    # re-read header
    copywcsinfo(qr/^DETX\s*$/,$infile,$outfile);
    copywcsinfo(qr/^DETY\s*$/,$infile,$outfile);
    if ($$opt{COPY_XY}) {
	copywcsinfo(qr/^X\s*$/,$savehdr,$outfile);
	copywcsinfo(qr/^Y\s*$/,$savehdr,$outfile);
    }
}


# sub no longer needed after esky2det bug correction
# left in case may be needed in the future
#
# sub splitfits {
#     use PDL::NiceSlice;

#     my $fits = shift;
#     my $size = shift;

#     my $nchunks = sclr(ceil( $size / $EVTCHUNKMAXSIZE ));

#     my $cmd = sprintf("%s ftstat '%s[events][col TIME]'",$pf,$fits);
#     my @stats = `$cmd`;

#     my $min = extract(qr/min/,@stats);
#     my $max = extract(qr/max/,@stats);

#     my $seq = $min + ($max-$min)/$nchunks * sequence($nchunks+1);
#     my $mins = $seq(0:-2);
#     my $maxs = $seq(1:-1);

#     my @pieces;
#     for my $c (1..$nchunks) {
# 	my $op = $c == $nchunks ? '<=' : '<';
# 	my $name = $fits;
# 	$name =~ s/.fits/_chunk$c.fits/;
# 	push(@pieces,$name);

# 	$cmd = sprintf("%s evselect table=%s:EVENTS withfilteredset=yes filteredset=%s expression='(TIME>=%s) && (TIME%s%s)'",
# 		      $pf,$fits,$name,sclr($mins($c-1)),$op,sclr($maxs($c-1)));
# 	system($cmd);
#     }

#     return @pieces;
#     no PDL::NiceSlice;
# }

# sub extract {
#     my $qr = shift;
#     my @lines = @_;

#     my @ok = grep(/$qr/,@lines);
#     my @cols = split(' ',$ok[0]);
#     return $cols[1];
# }


# sub no longer needed after esky2det bug correction
# left in case may be needed in the future
# sub rejoin {
#     my $outfile = pop;
#     my @pieces = @_;

#     my $tmp1 = mktemp('ftoolstmp-XXXXXX');
#     my $tmp2 = mktemp('ftoolstmp-XXXXXX');
#     my $first = shift @pieces;
#     system("mv $first $tmp1");

#     while (@pieces >= 1) {
# 	my $next = shift @pieces;
# 	mergevt($tmp1,$next,$tmp2);
#     }

#     system("mv $tmp1 $outfile");
#     #unlink($tmp1,$tmp2);
# }

# sub mergevt {
#     my ($one,$two,$tmp) = @_;

#     my $cmd = sprintf("merge set1=%s set2=%s outset=%s",
# 		      $one,$two,$tmp);
#     say $cmd;
#     system($cmd);
#     if ($?) {
# 	croak "Died on $cmd";
#     }
#     system("mv $tmp $one");
# }


sub esky2det {
    my $outfile = shift;

    # now we have to add the RAWX/Y, DETX/Y, CCDNR and FLAG
    # this is a task for esky2det, which can add the needed columns to
    # the file
    # first call, adding the RAWX,RAWY, CCDNR and FLAG
    # print(`esky2det datastyle=set intab=sim.fits:EVENTS witherrorcol=no withouttab=no outunit=raw calinfostyle=set calinfoset=sim.fits `);
    # # second call, adding DETX,DETY
    my $command = "$pf esky2det datastyle=set intab=$outfile:EVENTS witherrorcol=no withouttab=no outunit=det calinfostyle=set calinfoset=$outfile";
    say $command;
    print(`$command`);
    if ($?) {
	sleep(rand());
	say "Retry $command";
	print(`$command`);
	if ($?) {
	    die 'Giving up with esky2det';
	}
    }
}


sub copywcsinfo {
    my $key = shift;
    my $infile = shift;
    my $outfile = shift;

    my $in = ref($infile) eq 'HASH' ? $infile : rfitshdr($infile);
    my $out= rfitshdr($outfile);

    my $in_dx  = findTTYPE($in, $key, 'in');
    my $out_dx = findTTYPE($out,$key, 'out');


    for my $k (qw/TTYPE TUNIT TLMIN TDMIN TLMAX TDMAX TCRPX TCTYP TCRVL TCDLT TCUNI/) {

	if (exists($in->{${k}.${in_dx}})) {
	    # remove the spaces which are sometimes found in TUNIT (pn files)
	    if ($k =~ m/TUNIT/) {
		$in->{${k}.${in_dx}}  =~ s/\s/_/g;
	    }

	    system("$pf fthedit $outfile+1 ${k}${out_dx} add ".$in->{${k}.${in_dx}});
	    #	system("$pf fthedit $outfile+1 ${k}${out_dy} add $in->{$k.$in_dy}");
	}

    }

}


sub findTTYPE {
    my $hdr = shift;
    my $key = shift;
    my $where = shift;

    my @ttype = grep(/^TTYPE/, keys %$hdr);

    for my $tt (@ttype) {

	if ($hdr->{$tt} =~ m/$key/) {
	    # get the number
	    $tt =~ m/^TTYPE(\d+)/;
	    return $1;
	}else{
	    print"$tt $hdr->{$tt}\n";
	}
    }
    die "FixFitsHdr error: could not find $key among the TTYPEs in $where\n";
}





1;
