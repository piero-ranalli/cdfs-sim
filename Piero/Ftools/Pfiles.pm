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




package Pfiles;

=pod

NAME

package Pfiles

SYNOPSIS

my $pfiles = Pfiles->new;
system($pfiles->env." ftools_call");

DESCRIPTION

A temporary subdirectory named ftools-tmppfiles-XXXX is created in
the current dir. All files from $ENV{PFILES} are copied there.
env() returns the string "PFILES=tmpdir".
The object destructor removes the tmpdir at the end.

=cut

use File::Temp qw/tempdir/;
use Carp;

sub new {
    my $class = shift;
    my $self = {};

    my $tmpdir = tempdir( "ftools-tmppfiles-XXXX", DIR=>'.', CLEANUP=>0 );

    my $envpfiles = $ENV{PFILES};
    croak "cannot get \$PFILES variable\n" unless $envpfiles;

    chomp($envpfiles);
    my @pfdirs = split(/;/,$envpfiles);
    my $firstpf = shift @pfdirs;

    ###system("cp ".$firstpf."/* $tmpdir"); # actually not needed

    $self->{DIR} = $tmpdir;
    $self->{ENV} = join(';',($tmpdir,@pfdirs));


    bless($self,$class);

    return($self);
}

sub env {
    my $self = shift;
    return('PFILES="'.$self->{ENV}.'"');
}

sub DESTROY {
    my $self = shift;

    system('rm -rf '.$self->{DIR});
}

1;
