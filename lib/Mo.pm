package Mo;
$VERSION=0.38;

# Mo.pm -- (minimalist object constructor)
# https://metacpan.org/pod/Mo
#
# Copyright (c) 2011-2013. Ingy döt Net.
#
# This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.
#
# See http://www.perl.com/perl/misc/Artistic.html
#
# see github repo for real (non-golfed) source: https://github.com/ingydotnet/mo-pm/tree/master/src

no warnings;my$M=__PACKAGE__.'::';*{$M.Object::new}=sub{my$c=shift;my$s=bless{@_},$c;my%n=%{$c.::.':E'};map{$s->{$_}=$n{$_}->()if!exists$s->{$_}}keys%n;$s};*{$M.import}=sub{import warnings;$^H|=1538;my($P,%e,%o)=caller.'::';shift;eval"no Mo::$_",&{$M.$_.::e}($P,\%e,\%o,\@_)for@_;return if$e{M};%e=(extends,sub{eval"no $_[0]()";@{$P.ISA}=$_[0]},has,sub{my$n=shift;my$m=sub{$#_?$_[0]{$n}=$_[1]:$_[0]{$n}};@_=(default,@_)if!($#_%2);$m=$o{$_}->($m,$n,@_)for sort keys%o;*{$P.$n}=$m},%e,);*{$P.$_}=$e{$_}for keys%e;@{$P.ISA}=$M.Object};
