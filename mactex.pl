use strict;
use 5.10.0;

while(<>){
	chomp;
	say if (/\$\$/../\$\$/);
}
