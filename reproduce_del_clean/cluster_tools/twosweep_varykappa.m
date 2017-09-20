#!/usr/bin/perl
use Class::Struct;
#use feature qw(switch);
#use feature "switch";
#use Switch;

my $numpointsa = $ARGV[0];  # number of points first param is varied
my $numpointsb = $ARGV[1];  # number of steps to vary second param
my $basename = $ARGV[2];    # name of script to generate... will generate $basename1.sh, $basename2.sh, ...
my $commandname = $ARGV[3]; # name of the simulation executable

# THIS CODE USES TYPICAL SETTINGS FOR A GRID ENGINE JOB
# PREAMBLE - ADD YOUR EMAIL FOR NOTIFICATIONS, ETC...
my $preamble1 = '#!/bin/sh
#$ -cwd 
#$ -j y
#$ -M YOUR_EMAIL@EMAIL.EDU
#$ -m sea
'; 

# IN THIS CODE BLOCK, EDIT THE LINE MCRROOT=/home/username/MCR/v80
# SO THAT IT POINTS TO YOUR MCR DIRECTORY
my $preamble2 = '#$ -pe mpich 2
# wall clock time:
#$ -l h_rt=48:00:00
MCRROOT=/home/username/MCR/v80
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export LD_LIBRARY_PATH;
export XAPPLRESDIR;
echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
';

#print "$numruns \n";
#print "$basename \n";
#print $preamble1;
#print $preamble2;

struct VarIn =>
{
    name  => '$',
    range => '@',
    step  => '$',
    factor => '$',
};


# USE THIS SECTION TO CHOOSE WHICH PARAMETERS ARE VARIED
# PARAMETERS WITH STEP=>0 ARE KEPT CONSTANT
# PARAMETERS WITH STEP=>1 or STEP=>3 ARE VARIED
# THE FILENAME BASE, WITH STEP=>2, HAS THE RUN NUMBER APPENDED
my @va;
$va[0] = VarIn->new( name => 'savename',    range=> ['vary_kappa_tavg_'], step =>2);  
$va[1] = VarIn->new( name => 'Nits',     range=> [600], step =>0);  
$va[2] = VarIn->new( name => 'discardfactor', range=> [2], step =>0);  
$va[3] = VarIn->new( name => 'N', range=> [50], step =>0);  
$va[4] = VarIn->new( name => 'dt', range=> [0.02], step =>0);  
$va[5] = VarIn->new( name => 'steps', range=> [20000], step =>0);  
$va[6] = VarIn->new( name => 'L', range=> [15], step =>0);  
$va[7] = VarIn->new( name => 'plotevery', range=> [20], step =>0);  
$va[8] = VarIn->new( name => 'tau', range=> [1], step =>0);  
$va[9] = VarIn->new( name => 'svar', range=> [0.3], step =>0);  
$va[10] = VarIn->new( name => 'kappa', range=> [1,8], step =>3);  
$va[11] = VarIn->new( name => 'lzero', range=> [1], step =>0);
$va[12] = VarIn->new( name => 'tavg',    range=> [0.2,100], step =>1);  
$va[13] = VarIn->new( name => 'Dtheta', range=> [1], step =>0);  
$va[14] = VarIn->new( name => 'Sgrad', range=> [0.025], step =>0);  

print "Varying 2 variables, $numpointsa x $numpointsb \n";

$numruns = $numpointsa*$numpointsb;


my $iplus = 0;

for($p=0; $p < $numpointsa; $p++)
{
    for($q=0; $q < $numpointsb; $q++) 
    {
    $iplus = $iplus + 1;   # count for number of shell files

    open(FOUT,'>'."$basename$iplus.sh");
    print FOUT $preamble1;
    print FOUT '#$ -e t'.$iplus.".err \n";
    print FOUT '#$ -o t'.$iplus.".out \n";
    print FOUT $preamble2;
    my $counter_fac = 0;
    for($i = 0; $i < @va; $i++) {
	my $vv = $va[$i];
	my $ss = $vv->step;
	my $nn = $vv->name;

	if($ss==0) { 
	    my $rr = $vv->range->[0];
	    print FOUT "$nn=$rr \n"; 
	}
	elsif($ss==1) {
	    my $rr = $vv->range->[0] + $p*($vv->range->[1]-$vv->range->[0])/($numpointsa-1);
	    print FOUT "$nn=$rr \n"; 
	}
	elsif($ss==2) {
	    my $rr = $vv->range->[0];
	    print FOUT "$nn=\'$rr$iplus\' \n";
	}
        elsif($ss==3) {
	    my $rr = $vv->range->[0] + $q*($vv->range->[1]-$vv->range->[0])/($numpointsb-1);
	    print FOUT "$nn=$rr \n";     
	}
}

my $fs = "./$commandname ";
#print FOUT './' ,$commandname," ";
for($i = 0; $i < @va; $i++) {
    my $vv = $va[$i];
    my $ss = $vv->step;
    my $nn = $vv->name;
    $fs.="\$$nn ";
    #fs.="$name ";
}
print FOUT "$fs \n";
close(FOUT)

}  # the loop over 


}  #
