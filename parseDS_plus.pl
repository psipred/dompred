#!/usr/bin/perl

###################################
#	DomPred code that
#	read-in arguments from Java code
#	output results page
#	and alot more between
#
# Modification Record
# 23/02/2004 R. Marsden Created.
# 30/04/2004 K. Bryson  Added configuration parameters so that the
#                       script can be used for a 'stand-alone' test-bed.
# 10/05/2004 K. Bryson  Fixed output problem for 'Putative domain boundaries' table.
# 10/05/2004 K. Bryson  Fixed scale problem on graph
#                       (normalized smhash values used, not plothash values,
#                        so biggest_plothash value should no longer be used).
# 14/05/2004 K. Bryson  Fixed formatting problem for domain boundaries in the DomSSEA
#                       results table (re-positioned <TD>..</TD> tags).
#                       Eliminated output of empty rows when number of hits is less than 10.
# 14/05/2004 K. Bryson  Commented out tests on alpha/beta content for the query/match sequences.
#                       According to Liam these are not needed for DomSSEA output.
# 19/05/2004 K. Bryson  Output of DPS multi-domain boundaries in ascending order.
#                       Output boundary predictions by DomSSEA which are close to
#                       the N-/C-termini as values in red, rather than '?' symbols.
#                       Changed the 'flattening' around detected boundaries from +/-20aa to +/-40aa.
# 20/05/2004 K. Bryson  Output overall summary of predictions by DPS, even when only 1 domain is predicted.
#                       This simplifies parsing of the results for the email server and also
#                       clarifies the overall result produced by DPS for a novice user.
# 28/05/2004 K. Bryson  Changed 'flattening' to be +/- 80aa now.
#                       Fixed bug in output of number of domains
#                       (Uses $num_peaks now instead of $peaks_detected).
# 04/06/2004 K. Bryson  Added configurational parameter for the database used.
#                       Created a new uniref90_Pfam_A_nofrag database.
# 07/06/2004 K. Bryson  After examining the results of some benchmarking,
#                       put back the old database ! Probably optimally configured for this.
# 07/06/2004 K. Bryson  Extended the end-effect flattening to +/-70aa (plus another 10aa possibly ?).
# 15/06/2004 K. Bryson  Added '0' to pdb identifiers for SCOP file indexing.
# 17/06/2004 K. Bryson  Fixed bug calculating the variance in the z-score section (crashed CASP6 T203 target).
#                       TEMP FIX - Just take the abs(variance) value --- using the correct value
#                       'throws out' all the other calibration and domain prediction becomes unreliable.
# 17/06/2004 K. Bryson  Sorted problems with graphing code (flattening smhash, using '-' instead of 0).
# 17/06/2004 K. Bryson  Tidied up code by removing completely unused sections of code.
# 17/06/2004 K. Bryson  Took out 'helix' and 'strand' profiles from graph since these are distracting from
#                       the real graph of aligned termini profile and are undocumented in the help file.
#                       Replaced them with raw output of helix and sheet regions as points, similar to that for coil.
#                       This is cleaner, less distracting, and provides similar information for domain prediction.
# 18/06/2006 K. Bryson  Updated for CASP7 to use the latest uniref90_PfamA databsae.
#
# 06/04/2008 D. Buchan  Additionally added an output file for the DOMssea lines that are selected
#						Took out that nonsense of wrapping the blast files in <pre> tags so that they are "html"
#						(note that .blastdom and blast.bls are identical.  blast.bls is then gzipped for passing into the
#						rails db).
#						Changed the output image to a 16bit png and removed the use of jpegtran. The image is now rotated 
#						by mogrify in a single rescaling and rotating step. The png is clearer.  it might be worth
#						playing with the rescaling filter to make the image look nicer.
#						Added a couple of output flat files - .boundary for the DPS info
#															.domsseatable for the domssea results that will be reported on the webpage
#						Removed all code that writes out to any html files. This should be handled by the view in Rails not at the backend.
#						PSI-Blast is no longer run within this script, the rails backend now handles this (possible candidate for command line switch though).
#						Added too many commandline options to remove any hardcoded paths.
# 23/09/2012 D. Buchan Assorted minor changes to the .blastdom parsing to handle the changes in blast+ output
#
##################################
#
# TODO:
#	Suppress the calls to gs mogrify and have the java backend run those - annoyingly it does all sorts of stuff about testing
#	where the ps and png images are which would also need to be written out.
# 	Have the input data for the gnuplots output to files and have the java backend run gnuplot
#
##################################


#use strict;

###
#       Configuration parameters.
###
# Original configuration for web installation.
#
#$OUT_DIR       ="/webdata/tmp/NewPredServer/";
#$DOMSSEA_HOME  ="/webdata/data/current/dompred/";
#$GNUPLOT_EXEC  ="/usr/bin/gnuplot";
#$BLAST_EXEC    ="/webdata/binaries/current/blast/bin/blastpgp";
#$JPEGTRAN_EXEC ="/usr/bin/jpegtran";

#$DATABASE_NAME ="uniref90_PfamA"; # ORIGINAL DATABASE.
# $DATABASE_NAME ="nrdb90_Pfam_A_nofrag"; # ORIGINAL DATABASE.
# $DATABASE_NAME ="uniref90_Pfam_A_nofrag";

# Test-bed configuration.
#$OUT_DIR        ="/scratch1/NOT_BACKED_UP/ucackxb/projects/Russell/WORK/CASP4_T0087/REPEAT_CASP4_TO087";
#$DOMSSEA_HOME   ="/scratch1/NOT_BACKED_UP/ucackxb/projects/Russell/software/others/marsden/DomSSEA";
#$GNUPLOT_EXEC   ="/scratch1/NOT_BACKED_UP/ucackxb/projects/Russell/software/others/marsden/gnuplot-3.7.1/gnuplot";
#$BLAST_EXEC     ="/scratch1/NOT_BACKED_UP/ucackxb/projects/Russell/BLAST/blastpgp";
#$JPEGTRAN_EXEC  ="/scratch1/NOT_BACKED_UP/ucackxb/projects/Russell/JPEGTRAN/jpegtran";
#$DATABASE_NAME ="nrdb90_Pfam_A_nofrag";

# Parameter for the structure DB being used.
$STRUCTURE_DB  ="SCOP";
$helix_count = 0;
$extended_count = 0;
$coil_count = 0;
	

###
#	get name of input sequence
###
$input_path=$ARGV[0] or die "Usage parseDS <path/sequence_name>\n";# path of query sequence
$do_seq=$ARGV[1];			#1 	seqyes seqno (seq dom align)
$do_database=$ARGV[2];		#2	nonrucl pfama
$do_evalue=$ARGV[3];		#3	0.001 (evalue)
$do_its=$ARGV[4];			#4	5	(iterations)
$do_domssea=$ARGV[5];		#5	IGNORE
$do_ppred=$ARGV[6];			#6	ppyes ppno
$do_secprofile=$ARGV[7];	#7	secproyes secprono
$sequence_name=$ARGV[8];	#8	sequencename
$DATABASE_NAME=$ARGV[9];	#9	The blast database of domains that the predictions are based on (uniref90_PfamA)
$OUT_DIR=$ARGV[10];			#10 the temp directory for all files.
$DOMSSEA_HOME=$ARGV[11];	#11 The location of the Domssea code
$GNUPLOT_EXEC=$ARGV[12];	#12 The path to the installation of gnuplot


if($sequence_name eq ""){$sequence_name="?";}	 	#if no given name, call ?
$sequence_name2=substr($sequence_name,0,15);	#truncate nme to <=15 chars
###
#	add below security to front page code, and ask for more data if invalid
###
if($do_evalue==0) {$do_evalue=0.001;}						#default e_value of 0.001
$do_evalue=~s/\s+//;														#remove spaces
if(($do_evalue=~/^\d+\.\d+$/)||($do_evalue=~/^\d+$/)||($do_evalue=~/^\.\d+$/)) #if not an expected number type
	{}
else{$do_evalue=0.001;}			#make default
$do_its=~s/\s+//;						#remove spaces
if($do_its=~/\w+/g){$do_its=5;} #ignore non number chars
if($do_its<1) {$do_its=1;}		#make smallest it can be
if($do_its>10) {$do_its=5;}		#default iterations

if($input_path=~/(.*)\/(\S+)\.domssea$/) #parse out query name from given path
{
	$input_seq=$2;
}
$target=$input_seq;		#for use by difference in length?
###
#	this is the main output file
###
open (BOUNDARYOUT, ">$OUT_DIR/$input_seq.boundary") or die "Can't create output boundary file\n";

###																  eg >/www/htdocs/domout/17_18_59_238.pred
#	raw secondary struc output file		ACCCCCSCSSDDHHGHHSHSHSHSHAHAHAHAH
###																	CEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCC
open(SECSTRUC, "$OUT_DIR/$input_seq.pred") or die "can't open sec_pred $OUT_DIR $input_seq\n";
my $sec_pred=<SECSTRUC>; #first line
my $sec_pred=<SECSTRUC>; #second line
$sec_pred=<SECSTRUC>;		 #third line
close(SECSTRUC);
@sec=split(//,$sec_pred);
$length=scalar(@sec);	 #can get length of sequence from this
$ita=0;
foreach $aminoacid (@sec)	#this sets up the coil, helix, beta plot for gnuplot
	{
	if($sec[$ita]=~/C/){$coilhash{$ita}=0.1;}
	else{$coilhash{$ita}="-";}

	if($sec[$ita]=~/H/){$helixhash{$ita}=0.1;}
	else{$helixhash{$ita}="-";}

	if($sec[$ita]=~/E/){$extendedhash{$ita}=0.1;}
	else{$extendedhash{$ita}="-";}

	$ita++;
	}
###
#	open alpha beta percentage data and read in
###	
open(SECSTRUC, "$DOMSSEA_HOME/2struc_ab_percent") or die;
while ($line=<SECSTRUC>)
	{
	if($line=~/(\S+)\s+(\d+.\d+)\s+(\d+.\d+)/)
		{
		$alphahash{$1}=$2;
		$betahash{$1}=$3;
		}
	}
close(SECSTRUC);
###################################
#	Sequence profile code
# will also plot as gnuplot 
###################################

#if($do_seq eq "seqyes") #need to sort this
#	{
	&sequenceX;
#	}
###
#run sequence profile algorithm
###
sub sequenceX 
	{
###
#do secondary structure plot
###
###
#	calculate percentage of alph and beta residues
###
	for ($abperc=1; $abperc<=(scalar(@sec)); $abperc++)
  	{
		if(@sec[$abperc]=~/H/){$ares++;}
		if(@sec[$abperc]=~/E/){$bres++;}
		}
	$freq_alpha=$ares/$length;
	$freq_beta=$bres/$length;
	
	&sec_struc; 
	sub sec_struc
		{
		#start with 30 res window from start to end-30
		for ($secj=1; $secj<=(scalar(@sec) - 29); $secj++) #$C - 59....
  		{
  		$counterx = $countery = $counterz = 0;
			for ($secl=$secj; $secl<=($secj+29); $secl++)
    		{
    		$element = @sec[$secl-1];
    		if ($element =~ /E/){ $counterx++; }		#beta
    		elsif ($element =~ /H/){ $countery++; }	#alpha
    		elsif ($element =~ /C/){ $counterz++; } #coil
    		elsif ($element =~ //) { $counterz++;	} #nothin due to missin bit
    		}
			$hashb{$secj+14} = sprintf("%.2f",($counterx/30));#beta
			$hasha{$secj+14} = sprintf("%.2f",($countery/30));#alpha
			#$hashc{$secj+14} = sprintf("%.2f",($counterz/30));#coil
			}
		}
		#collect test sequence in a cgi way
		#run seq against psi blast, num iterations and cutoff params
		#filter fragments and output plot.gif

	$num_iterations=$do_its; #number of given itterations
###
#	run psiblast
###		
	#system ("$BLAST_EXEC -i  $OUT_DIR/$input_seq -j $num_iterations -a 2 -m 0 -b 1000 -d $DOMSSEA_HOME/$DATABASE_NAME  > $OUT_DIR/$input_seq.blastdom");
	
	$with_eval=1; #forgot what this is for
	$output_blast=$input_seq .".blast";	#thats the name for the blastout html file
	$output_blastsee=$input_seq .".blastsee";	#name for the parsed blast output data
	$output_blastPfamA=$input_seq ."PfamA";
	open (SEE, ">$OUT_DIR/$output_blastsee.txt") or die "Cannot open blastsee_html\n";
	&printtoparsed;#start table for parsed blast
	open (TEXT, ">$OUT_DIR/$output_blast.bls") or die "Cannot open blast_html\n";
	#print TEXT "<PRE>\n";
	open (PFAMA, ">$OUT_DIR/$output_blastPfamA.html") or die "Cannot open blast_html\n";
	#print PFAMA "<PRE>\n";
	$interval=1;
###
#	read in blast out and parse it
###		
	print "opening"."$OUT_DIR/$input_seq.blastdom\n";
	open (INFIL, "$OUT_DIR/$input_seq.blastdom") or die "Cannot open $input_seq\n";
	while ($line=<INFIL>)
		{
		print TEXT "$line";
		if($line=~/\((\d+)\s+letters\)/)
			{
			$length=$1;
			}
		if($line=~/Searching/)
			{
			#open (SEE, ">$OUT_DIR/$output_blastsee.html") or die "Cannot open blastsee_html\n";
			#&printtoparsed;#start table for parsed blast
			$iteration++;
			$curr_it=$1;
			$fragment_seen=1; #use this to stop previous info being used
			}
		if($line=~/CONVERGED!/)
			{
			#open (SEE, ">$OUT_DIR/$output_blastsee.html") or die "Cannot open blastsee_html\n";
			#&printtoparsed;#start table for parsed blast
			$curr_it=0;	
			$fragment_seen=1; #use this to stop previouse info being used
			}
		if(($line=~/^>(.*)/)||($line=~/Matrix: BLOSUM62/)) #next alignment, print out last if seen
			{
			if(($initiated==1)&&($fragment_seen==0)) #ie seen a > and last align wasn't a fragment
				{
				if($pfamA_hit==1)
					{
					#$pfamA_hit=0;
					$colour="blue";
					}
				for($com=1;$com<=$score_see;$com++)
					{
					$querylist="querylist" . $com; #for query ends
					$querylistx="querylistx" . $com;	#for db-hit ends
					$hitlength=shift @$querylist;
					$e_value=shift @$querylist;
					if($e_value<=$do_evalue)				#print out stats
						{
						$howmany{$iteration}++;
						$name2=substr($name,0,35);		#name of last one seen
						
						if(length($name2)<35)
							{
							until (length($name2)==35)
								{
								$name2=$name2 . " ";
								}
							}
						
						if($e_value==0.0){$e_value = "1e-200";}
						$mod_e_value=(log(1/$e_value));#/log(10); #rm if not want to log eval
						#print SEE "<TR>\n";
						#print SEE "<TD>$name2</TD>\n";
						#print SEE "<TD>$length</TD>\n";
						#print SEE "<FONT color = \"$colour\">\n";
						print SEE "$name2	";
						print SEE "$length	";
						$startnumber=shift @$querylist;
						$endnumber=pop @$querylist;
						$hitstart=shift @$querylistx;
						$hitend=pop @$querylistx;
						#print SEE "<TD>$startnumber</TD>\n";
						#print SEE "<TD>$endnumber</TD>\n"; 
						#print SEE "<TD>$hitlength</TD>\n";
						#print SEE "<TD>$hitstart</TD>\n";	
						#print SEE "<TD>$hitend</TD>\n";	
						#print SEE "<TD>$e_value</TD>\n";
						#print SEE "</TR>\n";
						
						print SEE "$startnumber	";
						print SEE "$endnumber	"; 
						print SEE "$hitlength	";
						print SEE "$hitstart	";  
						print SEE "$hitend	";  
						print SEE "$e_value	";
						print SEE "\n";
						#print SEE "</FONT>";
						if($pfamA_hit==1)
							{
							#print PFAMA "<FONT color = \"$colour\">\n";
							print PFAMA "$name2	";
							print PFAMA "$length	";
							print PFAMA "$startnumber ";
							print PFAMA "$endnumber "; 
							print PFAMA "$hitlength ";
							print PFAMA "$hitstart  ";  
							print PFAMA "$hitend  ";  
							print PFAMA "$e_value ";
							print PFAMA "\n";
							#print PFAMA "</FONT>";
							}
						$pfamA_hit=0;
						$colour="black";
						if(($lengthhit-$hitend)<10)	#ie C term
							{
							$hash{$iteration}{$endnumber}++;#=$mod_e_value;
							$hashC{$iteration}{$endnumber}++;
							}	
						if($hitstart<10)
							{
							$hash{$iteration}{$startnumber}++;#=$mod_e_value;
							$hashN{$iteration}{$startnumber}++;
							}
						$hash{$iteration}{$startnumber}++;
						$hash{$iteration}{$endnumber}++;
						$hashN{$iteration}{$startnumber}++;
						$hashC{$iteration}{$endnumber}++;						
						@$querylist=();
						@$querylistx=();
						}
					}
				$score_see=0;
				}
			if($line=~/^>(.*)/) #get name of this one
				{					
				$name=$1;
				if($line=~/^>(\S+)\s+PF\d+\;/)
					{
					$pfamA_hit=1;
					$pfamA_seen=1;
					}
				for($com=1;$com<=$score_see;$com++)
					{
					$querylist="querylist" . $com;
					$querylistx="querylistx" . $com;
					@$querylist=();@$querylistx=();
					}
				$score_see=0;
				}
			$initiated=1;
			if($line=~/(FRAGMENT)|(fragment)|(Fragment)/)
 				{
				$fragment_seen=1;
				next;
				}
			if($line!~/(FRAGMENT)|(fragment)|(Fragment)/)
 				{
				$fragment_seen=0;  #reinitialise
				}
			}
			
		if($with_eval==1)
			{
				#D.Buchan: Assorted minor changes to parse the new blast+ output
				if($line=~/Length=(\d+)/)
				#if($line=~/Length\s+=\s+(\d+)/)
				{
					$hitlength=$1;
				}
				if($line=~/Score\s+=\s+(.*)\s+bits/)
				{
###
#matching more than once to query --> $score_see >1
###					
					$score_see++;
					$querylist="querylist" . $score_see;
					$querylistx="querylistx" . $score_see;
					$bits=$1;
					push @$querylist,$hitlength;
				}
				if($line=~/Expect\s+=\s+(\d+\.\d+)/)
				{
					$expecting++;
					$e_value=$1;
					push @$querylist,$e_value;
				}
				if($line=~/Expect\s+=\s+(e\S+)/)
				{
					$expecting++;
					$e_value="1" . $1;
					push @$querylist,$e_value;
				}	
				if($line=~/Expect\s+=\s+(\d+e\S+)/)
				{
					$expecting++;
					$e_value=$1;
					push @$querylist,$e_value;
				}	
				if($line=~/^Query\s\s+(\d+).*\s+(\d+)$/)
				{
					$start=$1; $end=$2;
					push @$querylist, $start; push @$querylist, $end;
				}
				if($line=~/^Sbjct\s\s+(\d+).*\s+(\d+)$/)
				{
					$start=$1; $end=$2;
					push @$querylistx, $start; push @$querylistxs, $end;
				}
			}
		}
	#print SEE "</PRE>\n";
	close (SEE);
	#print TEXT "</PRE>\n";
	close (TEXT);
	#print PFAMA "</PRE>\n";
	close (PFAMA);
###
#	blast has been read in and parsed
#
#	output / calculate plot data (need to sort)	
###
# DON'T NEED THIS AS FAR AS I CAN TELL??
# THIS IS LIKE DONE BEFORE EVERYTHING IS FLATENED, 
 $biggest_plothash=0;
 for($i=1;$i<=($length-14);$i++)
   {                             #so plot from +45 to -55
   for($j=$i;$j<($i+14);$j++)     #or do 36 to 45
     {
     $plothash{$i+7}+=$hash{$iteration}{$j};    #get total in window -> middle res
     }
   $plothash{$i+7}=$plothash{$i+7}/15;                    #average smoothing
   if($plothash{$i+7}>$biggest_plothash){$biggest_plothash=$plothash{$i+7};}
   }
 if($biggest_plothash==0){$biggest_plothash=1;}

###
#	analyse psiblast profile
###

### KB (07/06/2004) - Extend this to +/-70aa at boundaries (plus another 10aa possibly ?).

	$bigpoint=0;
	for($i=1;$i<70;$i++)	
		{
		$hash{$iteration}{$i}=0;
		$hashN{$iteration}{$i}=0;
		}
	for($i=($length-70);$i<=$length;$i++)
		{
		$hash{$iteration}{$i}=0;
		$hashC{$iteration}{$i}=0;
		}
	#if 'hits' on boundaries rm.... don't want to incorparate into Zscore 
	for($yi=70;$yi<=80;$yi++) #10 or 20?
		{
		if($hashC{$iteration}{$yi}>0)
			{
			$seenC=1;#have seen a C-term hit here, (so keep)
			}
		}
	if($seenC!=1)	#not = 1 ie not seen any c-term hit
		{
		for($yi=70;$yi<=80;$yi++) #10 or 20?
			{
			$hash{$iteration}{$yi}=0;
			$hashN{$iteration}{$yi}=0;
			}
		}
	for($yi=($length-70);$yi>=($length-80);$yi--)
		{
		if($hashN{$iteration}{$yi}>0)
			{
			$seenN=1;#have seen a N-term hit here, so keep)
			}
		}
	if($seenN!=1)
		{
		for($yi=($length-70);$yi>=($length-80);$yi--)
			{
			$hash{$iteration}{$yi}=0;
			$hashC{$iteration}{$yi}=0;
			}		
		}
###
#	smoothing window of 15 residues for N and then C hits
###	
	for($i=1;$i<=$length-14;$i++)
		{															#so plot from +55 to -65
		for($j=$i;$j<($i+15);$j++)			
			{
			$smhashN{$i+7}+=$hashN{$iteration}{$j};		#get total in window -> middle res
			}
		$smhashN{$i+7}=$smhashN{$i+7}/15; 					#average window value
		}
	for($i=1;$i<=$length-14;$i++)
		{															#so plot from +55 to -65
		for($j=$i;$j<($i+15);$j++)			
			{
			$smhashC{$i+7}+=$hashC{$iteration}{$j};		#get total in window -> middle res
			}
		$smhashC{$i+7}=$smhashC{$i+7}/15; 					#average window value
		}
###
#check out N and C matching in above smoothed profile
###
	for($i=8;$i<=$length-7;$i++) #8 cos thats where 15 window starts
		{	
		if(($smhashN{$i}*$smhashC{$i})==0)
			{	
			$smhash{$i}=$smhashN{$i}+$smhashC{$i};
			}
		if(($smhashN{$i}*$smhashC{$i})>0)
			{	
			if($smhashN{$i}>$smhashC{$i}){$extra=$smhashC{$i}/2;}
			if($smhashC{$i}>=$smhashN{$i}){$extra=$smhashN{$i}/2;}
			$smhash{$i}=$smhashN{$i}+$smhashC{$i}+$extra;
			}
		}
###
#	calculate Zscores
###	
	foreach $x (sort {$a<=>$b} keys %smhash)
		{
		$total_res+=$smhash{$x};
		$total_seen++;
		$v2=($smhash{$x}*$smhash{$x});
		$vtot2+=$v2;
		if($smhash{$x}>$bigpoint){$bigpoint=$smhash{$x};}	#get higherst peak
		}
	if($total_seen==0){next;} #important, ie no blast profile
	$mean=$total_res/$total_seen;
	$sigx2=(($total_res*$total_res)/$total_seen);

        # 17/06/2004 KB Fixed bug - division by $total_seen-1,
        #                           rather than $total_res-1;
        # TEMP FIX TEMP FIX - Just use abs(variance).
	$variance=($vtot2-$sigx2)/($total_res-1);
	$sd=sqrt abs($variance);
	foreach $e (sort {$a<=>$b} keys %smhash)
		{
		if(($smhash{$e}==0)&&($mean==0)&&($sd==0))
			{
			$Zhash{$e}=0;
			}
		else
			{
			$Zscore=($smhash{$e}-$mean)/$sd;
			$Zhash{$e}=$Zscore;
			}
		}

###	
#	iterate over Zscored profile
###
	$bigZ=$g=1.5; #Zscore of 1.5
	$peak=1;
	$num_peaks=0;
	$reallybig=0;
	while($peak==1)
		{
		$bigZ=$bigZpeak=$g;
		$peak=0;
		foreach $x (sort {$a<=>$b} keys %smhash)
			{
			if($Zseen{$x}==1){next;}		#so don't see it on next iteration
			if($bigZpeak<$Zhash{$x}) 		#find highest peak first
				{
				$bigZpeak=$Zhash{$x};	
				if($bigZpeak>$reallybig) 		
					{
					$reallybig=$bigZpeak;
					}
				$peak=1;
				$Zres=$x;
				}
			}
		if($peak==1)
			{
			### sort this to store predicted boundaries
			$Zseen{$Zres}=1;
			$num_peaks++;
			push @termini_preds, $Zres;

			# Flatten the smhash/Zhash with +/-80aa around detected boundary.
			for($flat=($Zres-80);$flat<=($Zres+80);$flat++)
				{
				if($flat==$Zres){next;}

				# 17/06/2004 KB Surely should not flatten smhash,
				# otherwise graph shown will have data flattened.
				# $smhash{$flat}=0;
				$Zhash{$flat}=0;
				}
			}
		}	
	if($num_peaks==0)
		{
		}
	if($num_peaks>0)
		{
		$peak_detected++;
		# $num_peaks=0; # Use $num_peaks.
		}
###
#	to normalise this plot, divide all points by highest value!#
###
	if($bigpoint==0)	#ie no psiblast hits make all plot 0
		{
		foreach $x (sort {$a<=>$b} keys %smhash)
			{
			$smhash{$x}=0;
			}
		}
	if($bigpoint>0)
		{
		foreach $x (sort {$a<=>$b} keys %smhash)
			{
			$smhash{$x}=$smhash{$x}/$bigpoint;
			}
		}
        }
###
#	output profile plot and secondary structure plot, here use smhash (which is raw data?)
#	could also add hashN and hashC and smhash
###


&shownorm;

###
#    plot, again with smhash.
#    17/06/2004 KB Rewrote this section.
###
	#$helix_count = 0;
	#$extended_count = 0;
	#$coil_count = 0;
	
	sub shownorm
	{

	    open (OUT, ">$OUT_DIR/$input_seq.graph") or die "can't create outfile\n";
	    foreach $x (sort {$a<=>$b} keys %smhash)
		{
			print OUT "$x $smhash{$x} $helixhash{$x} $extendedhash{$x} $coilhash{$x}\n";
			#print "THINGY: ".$helixhash{$x}."\n";
			if($helixhash{$x} =~ /\d/)
			{
				$helix_count++;
			}
			if($extendedhash{$x} =~ /\d/)
			{
				$extended_count++;
			}
			if($coilhash{$x} =~ /\d/)
			{
				$coil_count++;
			}
	    }
		close(OUT);
	}
##
#   End of shownorm subroutine.
##


	if($howmany{$iteration} eq "")
	{
		$howmany{$iteration}=0; #eh? should do this sooner?
	}
		#$howmany{$iteration}=1;
	if(($do_secprofile eq "secproyes")&&($howmany{$iteration}>0)) #$do_secprofile=$ARGV[7];	#7	secproyes secprono
	{
		&gnuplot1; #plot 2' and psiblast
		&pstogif;
	 }
	if(($do_secprofile eq "secprono")&&($howmany{$iteration}>0))
	{
	 	&gnuplot2;	#plot only psiblast
	 	&pstogif;
	}
	if(($do_secprofile eq "secproyes")&&($howmany{$iteration}==0))
	{
	 	&gnuplot3;	#plot only 2'
	 	&pstogif;
	}
###
#	gnuplot commands
###

	#7/4/2009  D.Buchan added some globals that track how many helix, strand or coil residues are to be plotted and if there are none that line is NOT
	#added to the gnuplot script
	sub gnuplot1 
	{
 		### save tmp GNUPLOT file
 		open (GNUPLOT, ">$OUT_DIR/$input_seq.gnu") or print("Couldn't save gnuplot file: $gnufile");
 
 	 	print GNUPLOT ("set terminal postscript color solid \"Times-Roman\" 18\n");
		$myps=$input_seq . "myps";
  		print GNUPLOT ("set output \"$OUT_DIR/$myps.ps\"\n");
		#print GNUPLOT ("set yrange [0:]\n");
		
		print GNUPLOT ("set format y \"\"\n");
		#print GNUPLOT ("set format y2 \"\"\"\n");

		# print GNUPLOT ("set y2range [0:$biggest_plothash]\n");
		print GNUPLOT ("set y2range [0:1.2]\n");
 	 	print GNUPLOT ("plot \"$OUT_DIR/$input_seq.graph\" using 1:2 axes x1y2 title 'aligned termini profile' with lines lt 3");

		# 17/06/2004 KB Changed the alpha-helix, beta-sheet plots into output similar to coil.
		if($helix_count != 0)
		{
			print GNUPLOT (", \\\n\"$OUT_DIR/$input_seq.graph\" using 1:3 axes x1y2 title 'helix-residues'  with points lt 2 pt 5");
		}
		if($extended_count != 0)
		{
			print GNUPLOT (", \\\n\"$OUT_DIR/$input_seq.graph\" using 1:4 axes x1y2 title 'strand-residues' with points lt 1 pt 5");
		}
		if($coil_count != 0)
		{
			print GNUPLOT (", \\\n\"$OUT_DIR/$input_seq.graph\" using 1:5 axes x1y2 title 'coil-residues'   with points lt -1 pt 5");
  		}
		print GNUPLOT ("\n");
		print GNUPLOT ("exit\n");
  		close(GNUPLOT);
 		$command="$GNUPLOT_EXEC $OUT_DIR/$input_seq.gnu";
  		system("$command");
	}	
	
	sub gnuplot2 
	{
 		### save tmp GNUPLOT file
 		open (GNUPLOT, ">$OUT_DIR/$input_seq.gnu") or print("Couldn't save gnuplot file: $gnufile");
 
 	 	print GNUPLOT ("set terminal postscript color solid \"Times-Roman\" 18\n");
		$myps=$input_seq . "myps";
  		print GNUPLOT ("set output \"$OUT_DIR/$myps.ps\"\n");
 	 	print GNUPLOT ("plot \"$OUT_DIR/$input_seq.graph\" using 1:2 title 'aligned termini profile' with lines lt 3\n");
  		print GNUPLOT ("exit\n");
  		close(GNUPLOT);
 		$command="$GNUPLOT_EXEC $OUT_DIR/$input_seq.gnu";
  		system("$command");
	}	
	
	sub gnuplot3 
	{
 		### save tmp GNUPLOT file
 		open (GNUPLOT, ">$OUT_DIR/$input_seq.gnu") or print("Couldn't save gnuplot file: $gnufile");
 
 	 	print GNUPLOT ("set terminal postscript color solid \"Times-Roman\" 18\n");
		$myps=$input_seq . "myps";
  		print GNUPLOT ("set output \"$OUT_DIR/$myps.ps\"\n");
 	 	#print GNUPLOT ("plot \"$OUT_DIR/$input_seq.graph\" using 1:2 title 'aligned termini profile' with lines lt 3, \\\n");

		# 17/06/2004 KB Changed the alpha-helix, beta-sheet plots into output similar to coil.
		print GNUPLOT ("plot \"$OUT_DIR/$input_seq.graph\" using 1:3 axes x1y2 title 'helix-residues'  with points lt 2 pt 5, \\\n");
		print GNUPLOT ("\"$OUT_DIR/$input_seq.graph\" using 1:4 axes x1y2 title 'strand-residues' with points lt 1 pt 5, \\\n");
		print GNUPLOT ("\"$OUT_DIR/$input_seq.graph\" using 1:5 axes x1y2 title 'coil-residues'   with points lt -1 pt 5 \n");
  		print GNUPLOT ("exit\n");
  		close(GNUPLOT);
 		$command="$GNUPLOT_EXEC $OUT_DIR/$input_seq.gnu";
  		system("$command");
	}	
###
#	make plot correct size / format
###
	sub pstogif
		{
		system("gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -sOutputFile=$OUT_DIR/$myps.png -c save pop -f $OUT_DIR/$myps.ps	");
		$inputran=$input_seq . "tran";
		#system("$JPEGTRAN_EXEC -rot 90 $OUT_DIR/$myps.jpg > $OUT_DIR/$inputran.jpg");
		system("mogrify -resize 540x700 -rotate 90 -crop 0x0 -quality 100 $OUT_DIR/$myps.png");
		}#-crop 6%
	#http://irl.eecs.umich.edu/jamin/pointers/gnuplot_quick.html
	#}
# probably don't need to do this, but will just incase
undef %smhash; #average smoothing
undef %Zhash;
undef %hash;
undef %Zseen;
$total_res=0;
$total_seen=0;
$v2=0;
$vtot2=0;


#	Need to sort this out, to make look better / more reliable 0 plot stuff....
#
###
if(($do_secprofile eq "secproyes")||($howmany{$iteration}>0)) #ie can show graph
	{
	&graphpage;
	}
if(($do_secprofile eq "secprono")&&($howmany{$iteration}==0)) #ie can't show graph
	{
	#&graphpage2;
	}

###################################
#	read in DomSSEA out put
#	If peak_detected>0; means 
#	upreg multi domssea preds.......
#	
###################################
@score_array=(0,0,0,0,0,0,0,0,0,0); #initialise @score_array
if(($peak_detected>0))#&&($length>300))	#so saying PsBlast says multi and has to be >300 residues..
	{
	#&domsseapage;	#output domssea bit
	open (PAIRS, "$input_path") or die "Can't open domlist\n";
	while ($pair = <PAIRS>)
		{
		if($pair=~/^(\d+)(\.\d+)?\s+(\S+)\s+(\d+)\s+(.*)\"(.*)\"/) #open up domssea output file
			{
			$score=$1.$2;				#get score
			if($score>1){next;}		#if greater than 1, is wrong...
			}
		while($pair=~/\s+(\S+)\s+(\d+)\s+(.*)\"(.*)\"/g)	#all the stuff after the score
			{
			#$target=$1;			#global?
			$match=$1;		#matched db hit
			$num_doms=$2;	#global?		
			$rest=$3;			#ie all the predicted boundaries
			$sec_align{$match}=$4;	#the 2' structure alignment of this match (could find % a and b)
			$cant=0;			#
			#if(dothis==1)	#not intitialised so not doing
			#	{								#this checks boundary cuts that fall in first/last 40 residues
			#	if($num_doms>1)
			#		{
			#		while($rest=~/(\d+)/g)
			#			{
			#			if(($1<40)||($1>($length-40)))
			#				{				
			#				$cant=1;			
			#				}
			#			}
			#		}
			#	}
			#if($cant==1){next;}
###
#	check out a and b percentage matching...
###
			$freq_alpha_match=$alphahash{$match};
			$freq_beta_match=$betahash{$match};


			# KB - Commented out these tests, not really required for the DomSSEA input.

#			if(($freq_alpha<1)&&($freq_beta_match<1)&&($freq_alpha_match>1))	#target no alpha, match no beta......(so unless all coil..)
#				{
#				next;  # completely diff 2' structures !!!!!!!!!!!
#				}
#			if(($freq_beta<1)&&($freq_alpha_match<1)&&($freq_beta_match>1))	#target no alpha, match no beta......(so unless all coil..)
#				{
#				next;  # completely diff 2' structures !!!!!!!!!!!
#				}
#			if((($freq_alpha<10)&&($freq_alpha_match>10))||(($freq_alpha>10)&&($freq_alpha_match<10)))
#				{
#				next;
#				}
#			if((($freq_beta<10)&&($freq_beta_match>10))||(($freq_beta>10)&&($freq_beta_match<10)))
#				{
#				next;
#				}


			$numdom_hash{$match}=$num_doms;	#record the number of domains for this match
			# this is where single domains are excluded from the domssea output
			if(($peak_detected>0)&&($num_doms==1))
				{
				next;  # as seq pro has shown target to be most likely a multi-domain chain
				}
			if($num_doms>1)
				{
				$count=1;
				while($rest=~/(\d+)/g)	#still going through $rest, 
					{
					$cut_hash{$match}{$count}=$1;	#record domain boundaries
					$count++;
					}
				}
			if($score_hash{$score}==1)  #remember score, so appears only once in array
				{												#ie for multiple hits of same score
				push @{$score}, $match;	#@score keeps hits of that value
				next;
				}
			$score_hash{$score}=1;		#keep record that have seen that score
			@sorted=sort {$b<=>$a} @score_array;	#sort scores
			@score_array=@sorted;								 #sorted array
			if($score>=@score_array[9])	#add new score if bigger than last in list
				{
				push @{$score}, $match;
				@score_array[9]=$score;
				}
			}
		}
	&showdomssea;	#outputs table of hits (together with associated links, cath codes etc)
###	
#	IE do not output anyother DomSSEA data
###
	#$peak_detected=0; 
	#&clearbar;	#bar thingy?
	undef %score_hash;	#cos not benchmarking
	undef $cut_hash;		#just showing data
	}

if($peak_detected==0) #so PsiBlast has said single (ie no peaks detected)
	{
	#&domsseapage;	#output domssea bit
	@score_array=(0,0,0,0,0,0,0,0,0,0);
	@sorted=();
	open (PAIRS, "$input_path") or die "Can't open domlist\n";
	while ($pair = <PAIRS>)
		{
		if($pair=~/^(\d+)(\.\d+)?\s+(\S+)\s+(\d+)\s+(.*)\"(.*)\"/)
			{
			$score=$1.$2;				#get score
			if($score>1){next;}		#if greater than 1, is wrong...
			}
		while($pair=~/\s+(\S+)\s+(\d+)\s+(.*)\"(.*)\"/g)
			{
			#$target=$1;			#global?
			$match=$1;		#matched db hit
			$num_doms=$2;	#global?
			$rest=$3;			#ie all the predicted boundaries
			$sec_align{$match}=$4;	#the 2' structure alignment of this match (could find % a and b)
			$cant=0;
			#if(dothis==1)	#not intitialised so not doing
			#	{							#this checks boundary cuts that fall in first/last 40 residues
			#	if($num_doms>1)
			#		{
			#		while($rest=~/(\d+)/g)
			#			{
			#			if(($1<40)||($1>($length-40)))
			#				{					
			#				$cant=1;			
			#				}
			#			}
			#		}
			#	}
			#if($cant==1){next;}
			
			$freq_alpha_match=$alphahash{$match};
			$freq_beta_match=$betahash{$match};


			# KB - Commented out these tests, not really required for the DomSSEA input.

#			if(($freq_alpha<1)&&($freq_beta_match<1)&&($freq_alpha_match>1))	#target no alpha, match no beta......(so unless all coil..)
#				{
#				next;  # completely diff 2' structures !!!!!!!!!!!
#				}
#			if(($freq_beta<1)&&($freq_alpha_match<1)&&($freq_beta_match>1))	#target no alpha, match no beta......(so unless all coil..)
#				{
#				next;  # completely diff 2' structures !!!!!!!!!!!
#				}
#			if((($freq_alpha<10)&&($freq_alpha_match>10))||(($freq_alpha>10)&&($freq_alpha_match<10)))
#				{
#				next;
#				}
#			if((($freq_beta<10)&&($freq_beta_match>10))||(($freq_beta>10)&&($freq_beta_match<10)))
#				{
#				next;
#				}
#			if(($num_doms==1)&&($length>400)&&($score<7))
#				{
#				next;
#				}

			$numdom_hash{$match}=$num_doms;	#record the number of domains fro this match

			# this is where single domains are excluded from the domssea output
			if($num_doms>1)
				{
				$count=1;
				while($rest=~/(\d+)/g)	#still going through $rest, 
					{
					$cut_hash{$match}{$count}=$1;	#record domain boundaries
					$count++;
					}
				}
			if($score_hash{$score}==1)  #remember score, so appears only once in array
				{												#ie for multiple hits of same score
				push @{$score}, $match;	#@score keeps hits of that value
				next;
				}
			$score_hash{$score}=1;		#keep record that have seen that score
			@sorted=sort {$b<=>$a} @score_array;	#sort scores
			@score_array=@sorted;								 #sorted array
			if($score>=@score_array[9])	#add new score if bigger than last in list
				{
				push @{$score}, $match;
				@score_array[9]=$score;
				}
			}
		}
	&showdomssea; #outputs table of hits (together with associated links, cath codes etc)
	#&clearbar;		#bar thingy?
	}
###
#	output table of DomSSEA hits
###
sub showdomssea
	{
	#open structure database code file and read into hash
	open (DOMSSEATABLE, ">$OUT_DIR/$input_seq.domsseatable") or die "Can't create domssea table\n";
	print DOMSSEATABLE "Score\tMatch\tSSEA\tNo. Doms\tBoundaries\tSCOP code\n";
	
	open (DB_CODES, "$DOMSSEA_HOME/${STRUCTURE_DB}_CODES.dat");
	while ($line2 = <DB_CODES>)
		{
		if($line2=~/^(\S+)\s+(.*)/)
			{
			$db_code_hash{$1}=$2;	#record structure db (i.e. CATH/SCOP) codes
			}	
		}
	close(DB_CODES);

	@sorted=sort {$b<=>$a} @score_array;	#make sure score array is sorted 
	@score_array=@sorted;
###
#	start generating html
###	

	# Array is always of size 10, because entries
        # are added at position 9 ...
	# To eliminate non-existent entries ...
	# exclude ones with scores of 0 ... crude.
	for($hit=0;$hit < $#score_array;$hit++)
		{
		$x=@score_array[$hit];		#$x is score
		$length_array=@{$x};		#@{$x} conrtains matched codes

                if ($x == 0)
		    {
		    next;
		    }
                

###
#	only one match for this score, ie this score seen only once
###

		if($length_array==1)	
			{
			#no $k required on this bit so changed them...
			#$x is score
			#$numdom_hash{@{$x}[0]} is num doms
			#$cut_hash{@{$x}[$k]}{?} are cuts
			$match_num_doms=$numdom_hash{@{$x}[0]};
			print DOMSSEATABLE $x."\t";
			$url_code=substr(@{$x}[0],0,4);	#code to link to PDBSum
			print DOMSSEATABLE @{$x}[0]."\t";
	
			print DOMSSEATABLE $sec_align{@{$x}[0]}."\t";
			print DOMSSEATABLE $numdom_hash{@{$x}[0]}."\t";
			
			if($match_num_doms==1) #ie if match is a single domain chain, no boundaries
			{
				print DOMSSEATABLE "-";
			}
			for($p=1;$p<$match_num_doms;$p++)	#$p is not 0, starts at 1
			{

				if($use_cons!=1)
				{

					if(($cut_hash{@{$x}[0]}{$p}<40)||($cut_hash{@{$x}[0]}{$p}>($length-40)))
					{	
						print DOMSSEATABLE $cut_hash{@{$x}[0]}{$p};
					}
					else
					{
						print DOMSSEATABLE $cut_hash{@{$x}[0]}{$p};
					}
					if ($p < $match_num_doms - 1)
					{
					    print DOMSSEATABLE ", ";
					}
				}
				
###
#	method to use consensus boundary prediction
###
				#$use_con=1;
				if(($use_con==1)&&($num_doms==2))
					{
					#needs length and numdoms and DomSSEA cuts	
					$cut_DS=$cut_hash{@{$x}[0]}{$p};
					#	&consensus;
					#consensus returns $consensus_cut;
					#check cut region
					if($cons_stat==1) #no consensus met
						{
						print "No consensus met, using DomSSEA prediction";
						}
					if(($consensus_cut<40)||($consensus_cut>($length-40))) #doesn't matter if cons or not
						{	
						#	print HTOUT " ? ";
						}
					else
						{
						#	print HTOUT " $cut_hash{@{$x}[$k]}{$p} ";
						}
					}
				#here need to also output consensus hit..
				#print"$cut_hash{@{$x}[0]}{$p} kk"; #prints domain cut
				}	
			
                        # 15/06/2004 KB Added '0' domain suffix to scop code.
			@carray=();
			@carray=split(/,/,$db_code_hash{@{$x}[0].'0'});
			$pair_of_codes=0;
			#print NOUT "<TR><TD><FONT FACE='helvetica' SIZE='-1'>
			$num_cat_codes=scalar(@carray);
			print DOMSSEATABLE "\t";
			foreach $caco (@carray)
			{
				if($pair_of_codes==0)
				{
					#print NOUT "<TR><TD><FONT FACE='helvetica' SIZE='-1'>\n";
				}
				$pair_of_codes++;
				print DOMSSEATABLE $caco;
				if(($pair_of_codes==1)&&($num_cat_codes>1))
				{
					print DOMSSEATABLE ", ";
				}
				if($pair_of_codes==2)
				{
					#print NOUT "</FONT></TD></TR>\n";
					$pair_of_codes=0;
				}
				
			}
			print DOMSSEATABLE "\n";
		}
###
#	more than one match for this score, ie this score seen more than once
###
		else 
		{
			for($k=0;$k<=$length_array-1;$k++)
			{
				print DOMSSEATABLE $x."\t";
				$url_code=substr(@{$x}[$k],0,4);
				print DOMSSEATABLE @{$x}[$k]."\t";
				print DOMSSEATABLE $sec_align{@{$x}[$k]}."\t";
				print DOMSSEATABLE $numdom_hash{@{$x}[$k]}."\t";
				$match_num_doms=$numdom_hash{@{$x}[$k]};	#num doms of match

				if($match_num_doms==1)
				{
					print DOMSSEATABLE "-";
				}
				for($p=1;$p<$match_num_doms;$p++)
				{
					if($use_cons!=1)
					{
						if(($cut_hash{@{$x}[$k]}{$p}<40)||($cut_hash{@{$x}[$k]}{$p}>($length-40)))
						{	
						     print DOMSSEATABLE $cut_hash{@{$x}[0]}{$p}."";
						}
						else
						{
							print DOMSSEATABLE $cut_hash{@{$x}[$k]}{$p}."";
						}

						if ($p < $match_num_doms - 1)
						{
						    print DOMSSEATABLE ", ";
						}
				}
				print DOMSSEATABLE "\t";
###
#	method to use consensus boundary prediction
###
					#$use_con=1;	#if want to use, make == 1
					if(($use_con==1)&&($num_doms==2)) #only for two dom predictions
						{
						#needs length and numdoms and DomSSEA cuts	
						$cut_DS=$cut_hash{@{$x}[$k]}{$p};
						#	consensus;	#consensus returns $consensus_cut;
						#check cut region
						if($cons_stat==1) #no consensus met
							{
							print"No consensus met, using DomSSEA prediction";
							}
						if(($consensus_cut<40)||($consensus_cut>($length-40))) #doesn't matter if consensus or not
							{	
							#print HTOUT " ? ";
							}
						else
							{
							#print HTOUT " $cut_hash{@{$x}[$k]}{$p} ";
							}
						}
					}
				### outputting structure db codes
                                # 15/06/2004 KB Added '0' domain suffix to scop code.
				@carray=();	#
				@carray=split(/,/,$db_code_hash{@{$x}[$k].'0'});
				$pair_of_codes=0;
				#print NOUT "<TR><TD><FONT FACE='helvetica' SIZE='-1'>
				$num_cat_codes=scalar(@carray);
				foreach $caco (@carray)
				{
					if($pair_of_codes==0)
					{#print NOUT "<TR><TD><FONT FACE='helvetica' SIZE='-1'>\n";
					}
					
					$pair_of_codes++;
					print DOMSSEATABLE $caco;
					
					if(($pair_of_codes==1)&&($num_cat_codes>1))
					{
						print DOMSSEATABLE ", ";
					}
					
					if($pair_of_codes==2)
					{#print NOUT "</FONT></TD></TR>\n"; $pair_of_codes=0;
					}
					
				}
				}
		}
		@{$x}=();
		}	
	}	#end of sub showdomssea

###################################
# Consensus Subroutine			###SORT THIS AS NOT REALLY SORTED IT...
#	if input is consensus..  
#	for each multi dom hit, 
#	need to make consensus pred
###################################
sub consensus
	{
	$con_stat=0; #if 1, no consensus met
	&DGS;	#domain guess by size
	&DIL;	#difference in length
	&equal;	#equal division
	@cut=();
	@cuts=(0,0,$cut_DS,$cut_dgs,$cut_dil,$cut_equal);
	my %cnt;
	my $cons;
	for(my $i=2;$i<6;$i++)			#go thru cuts
  	{
  	for(my $j=$i;$j<6;$j++)		#pairwise
  		{
  		if($i==$j){next;}				#not against themselves
  		my $a=@cuts[$i];						
  		my $b=@cuts[$j];
  		my $diff=$a-$b; if($diff<0){$diff=$diff*(-1);}
  		if($diff<=10)
  			{
  			$cnt{$i}++;						#$i is array position
  			push @{$i}, $b;				#those grouped near this position
  			}
  		}
  	} 
	my $element=100;
	my $oldnum=0;
	foreach my $y (sort keys %cnt)		#number grouped together
		{
		if($cnt{$y}>$oldnum)				#find cut with most agreements
			{
			$element=$y;							#remember this array position
			$oldnum=$cnt{$y};
			}
		}
	undef %cnt;
	if($element==100)					#no consensus, use DomSSEA default
		{
		my $consensus_cut=@cuts[2];
		$con_stat=1;
		}
	if(($element>(1))&&($element!=100)) #ie 0,0,cuts....
		{
		$cons++;
		$consensus_cut=@cuts[$element];			#get cut (not in @{cut list}
		#print"$consensus_cut x @{$element}\n";
		my $length_s=@{$element};						#number grouped near
		for(my $s=0;$s<$length_s;$s++)
			{
			$consensus_cut=$consensus_cut+${$element}[$s];
		#	print"${$element}[$s] $length_s f\n";
			}
		$consensus_cut=int ($consensus_cut/($length_s+1))
		}
	@{0}=@{1}=@{2}=@{3}=@{4}=@{5}=@{6}=@{7}=@{8}=();

###
# DGS Subroutine											
###
	sub DGS
		{
  	my %pdom;
		my %pnd;
		my %twodom;
		my %threedom;
		my %scorerec;
		my %seen;
		my $step;
		my $pnds_length;
		my $c;
		my $length5;
		my $length40;
		my $score;
 	 open(PROBS,"Probabilities") or die "Can't open Probabilities\n";;
 	 while($data = <PROBS>)
  	{
  	if($data =~ /^(\d+)\s+(\d+\.\d+)/)
  		{
  		for(my $i=$1;$i<($1+40);$i++)
  			{
  			$pdom{$i}=$2;
  			}
  		}
  	if($data =~ /^(\d+)\s+(\d+)\s+(\d+\.\d+)/)
  		{
  		$pnd{$1}{$2}=$3;
  		}
  	}
	my $length5=(int ($length/5))*5;				#times 5
	if((($length/5)-(int ($length/5)))>=0.5){$length5+=5;}#round up end cut length
	for($i=0;$i<=600;$i+=20)	
		{
		if(($length<($i+20))&&($length>=$i))
			{
			$length40=$i;
			if((($length/40)/(int $length/40))>1){ $pnds_length=(int $length/40)*40;}#p(no. doms |c length)
			else {$pnds_length=$length}
			if($length40<200){$step=20;}
			if(($length40>=200)&&($length40<400)){$step=20;}
			if($length40>=400){$step=25;}
			last;
			}
		}
	if($num_doms==2)		#DomSSEA says two dom
		{	
		my $ptot=0;
		@twodom=();
		$c=$pnd{2}{$pnds_length};															#prob of 2 doms for this length
		for ($i=40;$i<=($length40-40);$i=$i+$step)						#calc length,
			{
			push @twodom, $i;																		#add first dom length
			push @twodom, ($length5-$i);												#add second dom length
			push @twodom, (($pdom{$i})*($pdom{($length5-$i)}));	#add to array		
			$ptot+=(($pdom{$i})*($pdom{($length5-$i)}));				#total prob
			}
		for($i=0;$i<=(scalar@twodom-3);$i+=3)		#itterate through and calc scores
			{
			$score=($c*(@twodom[$i+2]/$ptot));		#overall score
			if($score==0){next;}
			if($score>0){$score=log $score;}			#log it
			if($seen{$score}==1){$score=$score+0.000001;}
			$scorerec{$score}="2 @twodom[$i] @twodom[$i+1]";	#save score
			$seen{$score}=1;
			}	
		undef %seen;
		}
	if($num_doms==3)		#DomSSEA says three+ dom
		{
		$ptot=0;
		@threedom=();
		$c=$pnd{3}{$pnds_length};
		for ($i=40;$i<=($length40-40);$i+=$step)
			{
			for (my $j=40;$j<=($length40);$j+=$step)
				{
				if(($length5-$i-$j)<40) {next};
				push @threedom, $i;
				push @threedom, $j;
				push @threedom, ($length5-$i-$j);
				push @threedom, ($pdom{$i}*$pdom{$j}*$pdom{($length5-$i-$j)});
				$ptot+=($pdom{$i}*$pdom{$j}*$pdom{($length5-$i-$j)});
				}
			}
		for($i=0;$i<=(scalar@threedom-4);$i+=4)
			{
			$score=($c*(@threedom[$i+3]/$ptot));
			if($score==0){next;}
			if($score>0){$score=log $score;}
			if($seen{$score}==1){$score=$score+0.000001;}
			$scorerec{$score}="3 @threedom[$i] @threedom[$i+1] @threedom[$i+2]";
			$seen{$score}=1;
			}
		undef %seen;
		}
###	
# output stuff to consensus subroutine
###	
	$top_ten=0;
	foreach my $poss (sort {$b<=>$a} keys %scorerec)
		{
		if($scorerec{$poss}=~/^(\d+)\s+\d+/)
			{
			my $domain_num=$1;
			if($domain_num>1) #which it should be as not predicting single
				{
				if($scorerec{$poss}=~/^\d+\s+(\d+)/)
					{
					$cut_dgs=$1;
					#print"$cut xxx\n";
					}
				if($scorerec{$poss}=~/^\d+\s+(\d+)\s+(\d+)/)
					{
					$dgs_hash{$el_A}{1}=$1;
					$dgs_hash{$el_A}{2}=$2;
					}#a global dgs_cut hash?				
				}
			}
			last; #got first prediction
			}
		}
###
# Difference in length Subroutine					 
###

	sub DIL
		{
                # KB - ONLY USED FOR TESTING - DO NOT THINK WE NEED TO CHANGE FOR DATABASE TYPE (CATH/SCOP).

		#Needs database of lengths and domain number and cuts (CATH definitions)
		my %len_hash;
		my %len_cut;
		#read in database (here 2struc_H)
		open (DOM, "$DOMSSEA_HOME/length_data") or die "Can't open domlist\n";
		while ($line = <DOM>)
			{
			if($line=~/^(\S+)\s+(\d+)\s+2\s+(\d+)/)  #THIS IS ONLY GETTING TWO DOMAINS!
				{
				$code=$1;
				push @list, $code;
				$len_hash{$code}=$2;
				$len_cut{$code}=$3;
				}
			}
		#need to separate library into 2 and 3 domain chains
		#no single needed..
		#at moment only two domain
		my %score_hash;
		my %match;
		for(my $i=0;$i<(scalar@list);$i++)				#iterate over library
			{
			my $lib_seq=@list[$i];									#get library sequence
			my $lib_length=$len_hash{$lib_seq};			#length lib sequence
			if($lib_length>$length)									#lib seq > target seq
				{
				my $score=1-(($lib_length-$length)/$lib_length);  #difference in length
				}
			if($length>$lib_length)
				{
				$score=1-(($length-$lib_length)/$length);
				}
			#if($score==$score_hash{$lib_seq})		#if score is equal to one seen before 
			#	{																		#take last one for  moment
			#	push @{$lib_seq},$el_B;							#add on to list of hits
			#	}
			if($score>=$score_hash{$target})			#if score is bigger than one seen before
				{
				#$nums{$lib_seq}=1;
				#@{$lib_seq}=();										#wipe list hits
				#push @{$lib_seq},$el_B;						#add new one
				$score_hash{$target}=$score;				#remember new top score
				$match{$target}=$lib_seq;						#what matched to (to get cuts)
				}
			}
		#so hit at moment is taken as
		$hit=$match{$target}; #$hit is $lib_seq ie whats matched so now get cut
		$cut_dil=$len_cut{$hit};
		#$cut_dil=170;
		#get cut data somehow from library.
		}
###
# equal division				
###
	sub equal
		{
		my $fragment=int ($length/$num_doms);
		for(my $i=1;$i<=$num_doms;$i++)
			{
			$j_plus++;
			$tot_over+=$fragment;
			$equal_hash{$lib_seq}{$j_plus}=$tot_over;
			}
		$cut_equal=167;
		}
	}#end of consensus subroutine

###
#	PSIPRED output control
###
$ls=$input_seq ."_ls";
$ls_seq=$input_seq; #to get the number of psipred figures
if($do_ppred eq "ppyes")
	{
	system ("ls $OUT_DIR/$ls_seq*.png > $OUT_DIR/$ls.dat");
	open (LS, "$OUT_DIR/$ls.dat") or die "can't open ls\n";
	while ($type=<LS>)
		{
		if($type=~/thumb|myps|tran/)
			{
			}
		else
			{
			$ls_count++;
			}
		}
	system ("rm $OUT_DIR/$ls.dat");
	for($jp=1;$jp<=$ls_count;$jp++)
		{
		$file_amalg=$input_seq . "_" . "$jp";
		$file_amalg2=$file_amalg. "thumb";
		system("cp $OUT_DIR/$file_amalg.png $OUT_DIR/$file_amalg2.png");
		system("mogrify -geometry 500x800! -crop 0x0 -quality 100 $OUT_DIR/$file_amalg2.png");
		}#-geometry 300x400#-crop 18x5%
	#&psipredbit; #output to html file
	}
#&endpage;	#output end of results page
#&rmrubbish;	#remove unwanted files
sub rmrubbish
 {
 system("rm $OUT_DIR/$input_seq.blast");
 system("rm $OUT_DIR/$input_seq.domssea");
 system("rm $OUT_DIR/$input_seq.horiz");
 system("rm $OUT_DIR/$input_seq.parseDS");
 system("rm $OUT_DIR/$input_seq.pred");
 system("rm $OUT_DIR/$input_seq.ss");
 system("rm $OUT_DIR/$input_seq.ss2");
 system("rm $OUT_DIR/$input_seq.gnu");
 system("rm $OUT_DIR/$input_seq.graph");
 system("rm $OUT_DIR/$input_seq.blastdom");
 system("rm $OUT_DIR/$myps.ps");
 system("rm $OUT_DIR/$myps.png");
 system("rm $OUT_DIR/magic*");
 }
###make files writeable executable by me
system("chmod a+rwx $OUT_DIR/$input_seq*");
 
###
#	All the html output subroutines
### 
sub printtoparsed
	{
	#print SEE "<PRE>\n";
	print SEE "Hit name	";
	print SEE "Query length	";
	print SEE "Query start	";
	print SEE "Query end	";
	print SEE "Hit length	";
	print SEE "Hit start	";
	print SEE "Hit end	";
	print SEE "E-value	";
	print SEE "\n";
	}


sub graphpage
{
	if ($peak_detected>0)
	{
   		print BOUNDARYOUT "Putative domain boundaries located in PSI-BLAST alignment profile:\n";
	}
	print BOUNDARYOUT "Number of predicted domains by DPS: ".($num_peaks+1)."\n";
	print BOUNDARYOUT "Domain Bounary locations predicted DPS: ";
	foreach $put_link (sort {$a<=>$b} @termini_preds)
	{
		print BOUNDARYOUT $put_link." ";
	}
}
