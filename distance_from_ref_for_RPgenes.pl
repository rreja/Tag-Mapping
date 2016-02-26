use strict;
use warnings;
use Getopt::Std;
use File::Basename;


#reading all the commans line options
my %opt;
getopt('irud',\%opt);

&help_message if (defined $opt{h});

####### Variable declaration ######
my $up_dist = 1000; my $down_dist = 1000; my $sort = 1;
if(defined $opt{s}) {$sort = 0;}

# reading the directory that has the index files in  it
my $idxdir = $opt{'i'};
$idxdir = check_dir($idxdir);
# reading the directory that has the reference files in  it
my $refdir = $opt{'r'};
$refdir = check_dir($refdir);

if (defined $opt{u}){
    $up_dist = $opt{'u'}
}

if (defined $opt{d}){
    $down_dist = $opt{'d'}
}



opendir IDX, $idxdir || die "Cannot open the directory";

unless(-d $idxdir."output_all/"){
    system("mkdir ".$idxdir."output_all/");
    }

while((my $filename = readdir(IDX))){
    # checking for .idx / .tab file in the folder
    my $gff;
    next if(check_extension($filename) < 1 || check_extension($filename) > 3);
    if((check_extension($filename) == 2) || (check_extension($filename) == 3)){
        $gff = 0
    }
    elsif(check_extension($filename) == 1){
        $gff = 1
    }
    my $base_name = return_basename($filename);
    
    #creating  a tmp folder to store tmp files
    unless(-d $idxdir."tmp/"){
    system("mkdir ".$idxdir."tmp/");
    }
    
    print "\nINFO: Index file detected as input! Abort if incorrect\n";
    print "\nINFO: Input tag/index file = ".$filename;
    open IN, $idxdir.$filename || die "File not found";
    # open up the index/tag file and divide it into fwd and rev tag GFF files
    print "\nINFO: Splitting the tags file into fwd and rev gff files\n";
    open OUT1,">".$idxdir."tmp/fwd.tmp" || die "File not found";
    open OUT2,">".$idxdir."tmp/rev.tmp" || die "File not found";
    if($gff == 1){
        while(<IN>){
            chomp($_);
            next if(($_ =~ /^#/) || ($_ =~ /^chrom/));
            my @cols = split(/\t/,$_);
            if($cols[6] eq "+"){
                print OUT1 $_."\n";
            }
            elsif($cols[6] eq "-"){
                print OUT2 $_."\n";
            }
            
        }
        
    }
    else{
        while(<IN>){
        chomp($_);
        next if(($_ =~ /^#/) || ($_ =~ /^chrom/));
        my @cols = split(/\t/,$_);
        if($cols[2] > 0 && $cols[1] > 0){
            print OUT1 $cols[0]."\t.\t.\t".$cols[1]."\t".($cols[1]+1)."\t".$cols[2]."\t+\t.\t.\n";
        }
        if($cols[3] > 0 && $cols[1] > 0){
            
            print OUT2 $cols[0]."\t.\t.\t".$cols[1]."\t".($cols[1]+1)."\t".$cols[3]."\t-\t.\t.\n";
        }
     }
    }

    close(IN);close(OUT1);close(OUT2);
    # read the index file. Moving on to the refernece directory.
    opendir REF, $refdir || die "Cannot open the directory";
    while( (my $refname = readdir(REF))){ # Splitting each reference file into sense and anti-sense files.
        next if(check_extension($refname) != 1);
        my $ref_basename = return_basename($refname);
        print "INFO: Reference file = ".$refname."\nINFO: Splitting it into sense and antisesne files\n";
        open IN, $refdir.$refname || die "File not found";
        open OUT1,">".$idxdir."tmp/sense.tmp" || die "File not found";
        open OUT2,">".$idxdir."tmp/antisense.tmp" || die "File not found"; 
        while(<IN>){
            chomp($_);
            next if(($_ =~ /^#/) || ($_ =~ /^chrom/));
            my @cols = split(/\t/,$_);
            # First create sesnse and antisense ref feature file and add/subtract the up/down-stream distance from the start co-ordinate.
            if($cols[6] eq "+"){
                if($cols[3] - $up_dist <= 0){
                    print OUT1 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t1\t".($cols[3]+$down_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3].":".$cols[4]."\t".$cols[8]."\n";
                }
                else{
                    print OUT1 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t".($cols[3]-$up_dist)."\t".($cols[3]+$down_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3].":".$cols[4]."\t".$cols[8]."\n";
                }
            
            }
            else{
                if($cols[4] - $down_dist <= 0){
                    print OUT2 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t1\t".($cols[4]+$up_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3].":".$cols[4]."\t".$cols[8]."\n";
                }
                else{
                    print OUT2 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t".($cols[4]-$down_dist)."\t".($cols[4]+$up_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3].":".$cols[4]."\t".$cols[8]."\n";
                }
            
            }
         # completed splitting the reference into sense and anti-sense gff files   
        }
        close(IN);close(OUT1);close(OUT2);
    # Running intersectBed on the fwd/sense, rev/sense, fwd/antisense and rev/antisense and then joining fwd/sense with rev/antisense and rev/sense with fwd/antisense.
    # Feature with no tags would have "0" as overlap and features with tags would have >0 as overlap.
    #print "/usr/local/bin/closestBed -a ".$idxdir."fwd.gff -b ".$refdir."sense.gff  >".$idxdir."output/fwd_sense.txt\n";
    print "INFO: Running intersectBed on the generated files\n";
    system("/usr/local/bin/intersectBed -wao  -a ".$idxdir."tmp/sense.tmp -b ".$idxdir."tmp/fwd.tmp  >".$idxdir."tmp/fwd_sense.tmp");
    system("/usr/local/bin/intersectBed -wao  -a ".$idxdir."tmp/sense.tmp -b ".$idxdir."tmp/rev.tmp  >".$idxdir."tmp/rev_sense.tmp");
    system("/usr/local/bin/intersectBed -wao  -a ".$idxdir."tmp/antisense.tmp -b ".$idxdir."tmp/fwd.tmp  >".$idxdir."tmp/fwd_antisense.tmp");
    system("/usr/local/bin/intersectBed -wao  -a ".$idxdir."tmp/antisense.tmp -b ".$idxdir."tmp/rev.tmp  >".$idxdir."tmp/rev_antisense.tmp");
    # Concatenate fwd/sense with rev/antisense and rev/sense with fwd/antisense
    print "INFO: Combining fwd/sense with rev/antisense and rev/sense with fwd/antisense\n";
    system("cat ".$idxdir."tmp/fwd_sense.tmp ".$idxdir."tmp/rev_antisense.tmp >".$idxdir."tmp/sense_distance.tab");
    system("cat ".$idxdir."tmp/rev_sense.tmp ".$idxdir."tmp/fwd_antisense.tmp >".$idxdir."tmp/antisense_distance.tab");
    # remove all the tmp file in output direcotry.
    #system("rm ".$idxdir."tmp/*.tmp");
    #system("rm ".$refdir."*.tmp"); #uncomment this to delete the tmp files generated in idx directory.
    # opening sense_distance.txt and antisesnse_distance.txt to calculate distance
    print "INFO: Finding the distance from reference feature\n";
    opendir COM, $idxdir."tmp/" || die "Cannot open the directory";
    my @order;
    while( (my $comname = readdir(COM))){
        next if(check_extension($comname) != 3);
        open IN1, $idxdir."tmp/".$comname || die "File not found";
	
        #open OUT1, ">".$idxdir."output_".$base_name."/".$ref_basename."_".return_basename($comname).".txt" || die "File not found";
	open OUT1, ">".$idxdir."output_all/".$base_name."_".return_basename($comname).".txt" || die "File not found";
	
	
	
        # create CDT folder.        
        unless(-d $idxdir."output_all/CDT/"){
              system("mkdir ".$idxdir."output_all/CDT");
            }
        open OUT2, ">".$idxdir."output_all/CDT/".$base_name."_".return_basename($comname).".cdt" || die "File not found";
        
        my %hash = (),my $k =0;
        if($comname eq "sense_distance.tab"){$k = 1;} # Make k as 1 when reading a sense file.
        while(<IN1>){
            chomp($_);
            my $dist; my $key;
            my @cols = split(/\t/,$_);
	    #if($cols[8] eq "."){
		#### Uncomment when you want key to be chr3:1234-2345
		#$key = $cols[0].":".$cols[7].":".$cols[6]; # concatenating the interval co-ordinate in the following format: chr3:1234-2345
		#### Uncomment below when you want to match the tag and motif together. Change the distance accordingly.
		####my @rtmp = split(/-/,$cols[7]);
		#####my $st = $rtmp[0] - 101;
		####my $en = $rtmp[1] + 100;
		####$key = $cols[0]."_".$st."_".$en."_".$cols[6];
		
	    #}
#            else{
		$key = $cols[8]; # comment this out if the 9th column of the reference feature is "."
#	    }
	    my $value;
            if($cols[12] == -1){ # If for a reference feature no tag exsists then make the distance as 0 and write the entry as it is.
		#print "hello".$cols[12];
                $dist = 0;
                my $value = "*:*";
                push(@{$hash{$key}},$value);
                print OUT1 $_."\t".$dist."\n"; 
            }
            else{
            my @start = split(/:/,$cols[7]);    
            $dist = $cols[12] - $start[0];
            if($cols[6] eq "+"){  # The distance be - if the input feature is smaller then + strand REF feature.
                $value = $dist.":".$cols[14];  # concatenating the co-ordinate and the tag number at the coordinate: ex., 1234:10
                push(@{$hash{$key}},$value);
                print OUT1 $_."\t".$dist."\n"; 
            }
            else{                 # The distance be - if the input feature is larger than - strand REF feature.
                $dist = (-1)*($cols[12] - $start[1]); # Flipping step.
                print OUT1 $_."\t".$dist."\n";
                $value = $dist.":".$cols[14];  # concatenating the co-ordinate and the tag number at the coordinate: ex., 1234:10
                push(@{$hash{$key}},$value);
                print OUT1 $_."\t".$dist."\n"; 
            }
             
            } # end of the else #1.
            
             
        }   # end of nested while .
        print "INFO: Generating CDT file in the CDT folder for  ".$comname."\n";
        my %hash1 = create_cdt(\%hash,\*OUT2,$k,$sort);
        print_header(\*OUT2);
        if($k == 1 && $sort == 1){
            print_hash_by_order(\%hash1,\@order,\*OUT2);
        }
        else{
            @order = print_hash(\%hash1,\*OUT2);
        }
        
        close(IN1);close(OUT1); close(OUT2);
        
        
    } # end of while

    print "INFO: Completed calculating distance for ".$filename." and ".$refname."\n\n";
    
  } # end of the loop over reference directory.
    print "INFO: Your output is present in ".$idxdir."output_all/".$base_name."/\n";
    system("rm -fr ".$idxdir."tmp/");
      
}
 


###################### Subroutine declarations ###################

sub create_cdt{
    my ($hash,$OUT2,$k,$sort) = @_;
    my %hash1 = ();
    foreach my $val (sort keys %$hash)
{
	
	my %tags; #hash of co-ordinate and the tag count.
	#my @array_to_print;
        my $string_to_print = ""; # remove the one from here in case you decide to print the sun total in the second column.
	my $total = 0;
	
	foreach my $tag_info (@{$$hash{$val}})
	{
		
		my @tmp3 = split(/:/,$tag_info);
		$tags{$tmp3[0]} = $tmp3[1];
		
	}
	for(my $i= (-1)*$up_dist;$i<=$down_dist;$i++)
	{
                
		if(exists($tags{$i}))
		{
			#push(@array_to_print,$tags{$i});
                        $string_to_print = $string_to_print."\t".$tags{$i};
			$total = $total+$tags{$i};
			
		}
		else
		{
			#push(@array_to_print,0);
                        $string_to_print = $string_to_print."\t0";
			
		}
               
	}
        if($k == 1 && $sort == 1){
            $hash1{$val} = $total.$string_to_print;   # Create a different hash value when the file is sense distance since no sorting is needed. Order is predetermined.
            #$hash1{$val} = $string_to_print;
        }
        else{
            $hash1{$val."\t".$total.$string_to_print} = $total; # Create the hash in a way that enables sorting.
        }
        
	#print $OUT2 $val."\t".$total;
	#print_array(\@array_to_print, $OUT2);
        
}
    return(%hash1);
    
}



sub print_hash{
    
    my($hash1,$OUT2) = @_;
    my @order;
    foreach my $value (sort {$$hash1{$b} <=> $$hash1{$a} }keys %$hash1)
    {
     print $OUT2 $value."\n";
     my @tmp = split(/\t/,$value);
     push(@order,shift(@tmp));
     #print $OUT2 $$hash1{$value}."\t".$value."\n";
    }
    return(@order);
}

sub print_hash_by_order{
    my($hash,$order,$OUT2) = @_;
    foreach my $id (@$order){
        print $OUT2 $id."\t".$$hash{$id}."\n";
        
    }
    
}
sub print_array{
    
	my ($array_to_print,$OUT2) = @_;
	foreach my $val (@$array_to_print)
	{
		print $OUT2 "\t".$val;
	}
	print  $OUT2 "\n";
}


sub print_header{
    my $OUT2 = shift;
    
    print $OUT2 "Uniqe ID\tGWEIGHT\t";
    for(my $i = (-1)*$up_dist;$i<=$down_dist;$i++){
        print $OUT2 $i."\t";
    }
    print $OUT2 "\nEWEIGHT\t\t";
    
    for(my $i = (-1)*$up_dist;$i<=$down_dist;$i++){
        print $OUT2 "1\t";
    }
    print $OUT2 "\n";
    
}


sub check_dir{
    my $dir = shift;
    my $tmp = substr($dir,-1,1);
    if($tmp eq "/"){
        return($dir);
    }
    else{
        return($dir."/");
    }
}


sub return_basename{
    my $filename = shift;
    my @suffixes = (".gff",".idx",".tab",".txt",".tmp",".tab");
    my $basename = basename($filename, @suffixes);
    return $basename;
    
}
sub check_extension{
    my $filename = shift;
    my @suffixes = (".gff",".idx",".tab",".txt");
    my ($fname,$path,$suffix) = fileparse($filename,@suffixes);
    if($suffix eq ".gff"){
        return 1;
    }
    elsif($suffix eq ".idx"){
        return 2;
    }
    elsif($suffix eq ".tab"){
        return 3;
    }
    elsif($suffix eq ".txt"){
        return 4;
    }
    else{
        return 0;
    }
}

sub help_message {
  print qq{
    Program: distance_from_ref.pl (Calculate distance from reference feature and generate CDT files)
    Contact: Rohit Reja <rzr142\@psu.edu>
    Usage:   distance_from_ref.pl -i <path_to_index_file> -r <path_to_REF_feature_files>

    Options: -i <path1>     path to the folder with index files in .idx or .tab format 
             -r <path2>     path to the folder with reference feature files in gff format.
             -u <number>    Upstream distance to go, defualt = 1000bp.
             -d <number>    Downstream distance to go, default = 1000bp.
             -s             If you want to arrange the sense CDT in the order of antisense CDT (Defualt). In case "-s" is used
                            both the CDT files will be sorted from high to low independent of each other.
             -h             help

    Example:
      perl distance_from_ref.pl -i  /usr/local/folder -r /folder1/ref/
      
    Output:
    Produces a sense and antisense distance file for each index file and then produces a sense and antisense
    CDT file for the same index file.
  
  };
  exit;
}







