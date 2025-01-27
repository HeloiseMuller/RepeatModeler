#!/usr/bin/perl

use Getopt::Long;
use Cwd qw(abs_path getcwd cwd);

#Arguments
my @getopt_args = ('-help',
                    'i=s',
                    '-RRDir=s',
                    '-LTRtmp=s',
                    '-database=s',
                    '-o=s',
                    'ninja_dir=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

my $genome = $options{'i'};
my $dirRR =  abs_path($options{'RRDir'});
my $dirLTR = abs_path($options{'LTRtmp'});
my $genomeDB = $options{'database'};
my $out = $options{'o'};
my $ninjaDir =  $options{'ninja_dir'};

#Create directory for tmp files
system( "mkdir -p $dirLTR" );

#Run LTRPipeline
system( "LTRPipeline -debug $genome -ninja_dir $ninjaDir -tmpdir $dirLTR" );

#Continue only if LTRPipeline found something
$check = qx(ls $dirLTR/ | grep LTR_ );
chomp($check);

if ( -s "$dirLTR/${check}/families.fa" ) {

# 1. Combine results from both pipelines into a file
system(
"cat $dirLTR/LTR*/families.fa $dirRR/consensi.fa > $dirLTR/combined.fa" );
system(
"cat $dirLTR/LTR*/families.stk $dirRR/families.stk > $dirLTR/combined.stk"
      );
      
# 2. Cluster results
my $cmd =
 "cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 "
 . "-i $dirLTR/combined.fa -o $dirLTR/cd-hit-out"
 . ">$dirLTR/cd-hit-stdout 2>&1";
system( $cmd);
      
# 3. Process clusters and remove redundancy
my %redundant_families;
my %putative_subfamilies;
my $ltrFamCnt = 0;
my $rrFamCnt  = 0;
if ( -s "$dirLTR/cd-hit-out.clstr" ) {
  open IN, "<$dirLTR/cd-hit-out.clstr"
    or die
"RepeatModeler: Could not open cd-hit-out.clstr for reading!\n";
  my @cluster = ();
  my $longest_ltr_id;
  my $longest_ltr_size;
  my $longest_rnd_id;
  my $longest_rnd_size;
  while ( <IN> ) {
    if ( /^>Cluster/ ) {
      if ( @cluster > 1 ) {

              # I started out by picking the longest LTRPipeline
              # derived family as the cluster rep.  Now we
              # consider all LTRPipeline candidates as dominant.
              # The only reason two more more LTRPipeline candidates
              # would show up in a cluster would be due to
              # overlap between LTR/INT sequences that didn't get
              # removed by LTR_retriever.
              #
              # Keeping the longest_ltr_id variable as an indicator
              # that there is at least one LTRPipeline candidate in the
              # cluster.
        if ( $longest_ltr_id ) {
          foreach my $id ( @cluster ) {

            # Now...remove only overlapping RECON/RepeatScout
            # candidates ( rnd-#_family-# ).
            if ( $id !~ /^ltr-\d+_family-\d+/ ) {

              # Use the longest as the "reason" why we are removing it.
              #print "Removing $id because it's redundant with $longest_ltr_id\n";
              $redundant_families{$id}++;
             }
           }
         }
         elsif ( $longest_rnd_id ) {
          foreach my $id ( @cluster ) {
            if ( $id ne $longest_rnd_id ) {

              #print "Labeling $id as putative subfamily of $longest_rnd_id\n";
              $putative_subfamilies{$id} = $longest_rnd_id;
             }
           }
         }
      }
      $longest_ltr_id   = undef;
      $longest_ltr_size = 0;
      $longest_rnd_id   = undef;
      $longest_rnd_size = 0;
      @cluster          = ();
      next;
    }
    if ( /^\d+\s+(\d+)nt,\s*>((rnd|ltr)-\d+_family-\d+)/ ) {
      my $size = $1; #size of the element in nt
      my $id   = $2; #eg: ltr-1_family-25
      my $type = $3; #ltr or rnd
      #those 3 variables for both rnd and ltr
      if ( $type eq "ltr" ) {
        $ltrFamCnt++;
        if ( $size > $longest_ltr_size ) {
          $longest_ltr_id   = $id;
          $longest_ltr_size = $size;
        }
      }
      if ( $type eq "rnd" ) {
              $rrFamCnt++;
              if ( $size > $longest_rnd_size ) {
                $longest_rnd_id   = $id;
                $longest_rnd_size = $size;
              }
            }
            push @cluster, $id;
          }
        }
        close IN;

        # TODO: Do we really want to keep this?
        #unlink("cd-hit-out.clstr" );
      }
      print "\n\n$rrFamCnt RepeatScout/RECON families\n";
      print "$ltrFamCnt LTRPipeline families\n"; 
      if ( keys( %redundant_families ) ) {
        system(
               "cp $dirRR/consensi.fa $dirLTR/consensi.fa.recon_rscout_only" );
        system( "mv $dirLTR/combined.fa $dirLTR/consensi.fa.with_redundancy" );

        # Filter consensi.fa and families.stk
        open IN, "$dirLTR/consensi.fa.with_redundancy"
            or die
"RepeatModeler: Could not open consensi.fa.with_redundancy for reading";
        open OUT, ">$dirLTR/consensi.fa"
            or die
            "RepeatModeler: Could not open consensi.fa for writing";
        my $id;
        my $data;
        while ( <IN> ) {
          if ( /^>(\S+)/ ) {
            my $tmpID = $1;
            if ( $data ) {
              if ( !exists $redundant_families{$id} ) {
                print OUT $data;
              }
            }
            if ( exists $putative_subfamilies{$tmpID} ) {
              my $tstr = $_;
              $tstr =~ s/[\n\r]//g;
              $data =
                  "$tstr [ putative subfamily of "
                  . $putative_subfamilies{$tmpID} . " ]\n";
            }
            else {
              $data = $_;
            }
            $id = $tmpID;
            next;
          }
          $data .= $_;
        }
        if ( $data ) {
          if ( !exists $redundant_families{$id} ) {
            print OUT $data;
          }
        }
        close IN;
        close OUT;
        system(
             "cp $dirRR/families.stk $dirLTR/families.stk.recon_rscout_only" );
        system(
               "mv $dirLTR/combined.stk $dirLTR/families.stk.with_redundancy" );
        open IN, "<$dirLTR/families.stk.with_redundancy"
            or die
"RepeatModeler: Could not open families.stk.with_redundancy for reading";
        open OUT, ">$dirLTR/families.stk"
            or die
            "RepeatModeler: Could not open families.stk for writing";
        $id   = "";
        $data = "";

        while ( <IN> ) {
          if ( /^#=GF\s+ID\s+(\S+)/ ) {
            $id = $1;
          }

          if ( /^#=GF\s+DE\s+(\S.*)/ ) {
            if ( exists $putative_subfamilies{$id} ) {
              my $tstr = $_;
              $tstr =~ s/[\n\r]//g;
              $data .=
                  "$tstr [ putative subfamily of "
                  . $putative_subfamilies{$id} . " ]\n";
            }
            else {
              $data .= $_;
            }
          }
          else {
            $data .= $_;
          }

          if ( /^\/\// ) {
            if ( !exists $redundant_families{$id} ) {
              print OUT $data;
            }
            $data = "";
          }
        }
        close IN;
        close OUT;

        my $removed = scalar( keys( %redundant_families ) );
        print "Removed $removed redundant LTR families.\n";
        my $count = ( $rrFamCnt + $ltrFamCnt ) - scalar( keys( %redundant_families ) );
        print "Final family count = $count\n"
      }

  else {
    #if no redundency, save the combined.fa 
    system( "mv $dirLTR/combined.fa $dirLTR/consensi.fa" );
    system ( "mv $dirLTR/combined.stk $dirLTR/families.stk" );
    my $count = ( $rrFamCnt + $ltrFamCnt );
    print "No redundant LTR families.\n";
    print "Final family count = $count\n";
  }
  
print "Working directory:  $dirLTR\n";
print "may be deleted unless there were problems with the run.\n";

system("RepeatClassifier "
  . "-consensi $dirLTR/consensi.fa -stockholm $dirLTR/families.stk -pa 10" );
system( "cp $dirLTR/consensi.fa.classified $out/$genomeDB-families.fa" )
  if ( -s "$dirLTR/consensi.fa.classified" );
system( "cp $dirLTR/families-classified.stk $out/$genomeDB-families.stk" )
  if ( -s "$dirLTR/families-classified.stk" );
  
  }
  
  #if no output with LTRPipeline, simply take RR as final output
  else { 
  system("RepeatClassifier "
  . "-consensi $dirRR/consensi.fa -stockholm $dirRR/families.stk -pa 10" );
  system( "cp $dirRR/consensi.fa.classified $out/$genomeDB-families.fa" )
  if ( -s "$dirRR/consensi.fa.classified" );
system( "cp $dirRR/families-classified.stk $out/$genomeDB-families.stk" )
  if ( -s "$dirRR/families-classified.stk" );
  }
  
print "\nThe results have been saved to:\n";
print  "$out/$genomeDB-families.fa  - Consensus sequences for each family identified.\n";
print "$out/$genomeDB-families.stk - Seed alignments for each family identified.\n";
