#!/usr/bin/perl

#I need this funciton:
sub log_print {
 my $string = shift;
 print "$string";
 print $LOG "$string";
}

# 3. Process clusters and remove redundancy
my %redundant_families;
my %putative_subfamilies;
my $ltrFamCnt = 0;
my $rrFamCnt  = 0;
if ( -s "cd-hit-out.clstr" ) {
  open IN, "<cd-hit-out.clstr"
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
      print $rrFamCnt\n;
      print $ltrFamCnt\n; #THose two print work, so why log_print does not work of ltreFamCnt
      log_print "       - $rrFamCnt RepeatScout/RECON families\n";
      log_print "       - $ltrFamCnt LTRPipeline families\n";
      if ( keys( %redundant_families ) ) {
        system(
               "mv consensi.fa consensi.fa.recon_rscout_only" );
        system( "mv combined.fa consensi.fa.with_redundancy" );

        # Filter consensi.fa and families.stk
        open IN, "consensi.fa.with_redundancy"
            or die
"RepeatModeler: Could not open consensi.fa.with_redundancy for reading";
        open OUT, ">consensi.fa"
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
             "mv families.stk families.stk.recon_rscout_only" );
        system(
               "mv combined.stk families.stk.with_redundancy" );
        open IN, "<families.stk.with_redundancy"
            or die
"RepeatModeler: Could not open families.stk.with_redundancy for reading";
        open OUT, ">families.stk"
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
        #### 1st log_print work but not the second
        log_print "       - Removed "
            . scalar( keys( %redundant_families ) )
            . " redundant LTR families.\n";
        #log_print "       - Final family count = "
        #    . (
        #    ( $rrFamCnt + $ltrFamCnt ) - scalar( keys( %redundant_families ) ) )
        #    . "\n";
      }

  else {
    log_print
"\nWARNING: Could not create input file for LTRPipeline from $genomeDB! Continuing using\nresults from RepeatScout/RECON pipeline only.\n";
  }
  log_print "LTRPipeline Time: " . elapsedTime( 1 ) . "\n";
