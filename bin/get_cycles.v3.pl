#################################################################################
###### 08-25-2020                                                          ######
###### Author: Georgette Femerling                                         ######
######                                                                     ######
###### Searches for cycles in the metabolism in order to find new functions######
###### Argument usage:                                                     ######
###### [0] - Input TF                                                      ######
###### [1] - Path to Regulatory Network. Genes separated by //             ######
###### [2] - Path to gu_library directory with gene_links and GenProtEC    ######
###### Output:                                                             ######
######       One Cycle report per Multifun term                            ######
######        (metabolism or cell processes) associated with input TF      ######
###### History:                                                            ######
######         v1 - takes gene regulators from Ecocyc database             ######
######         v2 - takes gene regulators from input network               ######
######            - changes output to 3 fields                             ######
######            - adds number of genes in each MF category in output     ######
######            - marks genes in rxs with a neg Delta G below -58.5      ######
######         v3 - No longer iterates on MF genes, but rather joins all   ######
######              compounds, to create a compound pool.                  ######
######              This does not alter the result                         ######
######            - Takes into account all MFs of input TF and of genes    ######
######              so it searches cycles for each TF Multifun term        ######
######            - Outputs to output directory creating a file            ######
#################################################################################

#Load modules
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#get arguments
$time=localtime();
$tf=$ARGV[0];
$network=$ARGV[1];
$type=$ARGV[2];
$gu_library=$ARGV[3];
$dir_out=$ARGV[4];
chomp($dir_out);

$gene_links=$gu_library . "gene_links.txt";
$Multifun=$gu_library . "gene_GenProtEC.txt";
if($type eq "Regulon"){
    $out_file=$dir_out . $tf . "_Cycles_Report_Regulon_o1_v3.txt";
}
else{
    $out_file=$dir_out . $tf . "_Cycles_Report_o1_v3.txt";
}

if(($network) && ($tf) && ($type) && ($gu_library) && ($dir_out)){

  #Open outfile
  open(OUT,">$out_file") || die "Cannot open OUT at $out_file.\n";
  print OUT"# From $0 on $time\n# Input Network= $network\n";
  print OUT"# Genes marked with * participate in a reaction with a DeltaG less than -58.5\n";
  # Get ID of Gene
  $tf_id = &IDS($tf);
  print OUT"# TF name:\t",$tf,"\n","# Gene ID:\t",$tf_id,"\n\n";

  # Sacar los genes regulados por el input TF
  open(NT,$network) || die "Cannot open NT at $network.\n";
  while(<NT>){
      if($_=~/^\#/){
          next;
      }
      if($_=~/^$tf\t([^\t]+)$/i){
          $genelist=$1;
          chomp($genelist);
          chop($genelist);
          chop($genelist);
          @allgenelist=split(/\/\//,$genelist);
      }
  }
  close(NT);

  # Parse Metabolism multifun and ID
  %MFs = &GET_MULTIFUN($tf);

  foreach $MFid(keys %MFs){

        $name = $MFs{$MFid};
        print OUT"# Multifun ID:\t",$MFid,"\n# Multifun:\t",$name,"\n";

        # Get instances (genes) of Multifun
        my @instances = ();
        @instances = $cyc -> get_class_all_instances($MFid);

        my @instances_name = ();
        foreach $id(@instances){
            $name = &NAMES($id);
            push(@instances_name, $name);
        }

        #Sacar interseccion de genes regulados y de genes del Multifun dominante
        my @reg_genes_name = ();
        foreach $reg(@allgenelist){
            if ($reg ~~ @instances_name){
                if ($reg =~/$tf/i){
                    next;
                }
                push(@reg_genes_name,$reg);
            }
        }
        print OUT"# Genes regulated by $tf that are also part of the Multifun:\t",join(",",@reg_genes_name),"\n";

        ##### Find Cycles ######
        $cycle = 0;

        # Sacar sustratos y productos de los genes del multifun que son reguladors por el TF original
        # Sacar reactions of genes
        my @allcompounds = ();
        my %compounds = ();
        foreach $gene(@reg_genes_name){
            $geneid = &IDS($gene);
            @rxs = $cyc -> reactions_of_gene($geneid);
            foreach $rx(@rxs){
                @cmps = $cyc -> substrates_of_reaction($rx);
                push(@allcompounds,@cmps);
            }
        }

        @compounds{ @allcompounds } = ();
        my @unique_comps = keys %compounds;

        # Filter compounds, removing common currency metabolites
        @unique_comps = grep(!/^PROTON$|^.MP$|^.TP$|^.DP$|^WATER$|^NAD.*$|^\|*Pi\|*$|^OXYGEN-MOLECULE$|^CARBON-DIOXIDE$|^CO-A$|^PPI$|^\|Donor-H2\|$|^\|Acceptor\|$|^\|Nucleoside-Triphosphates\|$/i, @unique_comps);

        # Get reactions that use the same compunds and then get genes that catalize those reactions
        my @allcat_genes = ();
        my %genes_deltaG = ();
        my %cat_gns_d = ();
        my @crxs = ();
        my @gns = ();
        foreach $cmp(@unique_comps){
            @crxs = $cyc -> reactions_of_compound($cmp);
            foreach $rx(@crxs){
                @gns = $cyc -> genes_of_reaction($rx);
                push(@allcat_genes,@gns);
                # Get reaction's Delta G value, and check if it passes the threshold of -58.5, if it does, each gene is annotated as deltaG negative
                $DG = $cyc -> get_slot_value($rx, "GIBBS-0");
                if ($DG < -58.5){
                    foreach $gn(@gns){
                        $gn_name = &NAMES($gn);
                        $genes_deltaG{$gn_name} = 1;
                    }
                }
            }
        }

        @cat_gns_d{ @allcat_genes } = ();
        my @cat_gns = keys %cat_gns_d;

        my @new_genes = ();
        # Get new genes that use the same compounds and that are not part of the original Multifun term
        foreach $gn(@cat_gns){
            if ($gn ~~ @instances){
                next;
            }
            else{
                #Genes that are not part of the Multifun and that share metabolites
                push(@new_genes,$gn);
            }
        }
        
        # Get regulators of the new genes that share metabolites with the original tf-regulated gene
        my @direct_cycle_genes = ();
        my %indirect_cycle_genes;


        foreach $gn(@new_genes){
            $gn_name = &NAMES($gn);
            @regs = &GET_REGULATORS($gn_name);
            if ($tf ~~ @regs){
                $cycle = 1;
                push(@direct_cycle_genes, $gn_name);
            }
            foreach $reg(@regs){
                if ($reg eq $tf){ # Nac is part of the direct regulators
                    next;
                }
                #@regs2 = $cyc -> genes_regulating_gene($reg); #Who is regulating the regulator
                @regs2 = &GET_REGULATORS($reg);
                if ($tf ~~ @regs2){ #If nac regulates the regulator
                    $cycle = 1;
                    push(@{$indirect_cycle_genes{&NAMES($gn)}}, $reg);
                }
            }
        }

        @direct_genes{ @direct_cycle_genes } = ();
        my @direct_genes = keys %direct_genes;

        #print "DeltaG genes\t",join(",", keys %genes_deltaG),"\n";
        @dg_genes = keys %genes_deltaG;
        if ($cycle){
            print OUT"# Compounds of $tf:\t",join(",",@unique_comps),"\n";
            print OUT"# Interactions and associated Multifun: \n";
            foreach $dir(@direct_genes){
                %mf =  &GET_MULTIFUN($dir);
                foreach $key(keys %mf){
                    if ($dir ~~ @dg_genes){
                    print OUT$tf,"\t",$dir,"*\t",$key,":",$mf{$key},"\n";
                    }
                    else {
                    print OUT$tf,"\t",$dir,"\t",$key,":",$mf{$key},"\n";
                    }
                }
            }
            foreach $newgene(keys %indirect_cycle_genes){
                %mf = &GET_MULTIFUN($newgene);
                foreach $reg(@{$indirect_cycle_genes{$newgene}}){
                    #print $reg,"\t",$newgene,"\t";
                    #%mf = &GET_MULTIFUN($reg);
                    foreach $key(keys %mf){
                        if ($newgene ~~ @dg_genes){
                            print OUT$reg,"\t",$newgene,"*\t",$key,":",$mf{$key},"\n";
                        }
                        else {
                            print OUT$reg,"\t",$newgene,"\t",$key,":",$mf{$key},"\n";
                        }
                    }
                }
            }
        }
        print OUT"\n";
    }
} 
else {
  die "\nArgument usage:\n\n[1] - Query TF. \n[2] - Path to Regulatory Network. Genes separated by //.\n[3] - Type of Network: Regulon for RegulonDB network.\n[4] - Path to gu_library dir with gene_links and MultiFunTerms (GenProtEC).\n[5] - Output directory.\n\n";
}


############## FUNCTIONS ###########################
sub IDS {
 
  local($name,$id)=($_[0],"");
  open(IDS, $gene_links) || die "Cannot open IDS at $gene_links.\n";
  while(<IDS>){
    if($_=~/^([^\t]+)\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t$name$/i){
      $id=$1;
      close(IDS);
      return($id);
    }
  }
  close(IDS);
  print"WARNING 1! Gene $name is not present in gene_links.\n";
  return($name);
}
##################################################

sub NAMES {
 
  local($id,$name)=($_[0],"");
  open(NAMES, $gene_links) || die "Cannot open NAMES at $gene_links.\n";
  while(<NAMES>){
    if($_=~/^$id\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]+)$/i){
      $name=$1;
      chomp($name);
      close(NAMES);
      return($name);
    }
  }
  close(NAMES);
  print"WARNING 2! Gene $id is not present in gene_links.\n";
  return($id);
}
##################################################

sub GET_MULTIFUN {

  local($tf,%MF,@instances)=($_[0],(),());
  open(MF, $Multifun) || die "Cannot open MF at $Multifun.\n";
  while(<MF>){
    if($_=~/^\#/){
      next;
    }
    if($_=~/^$tf\t[1|2|3|4|5|6|8|9] --> (metabolism|cell processes|cell structure|DNA sites|extrachromosomal|information transfer|regulation|transport).*\s(.*)\s-->\s(.*)$/i){
      $num = $2;
      $name = $3;
      $MFid = "BC-".$num;
      @mfinstances = $cyc -> get_class_all_instances($MFid);
      $name = $name . "(" . scalar(@mfinstances) . ")";
      $MF{$MFid} = $name;
    }
  }
  close(MF);
  return(%MF);
}
##################################################

sub GET_REGULATORS {

  local($gn,@Regulators)=($_[0],());
  open(NT, $network) || die "Cannot open NT at $network";
  while(<NT>){
    if($_=~/^\#/){
      next;
    }
    if($_=~/^([^\t]+)\t.*$gn\/\//i){
      $Regulator = $1;
      push(@Regulators,$Regulator);
    }
  }
  close(NT);
  return(@Regulators);
}
##################################################