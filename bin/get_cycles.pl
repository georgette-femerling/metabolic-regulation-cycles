#################################################################################
###### 08-25-2020                                                          ######
###### Author: Georgette Femerling                                         ######
######                                                                     ######
###### Gets compounds of pathways and gets pathways with a distance of 0   ######
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
#################################################################################


#Load modules
use perlcyc;
$cyc = perlcyc -> new("ECOLI");

#get arguments
$time=localtime();
$tf=$ARGV[0];
$network=$ARGV[1];
$gu_library=$ARGV[2];
chomp($gu_library);
#$dir_out=$ARGV[3];
#chomp($dir_out);

$gene_links=$gu_library . "gene_links.txt";
$Multifun=$gu_library . "gene_GenProtEC.txt";


if(($network) && ($tf)){

  # Parse Metabolism multifun and ID
  %MFs = &GET_MULTIFUN($tf);

  foreach $MFid(keys %MFs){

      # Get ID of Gene
      $tf_id = &IDS($tf);
      print "TF name:\t",$tf,"\n","Gene ID:\t",$tf_id,"\n";

      $name = $MFs{$MFid};
      print "Multifun ID:\t",$MFid,"\nMutifun:\t",$name,"\n";

      # Get instances (genes) of Multifun
      @instances = $cyc -> get_class_all_instances($MFid);

      foreach $id(@instances){
        $name = &NAMES($id);
        push(@instances_name, $name);
      }

      # Sacar los genes regulados por el input TF
      open(NT,$network) || die "Cannot open NT at $network.\n";
        while(<NT>){
          if($_=~/^\#/){
            next;
          }

          if($_=~/^$tf\t([^\t]+)$/){
            $genelist=$1;
            chomp($genelist);
            chop($genelist);
            chop($genelist);
            @allgenelist=split(/\/\//,$genelist);
          }
        }
      close(NT);

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
      print "Genes regulated by $tf that are also part of the Multifun:\t",join(",",@reg_genes_name),"\n\n\n";

      ##### Find Cycles ######
      $cycle = 0;

      # Sacar sustratos y productos de los genes del multifun que son reguladors por el TF original

      # Sacar reactions of genes
      my @allcompounds = ();
      foreach $gene(@reg_genes_name){
        $geneid = &IDS($gene);
        @rxs = $cyc -> reactions_of_gene($geneid);
        foreach $rx(@rxs){
          @cmps = $cyc -> substrates_of_reaction($rx);
          push(@allcompounds,@cmps);
        }

        @compounds{ @allcompounds } = ();
        my @unique_comps = keys %compounds;

        # Filter compounds, removing common currency metabolites
        @unique_comps = grep(!/^PROTON$|^.MP$|^.TP$|^.DP$|^WATER$|^NAD.*$|^\|*Pi\|*$|^OXYGEN-MOLECULE$|^CARBON-DIOXIDE$|^CO-A$|^PPI$/i, @unique_comps);

        # Get reactions that use the same compunds and then get genes that catalize those reactions
        my @allcat_genes = ();
        foreach $cmp(@unique_comps){
          @crxs = $cyc -> reactions_of_compound($cmp);
          foreach $rx(@crxs){
            @gns = $cyc -> genes_of_reaction($rx);
            push(@allcat_genes,@gns);
          }
        }
        @cat_gns{ @allcat_genes } = ();
        my @cat_gns = keys %cat_gns;

        my @new_genes = ();
        
        # Get new genes that use the same compounds and that are not part of the original Multifun term
        foreach $gn(@cat_gns){
          if ($gn ~~ @instances){
            next;
          }
          else{
            #print "not in instances:",&NAMES($gn), "\n";
            #Genes that are not part of the Multifun and that share metabolites
            push(@new_genes,$gn);
          }
        }
        
        # Get regulators of the new genes that share metabolites with the original tf-regulated gene
        my @direct_cycle_genes = ();
        my %indirect_cycle_genes;

        foreach $gn(@new_genes){
          @regs = $cyc -> genes_regulating_gene($gn); # Regulators of the new genes
          if ($tf_id ~~ @regs){
            $cycle = 1;
            #print "C1\t",$tf,"\t",&NAMES($gn),"\n"; # If Nac directly regulates a new gene
            push(@direct_cycle_genes, &NAMES($gn));
          }
          foreach $reg(@regs){
            if ($reg eq $tf_id){ # Nac is part of the direct regulators
              next;
            }
            @regs2 = $cyc -> genes_regulating_gene($reg); #Who is regulating the regulator
            if ($tf_id ~~ @regs2){ #If nac regulates the regulator
              $cycle = 1;
              #print "C2\t",$tf,"\t",&NAMES($reg),"\t",&NAMES($gn),"\n";
              push(@{$indirect_cycle_genes{&NAMES($gn)}}, &NAMES($reg));
            }
          }
        }

        @direct_genes{ @direct_cycle_genes } = ();
        my @direct_genes = keys %direct_genes;

        if ($cycle){
          print "Start gene: ",$gene,"\n";
          print "Compunds of gene:\t",join(",",@unique_comps),"\n";
          print "Genes directly regulated by TF that share compounds:\n";
          foreach $dir(@direct_genes){
            %mf =  &GET_MULTIFUN($dir);
            foreach $key(keys %mf){
              print $dir,"\t",$key,":",$mf{$key},"\n";
            }
          }
          print "Indirect regulation of genes that share compunds:\n";
          foreach $newgene(keys %indirect_cycle_genes){
            %mf1 = &GET_MULTIFUN($newgene);
            foreach $reg(@{$indirect_cycle_genes{$newgene}}){
              print $reg,"\t",$newgene,"\t";
              %mf = &GET_MULTIFUN($reg);
              foreach $key(keys %mf){
                print $key,":",$mf{$key},",";
              }
              print "\t";
              foreach $key(keys %mf1){
                print $key,":",$mf1{$key},",";
              }
              print "\n";
            }
          }
        }
        print "\n\n"
     }
  }
} 
else {
  die "\nArgument usage:\n\n[1] - Query TF. \n[2] - Path to Regulatory Network. Genes separated by //.\n[3] - Path to gu_library dir with gene_links and MultiFunTerms (GenProtEC).\n\n";
}


########################################
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
  print"WARNING 2! Gene $name is not present in gene_links.\n";
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

sub GET_MULTIFUN {

  local($tf,%MF)=($_[0],());
  open(MF, $Multifun) || die "Cannot open MF at $Multifun.\n";
  while(<MF>){
    if($_=~/^\#/){
      next;
    }
    if($_=~/^$tf\t[1|5] --> (metabolism|cell processes).*\s(.*)\s-->(.*)$/i){
      $num = $2;
      $name = $3;
      $MFid = "BC-".$num;
      $MF{$MFid} = $name;
    }
  }
  close(MF);
  return(%MF);
}