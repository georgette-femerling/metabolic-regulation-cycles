#!/bin/bash

# Arguments
TF_list=$1
New_nt=$2
Regulon_nt=$3
lib=$4
Output_dir=$5

mkdir $Output_dir $Output_dir/Reports $Output_dir/Tables

while read tf; do 
    echo -e "--------------$tf-----------------\n\n"
    perl bin/get_cycles.v3.pl $tf $New_nt New $lib $Output_dir/Reports/
    perl bin/get_cycles.v3.pl $tf $Regulon_nt Regulon $lib $Output_dir/Reports/
    python bin/get_cycles_output2.v2.py $tf $Output_dir/Reports/${tf}_Cycles_Report_o1_v3.txt $Output_dir/Reports/${tf}_Cycles_Report_Regulon_o1_v3.txt $Output_dir/Tables/
done < $TF_list


