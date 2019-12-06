#!/bin/bash

# This version computes the interaction and binding coordinates for core Enh/Pro and first neighbors (e.g. core 2 -> neighbor1 -> core 1) where first neighbor is also regulated by core TF

echo "In this version we also compute Accessible Enh and Pro interactions with -f 0.5 and -F 0.5 mean 50% of TF peak is within DNase peak OR 50% of DNase peak is within ChIP-seq peak"

# sh code_v2.sh CardiacMC

in_dir=/work/users/mali/Tasks/core_networks/Evan_networks/$1/
out_dir=/work/users/mali/Tasks/core_networks/Evan_networks/$1/
K27ac='_K27ac'
K4me3='_K4me3'
mRNA='_mRNA'
DNaseFile='_DNase'

cd $in_dir
pwd
for i in *_TF.txt
do
name=$(echo $i | cut -d "_" -f 1)
echo $name

# getting enhancer file
eFile=$name$K27ac
#echo $eFile
enhFile="$eFile*"
#echo $enhFile
enhancerFile=$(ls $enhFile)
echo $enhancerFile

# getting promoter file
pFile=$name$K4me3
proFile="$pFile*"
promoterFile=$(ls $proFile)
echo $promoterFile

# getting expressed TFs file
expTFs=$name$mRNA
expTFFile="$expTFs*"
exprFile=$(ls $expTFFile)
echo $exprFile

# getting DNAse file
acc=$name$DNaseFile
accFile="$acc*"
DNase=$(ls $accFile)
echo $DNase

#: <<'end_long_comment'

cd ${out_dir}

awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' ${in_dir}$i /work/users/mali/Tasks/core_networks/Evan_networks/input_files/Processed_GeneHancer_GRCh38.bed > "${name}_Enhancers.tmp"
wc -l "${name}_Enhancers.tmp"
awk '(NR>1) && ($8 > 10 ) ' "${name}_Enhancers.tmp" > filtered_Enh.tmp # Remove those enhancers-gene association where score is less than 10
sort -k1,1 -k2,2n filtered_Enh.tmp > "sorted_${name}_Enhancers.tmp"
wc -l "sorted_${name}_Enhancers.tmp"

bedtools intersect -a "sorted_${name}_Enhancers.tmp" -b ${in_dir}$enhancerFile > "${name}_Act_Enhancers.tmp"
wc -l "${name}_Act_Enhancers.tmp"
sort -k1,1 -k2,2n "${name}_Act_Enhancers.tmp" > "sorted_${name}_Act_Enhancers.tmp"

# -F 1.0 will look for "entire" TFBS peak falling completely within Active Enhancer region.
bedtools intersect -a "sorted_${name}_Act_Enhancers.tmp" -b /work/users/mali/Tasks/core_networks/Evan_networks/input_files/Cistrome_ChIPseq_sorted.bed -F 1.0 -wo -sorted > "${name}_All_TFs_In_Core_Enh.tmp"
wc -l "${name}_All_TFs_In_Core_Enh.tmp"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$i "${name}_All_TFs_In_Core_Enh.tmp" > "${name}_CoreTFs_In_CoreEnh.tmp"
wc -l "${name}_CoreTFs_In_CoreEnh.tmp"

#### Getting first neighbours

awk 'FNR==NR {a[$1]; next} !($(NF-2) in a)' ${in_dir}$i "${name}_All_TFs_In_Core_Enh.tmp" > NC_TF_In_CoreEnh.tmp
wc -l NC_TF_In_CoreEnh.tmp
cat NC_TF_In_CoreEnh.tmp | awk '{print $12}' | sort -n | uniq > NC_Uniq_TFs_Binding_CoreEnh.tmp
wc -l NC_Uniq_TFs_Binding_CoreEnh.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' NC_Uniq_TFs_Binding_CoreEnh.tmp /work/users/mali/Tasks/core_networks/Evan_networks/input_files/Processed_GeneHancer_GRCh38.bed > NC_TFs_Enh.tmp
awk '(NR>1) && ($8 > 10)' NC_TFs_Enh.tmp > Filtered_NC_TFs_Enh.tmp # Remove those enhancers-gene association where score is less than 10
wc -l NC_TFs_Enh.tmp
wc -l Filtered_NC_TFs_Enh.tmp

# get active non-core TFs enhancers which are first neigbors
sort -k1,1 -k2,2n Filtered_NC_TFs_Enh.tmp > Filtered_NC_TFs_Enh_Sorted.tmp
bedtools intersect -a Filtered_NC_TFs_Enh_Sorted.tmp -b ${in_dir}$enhancerFile > Filtered_NC_TFs_Act_Enh.tmp
wc -l Filtered_NC_TFs_Act_Enh.tmp
sort -k1,1 -k2,2n Filtered_NC_TFs_Act_Enh.tmp > Filtered_NC_TFs_Act_Enh_Sorted.tmp

# get TF bindings in non-core active enhancers
bedtools intersect -a Filtered_NC_TFs_Act_Enh_Sorted.tmp -b /work/users/mali/Tasks/core_networks/Evan_networks/input_files/Cistrome_ChIPseq_sorted.bed -F 1.0 -wo -sorted > "${name}_All_TFs_In_NonCore_Enh.tmp"
wc -l "${name}_All_TFs_In_NonCore_Enh.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$i "${name}_All_TFs_In_NonCore_Enh.tmp" > "${name}_Final_CoreTFs_In_NonCore_Enh.tmp"
wc -l "${name}_Final_CoreTFs_In_NonCore_Enh.tmp"

# get non-core TF names whose enhancers are bound by core enhancers
cat "${name}_Final_CoreTFs_In_NonCore_Enh.tmp" | awk '{print $7}' | sort -n | uniq > NC_Uniq_TF_Names_Bound_By_CoreTF.tmp
wc -l NC_Uniq_TF_Names_Bound_By_CoreTF.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${name}_All_TFs_In_Core_Enh.tmp" > "${name}_Final_NonCoreTFs_In_Core_Enh.tmp"
wc -l "${name}_Final_NonCoreTFs_In_Core_Enh.tmp"

########################
# Uncomment this chunk if you want to get interactions among non-core TFs as well (e.g. Neigbor1 -> Neighbor 3)
# Comment the chunk below in this case
# Chunk Start
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${name}_All_TFs_In_NonCore_Enh.tmp" > "${name}_Final_NonCoreTFs_In_NonCore_Enh.tmp"
wc -l "${name}_Final_NonCoreTFs_In_NonCore_Enh.tmp"
# concatenating all into one
cat "${name}_CoreTFs_In_CoreEnh.tmp" "${name}_Final_NonCoreTFs_In_Core_Enh.tmp" "${name}_Final_CoreTFs_In_NonCore_Enh.tmp" "${name}_Final_NonCoreTFs_In_NonCore_Enh.tmp" > shell1_Enh.tmp
wc -l shell1_Enh.tmp
# Chunk End

################
# Uncomment this chunk if you DO NOT want to get interactions among non-core TFs (e.g. Neigbor1 -> Neighbor 3)
# comment the above chunk
# Chunk Start
# concatenating all into one
#cat "${name}_CoreTFs_In_CoreEnh.tmp" "${name}_Final_NonCoreTFs_In_Core_Enh.tmp" "${name}_Final_CoreTFs_In_NonCore_Enh.tmp" > shell1_Enh.tmp
#wc -l shell1_Enh.tmp
# Chunk End

################

# removing TFs which are not expressed in the given cell type
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$exprFile shell1_Enh.tmp > shell1_Enh_expTF.tmp
wc -l shell1_Enh_expTF.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-7) in a)' ${in_dir}$exprFile shell1_Enh_expTF.tmp > shell1_Enh_expTF_2.tmp
wc -l shell1_Enh_expTF_2.tmp

# processing the files to give output in the desired format (ChIP-peak, active enhancer region)
awk 'BEGIN {FS=OFS="\t"} {print $9,$10,$11,$12,$13,$1,$2,$3,$4,$5,$6,$7,$8,$14}' shell1_Enh_expTF_2.tmp > "${name}_Shell1_Enh.bed"
wc -l "${name}_Shell1_Enh.bed"

# getting interactions which are in active regulatory regions (accessibility not applied yet)
#cat "${name}_Shell1_Enh.bed" | awk '{print $4"\t"$12}' | sort -n | uniq > "${name}_Shell1_Enh_Int.txt"
#wc -l "${name}_Shell1_Enh_Int.txt"

# getting accessible bed files and interactions for shell1
# -f 0.5 means 50% of Enhancer Region (ChIP peak in Enh.) should be falling within the accessibility peak
bedtools intersect -a "${name}_Shell1_Enh.bed" -b ${in_dir}$DNase -f 0.50 -F 0.50 -e -wa > "${name}_Shell1_Enh_Acc.bed"
# remove non-TF rows from the bed file
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' /work/users/mali/Tasks/core_networks/Evan_networks/input_files/AnimalTFDB_TF.bed "${name}_Shell1_Enh_Acc.bed" > "${name}_Shell1_Enh_Acc_OnlyTF.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-10) in a)' /work/users/mali/Tasks/core_networks/Evan_networks/input_files/AnimalTFDB_TF.bed "${name}_Shell1_Enh_Acc_OnlyTF.tmp" > "${name}_Shell1_Enh_Acc_OnlyTF.bed"
wc -l "${name}_Shell1_Enh_Acc_OnlyTF.bed"
cat "${name}_Shell1_Enh_Acc_OnlyTF.bed" | awk '{print $4"\t"$12}' | sort -n | uniq > "${name}_Shell1_Enh_Int_Acc_OnlyTF.txt"
wc -l "${name}_Shell1_Enh_Int_Acc_OnlyTF.txt"

# remove for only the TF genes (originally I also had CRM and co-factors in TF gene list)
#awk 'FNR==NR {a[$1]; next} ($(NF) in a)' ./AnimalTFDB_TF.txt "${name}_Shell1_Enh_Int_Acc.tmp" > "${name}_Shell1_Enh_Int_Acc_onlyTF.tmp"
#awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' ./AnimalTFDB_TF.txt "${name}_Shell1_Enh_Int_Acc_onlyTF.tmp" > "${name}_Shell1_Enh_Int_Acc.txt"
#wc -l "${name}_Shell1_Enh_Int_Acc.txt"

# getting bed files and interactions for core TFs only
: <<'end_long_comment'
awk 'BEGIN {FS=OFS="\t"} {print $9,$10,$11,$12,$13,$1,$2,$3,$4,$5,$6,$7,$8,$14}' "${name}_CoreTFs_In_CoreEnh.tmp" > "${name}_core_Enh.bed"
wc -l "${name}_core_Enh.bed"
cat "${name}_core_Enh.bed" | awk '{print $4,$12}' | sort -n | uniq > "${name}_core_Enh_Int.txt"
wc -l "${name}_core_Enh_Int.txt"
end_long_comment

rm -f *.tmp

#######################
###### Promoters ######
#######################

#: <<'end_long_comment'
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$i /work/users/mali/Tasks/core_networks/Evan_networks/input_files/GRCh38_promoter.bed > "${name}_Promoters.tmp"
sort -k1,1 -k2,2n "${name}_Promoters.tmp" > "sorted_${name}_Promoters.tmp"
wc -l "sorted_${name}_Promoters.tmp"

# consider the entire promoter for ChIP-scan if only a single K4me3 peak is falling within this promoter region (-u option) 
bedtools intersect -a "sorted_${name}_Promoters.tmp" -b ${in_dir}$promoterFile -u > "${name}_Act_Pro.tmp"
wc -l "${name}_Act_Pro.tmp"

# -F 1.0 means entire ChIP-peak should fall within Active promoter region
bedtools intersect -a "${name}_Act_Pro.tmp" -b /work/users/mali/Tasks/core_networks/Evan_networks/input_files/Cistrome_ChIPseq_sorted.bed -F 1.0 -wo -sorted > "${name}_All_TFs_In_Core_Pro.tmp"
wc -l "${name}_All_TFs_In_Core_Pro.tmp"

awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$i "${name}_All_TFs_In_Core_Pro.tmp" > "${name}_CoreTFs_In_Core_Pro.tmp"
wc -l "${name}_CoreTFs_In_Core_Pro.tmp"


#### Getting first neighbours ####

# get unique non-core TFs bound in core promoters
awk 'FNR==NR {a[$1]; next} !($(NF-2) in a)' ${in_dir}$i "${name}_All_TFs_In_Core_Pro.tmp" > NC_TF_In_CorePro.tmp
wc -l NC_TF_In_CorePro.tmp
cat NC_TF_In_CorePro.tmp | awk '{print $10}' | sort -n | uniq > NC_Uniq_TF_Binding_CorePro.tmp
wc -l NC_Uniq_TF_Binding_CorePro.tmp

# get active promoter regions of non-core TFs
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Binding_CorePro.tmp /work/users/mali/Tasks/core_networks/Evan_networks/input_files/GRCh38_promoter.bed > NC_TFs_Promoters.tmp
wc -l NC_TFs_Promoters.tmp
sort -k1,1 -k2,2n NC_TFs_Promoters.tmp > sorted_NC_TFs_Promoters.tmp
bedtools intersect -a sorted_NC_TFs_Promoters.tmp -b ${in_dir}$promoterFile -u > NC_TFs_Act_Pro.tmp
wc -l NC_TFs_Act_Pro.tmp

# scan non-core promoters for core-TF bindings and consider only those which are bound by core-TFs.
bedtools intersect -a NC_TFs_Act_Pro.tmp -b /work/users/mali/Tasks/core_networks/Evan_networks/input_files/Cistrome_ChIPseq_sorted.bed -F 1.0 -wo -sorted > "${name}_All_TFs_In_NonCore_Pro.tmp"
wc -l "${name}_All_TFs_In_NonCore_Pro.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$i "${name}_All_TFs_In_NonCore_Pro.tmp" > "${name}_Final_CoreTFs_In_NonCore_Pro.tmp"
wc -l "${name}_Final_CoreTFs_In_NonCore_Pro.tmp"

# unique non-core TFs whose promoters are bound by core-TFs
cat "${name}_Final_CoreTFs_In_NonCore_Pro.tmp" | awk '{print $4}' | sort -n | uniq > NC_Uniq_TF_Names_Bound_By_CoreTF.tmp
wc -l NC_Uniq_TF_Names_Bound_By_CoreTF.tmp

# get core-TFs promoters bound by non-core TFs (which in turn are regulated by the core TFs)
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${name}_All_TFs_In_Core_Pro.tmp" > "${name}_Final_NonCoreTFs_In_Core_Pro.tmp"
wc -l "${name}_Final_NonCoreTFs_In_Core_Pro.tmp"

################

# Uncomment this chunk if you want to get interactions among non-core TFs as well (e.g. Neigbor1 -> Neighbor 3)
# Chunk Start
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' NC_Uniq_TF_Names_Bound_By_CoreTF.tmp "${name}_All_TFs_In_NonCore_Pro.tmp" > "${name}_Final_NonCoreTFs_In_NonCore_Pro.tmp"
wc -l "${name}_Final_NonCoreTFs_In_NonCore_Pro.tmp"
# concatenating all into one
cat "${name}_CoreTFs_In_Core_Pro.tmp" "${name}_Final_NonCoreTFs_In_Core_Pro.tmp" "${name}_Final_CoreTFs_In_NonCore_Pro.tmp" "${name}_Final_NonCoreTFs_In_NonCore_Pro.tmp" > shell1_Pro.tmp
wc -l shell1_Pro.tmp
# comment out the below chunk
# Chunk End

################

# Uncomment this chunk if you DO NOT want to get interactions among non-core TFs (e.g. Neigbor1 -> Neighbor 3)
# comment the above chunk
# Chunk Start
# concatenating all into one
#cat "${name}_CoreTFs_In_Core_Pro.tmp" "${name}_Final_NonCoreTFs_In_Core_Pro.tmp" "${name}_Final_CoreTFs_In_NonCore_Pro.tmp" > shell1_Pro.tmp
#wc -l shell1_Pro.tmp
# Chunk End
################

# removing TFs which are not expressed in the given cell type
awk 'FNR==NR {a[$1]; next} ($(NF-2) in a)' ${in_dir}$exprFile shell1_Pro.tmp > shell1_pro_expTF.tmp
wc -l shell1_pro_expTF.tmp
awk 'FNR==NR {a[$1]; next} ($(NF-8) in a)' ${in_dir}$exprFile shell1_pro_expTF.tmp > shell1_pro_expTF_2.tmp
wc -l shell1_pro_expTF_2.tmp

# processing the files to give output in the desired format (ChIP-peak, active enhancer region)
awk 'BEGIN {FS=OFS="\t"} {print $7,$8,$9,$10,$11,$1,$2,$3,$4,$5,$6,$12}' shell1_pro_expTF_2.tmp > "${name}_Shell1_Pro.bed"
wc -l "${name}_Shell1_Pro.bed"
# get interactions among active regulatory regions (accessibility not used yet)
#cat "${name}_Shell1_Pro.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${name}_Shell1_Pro_Int.txt"
#wc -l "${name}_Shell1_Pro_Int.txt"

#end_long_comment
# Remove inaccessible interactions
bedtools intersect -a "${name}_Shell1_Pro.bed" -b ${in_dir}$DNase -f 0.50 -F 0.50 -e -wa > "${name}_Shell1_Pro_Acc.bed"
wc -l "${name}_Shell1_Pro_Acc.bed"
# get output columns
#cat "${name}_Shell1_Pro_Acc.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${name}_Shell1_Pro_Int_Acc.tmp"
#wc -l "${name}_Shell1_Pro_Int_Acc.tmp"

# remove for only the TF genes (originally I also had CRM and co-factors in TF gene list)
awk 'FNR==NR {a[$1]; next} ($(NF-3) in a)' /work/users/mali/Tasks/core_networks/Evan_networks/input_files/AnimalTFDB_TF.bed "${name}_Shell1_Pro_Acc.bed" > "${name}_Shell1_Pro_Acc_OnlyTF.tmp"
awk 'FNR==NR {a[$1]; next} ($(NF-8) in a)' /work/users/mali/Tasks/core_networks/Evan_networks/input_files/AnimalTFDB_TF.bed "${name}_Shell1_Pro_Acc_OnlyTF.tmp" > "${name}_Shell1_Pro_Acc_OnlyTF.bed"
wc -l "${name}_Shell1_Pro_Acc_OnlyTF.bed"
cat "${name}_Shell1_Pro_Acc_OnlyTF.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${name}_Shell1_Pro_Int_Acc_OnlyTF.txt"
wc -l "${name}_Shell1_Pro_Int_Acc_OnlyTF.txt"

# remove for only the TF genes (originally I also had CRM and co-factors in TF gene list)
#awk 'FNR==NR {a[$1]; next} ($(NF) in a)' ./AnimalTFDB_TF.txt "${name}_Shell1_Pro_Int_Acc.tmp" > "${name}_Shell1_Pro_Int_Acc_onlyTF.tmp"
#awk 'FNR==NR {a[$1]; next} ($(NF-1) in a)' ./AnimalTFDB_TF.txt "${name}_Shell1_Pro_Int_Acc_onlyTF.tmp" > "${name}_Shell1_Pro_Int_Acc.txt"
#wc -l "${name}_Shell1_Pro_Int_Acc.txt"

# getting bed files and interactions for core TFs only
: <<'end_long_comment'
awk 'BEGIN {FS=OFS="\t"} {print $7,$8,$9,$10,$11,$1,$2,$3,$4,$5,$6,$12}' "${name}_CoreTFs_In_Core_Pro.tmp" > "${name}_core_Pro.bed"
wc -l "${name}_core_Pro.bed"
# get output columns (Source Target)
cat "${name}_core_Pro.bed" | awk '{print $4"\t"$9}' | sort -n | uniq > "${name}_core_Pro_Int.txt"
wc -l "${name}_core_Pro_Int.txt"
end_long_comment

rm -f *.tmp

echo "Interactions common between enhancers and promoters"
comm -12 "${name}_Shell1_Enh_Int_Acc_OnlyTF.txt" "${name}_Shell1_Pro_Int_Acc_OnlyTF.txt" | wc -l
echo "unique source in enh"
awk '{print $1}' "${name}_Shell1_Enh_Int_Acc_OnlyTF.txt" | sort -n | uniq | wc -l 
echo "unique source in pro"
awk '{print $1}' "${name}_Shell1_Pro_Int_Acc_OnlyTF.txt" | sort -n | uniq | wc -l

cd ${in_dir}

rm -f *_Int_Acc_OnlyTF.txt
rm -f *_Acc.bed

done
pwd
