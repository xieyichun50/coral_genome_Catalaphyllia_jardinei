## Get files to the directory
## In Orthofinder/Results_Feb19/
ls -l Species_Tree/SpeciesTree_rooted_node_labels.txt
ls -l MultipleSequenceAlignments/SpeciesTreeAlignment.fa
ls -l Orthogroups/Orthogroups.GeneCount.tsv
##cp to cafe working directory

samtools faidx SpeciesTreeAlignment.fa
cat SpeciesTreeAlignment.fa.fai

python2 /tools/CAFE5/python_scripts/cafetutorial_prep_r8s.py -i SpeciesTree_rooted_node_labels.txt -o r8s_ctl_file.txt -s 484467 -p 'Hydra_vulgaris,Pachyseris_speciosa' -c '617'
/tools/r8s/r8s1.81/src/r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > r8s.ultrametric.tre
cat r8s.ultrametric.tre

awk 'OFS="\t" {$NF=""; print}' Orthogroups.GeneCount.tsv > tmp && awk '{print "(null)""\t"$0}' tmp > cafe.input.tsv && sed -i '1s/(null)/Desc/g' cafe.input.tsv && rm tmp
python3 /tools/CAFE5/python_scripts/cafetutorial_clade_and_size_filter.py -i cafe.input.tsv -o cafe.input.filter.tsv -s

##run cafe
echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda1 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda1 -c 38

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda2 -k 2 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda2 -k 2 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda3 -k 3 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda3 -k 3 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda4 -k 4 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda4 -k 4 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda5 -k 5 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda5 -k 5 -c 180

echo "/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda6 -k 6 -c 180"
/tools/CAFE5/bin/cafe5 -i cafe.input.filter.tsv -t r8s.ultrametric.tre -o r8s_filter_lambda6 -k 6 -c 180
