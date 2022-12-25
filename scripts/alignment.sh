cat genefile | while read i;
do
echo "Aligning $i"
mafft --thread 38 --auto $i.fa > $i.align.fa
done
