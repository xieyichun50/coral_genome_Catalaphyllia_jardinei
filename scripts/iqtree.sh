##use iqtree model 3
##3. Combine ModelFinder, tree search, ultrafast bootstrap and SH-aLRT test:
##     iqtree2 -s example.phy --alrt 1000 -B 1000
## model search by default -m MFP
cat genefile | while read i;
do 
echo "iqtree -s $i.align.fa --seqtype AA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000"
iqtree -s $i.align.fa --seqtype AA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000 
done
