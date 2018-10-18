### aracne-ap tutorial ###
# https://sourceforge.net/p/aracne-ap/wiki/Home/

module load aracne_ap

#calculate threshold with a fixed seed
java -Xmx5G -jar $ARACNE_HOME/Aracne.jar -e inputs/development_data_150samples.txt -o outputs/output_MsigDB_geneID_development_150samples --pvalue 1E-8 --seed 1 --calculateThreshold

# run 100 bootstraps
for i in {1..100}
do
java -Xmx5G -jar $ARACNE_HOME/Aracne.jar -e inputs/development_data_150samples.txt -o outputs/output_MsigDB_geneID_development_150samples --tfs inputs/MsigDB_TFs.txt --pvalue 1E-8 --threads 100
done

# consolidate bootstraps in the output folder
java -Xmx50G -jar $ARACNE_HOME/Aracne.jar -o outputs/output_MsigDB_geneID_development_150samples --consolidate

