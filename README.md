# NeTFactor
NeTFactor is an R code that tries to identify transcription factors (TFs) (or regulators in a generic network) that are most
likely regulating a given set of biomarkers. In the most ideal scenario, NeTFactor uses a computationally inferred context-specific
gene regulatory network and applies statistical and optimization methods to this network to identify a sparse set of
regulator TFs. For more info please see [1] and cite if you use NeTFactor.

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [NeTFactor](#netfactor)
- [Reverse Engineer a Context-Specific Network (Optional)](#reverse-engineer-a-context-specific-network-optional)
- [Run the NeTFactor Pipeline](#run-the-netfactor-pipeline)
	- [Required Libraries](#required-libraries)
	- [Running the Algorithm](#running-the-algorithm)
- [Plots](#plots)
- [Citation](#citation)

<!-- /TOC -->

# Reverse Engineer a Context-Specific Network (Optional)
One of the pre-requisites to run NeTFactor is a Gene Regulatory Network (GRN). Such a GRN consists of directed edges denoting interactions between regulators (e.g. TFs) and their target(s) (e.g. gene(s) they regulate). NeTFactor utilizes the structure and constituents of such a GRN to identify the regulators, specifically TFs, that most significantly regulate the genes underlying the biomarker. For more accurate results, we suggest to reverse engineer a context-specific network if there is available contextual gene expression data. However, NeTFactor is a generic algorithm and the user can use any regulatory network available in the literature. For the companion application in [1], we used the ARACNe-AP implementation of the ARACNe algorithm (can be downloaded from https://github.com/califano-lab/ARACNe-AP) algorithm to create a context-specific GRN from a nasal RNAseq data of a case-control asthma cohort. For reproducibility we share the network created from this data (see ```data/input/aracne_network.txt```). Next, we want to give commands used to generate the network. First, to run the ARACNe algorithm one requires a gene-expression matrix as well as a set of TFs (or potential regulators) as input. The gene expression matrix is a tab-delimited file formatted such that genes are  rows and samples are columns. Moreover, the first row contains the gene names and the first column contains the sample names.
The TF list  file is also a tab-delimited file such that each row contains the gene symbol (should match with the names in the gene expression matrix). For the current application,
we obtained 221 putative TFs from the Molecular Signature Database (MSigDB) version 5.1 (accessed June 3rd, 2016). Of these, 132 were included in the nasal RNAseq dataset and we retained those for further analysis.
We share the expression data ```data/input/nasal_expression.txt``` and TF list  ```data/input/TF_list.txt``` used in our paper [1]. Next, the ARACNe inferred network can be inferred using the following command:

1. Calculate threshold with a fixed seed:
```
java -Xmx5G -jar $ARACNE_HOME/Aracne.jar -e input/nasal_expression.txt -o outfolder/ --pvalue 1E-8 --seed 1 --calculateThreshold
```
2.  Run 100 Bootstraps (can be modified)
```
for i in {1..100}
do
java -Xmx5G -jar $ARACNE_HOME/Aracne.jar -e input/nasal_expression.txt -o outfolder/ --tfs input/TF_list.txt --pvalue 1E-8 --threads 100
done
```
3. Consolidate bootstraps in the output folder
```
java -Xmx50G -jar $ARACNE_HOME/Aracne.jar -o outputs/ --consolidate
```
After these three simple codes one has the desired network in the outputs folder. For the ease of use, we share the network we used in [1] ```(see data/outfolder/aracne_network.txt)```.
Note that all these steps are optimized for running ARACNe on a unix based system and for other platforms please refer to the ARACNe-AP github [page](https://github.com/califano-lab/ARACNe-AP).

# Run the NeTFactor Pipeline

## Required Libraries
NeTFactor requires the following R packages to be installed:
* Viper
* RCVX


## Running the Algorithm
After we have the main ingredient for the NeTFactor (i.e. the GRN), we can continue running the main pipeline using the following simple command:

# Plots


# Citation
[1] **Network analyses identify transcription factor regulators of a novel gene expression-based asthma biomarker,** Mehmet Eren Ahsen, Alexander Grishin,Yoojin Chun, Galina Grishina, Gustavo Stolovitzky, Gaurav Pandey and Supinda Bunyanovich, *Journal of Cool Biology*, 2018 (hopefully).
