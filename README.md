# NeTFactor
NeTFactor is an R package for identifying transcription factors (TFs) (or regulators in a generic network) that are most
likely regulating a given set of biomarkers. In the most ideal scenario, NeTFactor uses a computationally inferred context-specific
gene regulatory network and applies statistical and optimization methods to this network to identify a sparse set of
regulator TFs. 

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [NeTFactor](#netfactor)
- [Reverse engineering a context-specific gene regulatory network](#reverse-engineering-a-context-specific-gene-regulatory-network)
- [Run the NeTFactor pipeline](#run-the-netfactor-pipeline)
	- [Required libraries](#required-libraries)
	- [Running the algorithm](#running-the-algorithm)
- [Citation](#citation)

<!-- /TOC -->


# Reverse engineering a context-specific gene regulatory network

One of the pre-requisites to run NeTFactor is a Gene Regulatory Network (GRN). Such a GRN consists of directed edges denoting interactions between regulators (e.g. TFs) and their target(s) (e.g. gene(s) they regulate). NeTFactor utilizes the structure and constituents of such a GRN to identify the regulators, specifically TFs, that most significantly regulate the genes underlying the biomarker. We used the ARACNe-AP implementation of the ARACNe algorithm (can be downloaded from https://github.com/califano-lab/ARACNe-AP) to create a context-specific GRN from a nasal RNAseq data of a case-control asthma cohort. For reproducibility we share the network created from this data (see ```data/input/aracne_network.txt```). Next, we want to give commands used to generate the network. First, to run the ARACNe algorithm one requires a gene-expression matrix as well as a set of TFs (or potential regulators) as input. The gene expression matrix is a tab-delimited file formatted such that genes are  rows and samples are columns.
Moreover, the first row contains the gene names and the first column contains the sample names. Here is a simple example data:

|         | Sample1          | Sample2  |
| ------------- |:-------------:| -----:|
| Gene1     | 1.3 | 2|
| Gene2     | 6      |   10 |
| Gene3| 4      |    15|


The TF list  file is also a tab-delimited file such that each row contains the gene symbol (should match with the names in the gene expression matrix). Here is a simple example:


|    TF1     |
| ------------- |
| TF2   |
| TF2     |
| TF3|

For the current application,
we obtained 221 putative TFs from the Molecular Signature Database (MSigDB) version 5.1 (accessed June 3rd, 2016). Of these, 132 were included in the nasal RNAseq dataset and we retained those for further analysis.
We have already shared the expression data [here](https://www.synapse.org/#!Synapse:syn9878922/files/)  and the TF list is shared with this github directory (see ```data/input/TF_list.txt```). Next, the ARACNe inferred network can be inferred using the following command:

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
After these three simple codes one has the desired network in the outputs folder. For the ease of use, we share the network we used in  (see``` data/input/aracne_network.txt)```.
Note that all these steps are for running ARACNe on a unix based system and for other platforms please refer to the ARACNe-AP github [page](https://github.com/califano-lab/ARACNe-AP).

# Run the NeTFactor pipeline

## Required libraries
NeTFactor requires the following R packages to be installed:
* [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html)
* [CVXR](https://cran.r-project.org/web/packages/CVXR/index.html)
* [mixtools](https://cran.r-project.org/web/packages/mixtools/index.html)
* [Biobase](https://bioconductor.org/packages/release/bioc/html/Biobase.html)
* [BiocGenerics](https://bioconductor.org/packages/release/bioc/html/BiocGenerics.html)

All the above packages are available either on [Bioconductor](https://www.bioconductor.org/) or [CRAN](https://cran.r-project.org/) and can be easily installed freely. Please refer to the package websites for detailed installation instructions.






## Running the algorithm
Apart from the network and gene expression, we need a phenotype table which lists the phenotype for the samples in the gene expression data as well as the list of biomarkers genes. The phenotype table is a tab-delimited file which two columns, where the first column have the sample names (this should match with the sample names in the gene expression file) and the second column contains phenotype information. We shared the phenotype data used in our application within this github folder (see ```data/input/pheno_for_netfactor.txt```) and a simple example file should look like:



|      Sample Name   | Phenotype
| ------------- |:-------------:|
| Sample1    | A |
| Sample2     | N      |
| Sample3| A     |


The biomarker file (see ```data/input/biomarker_genes.txt``` for the one used in our application) is a tab-delimited file with a single column that contains the names of the biomarker genes. A simple example file should look like this:



|      Gene_Name  |
| ------------- |
| Biomarker_gene_1   |
| Biomarker_gene_2     |
| Biomarker_gene_3 |


With all the data needed we can run the NeTFactor pipeline with the following simple commands:

```
source("run_netfactor.R")
path_to_expression="data/input/nasal_expression.txt"
path_to_network="data/input/aracne_network.txt"
path_to_pheno="data/input/pheno_for_netfactor.txt"
path_to_biomarker="data/input/biomarker_genes.txt"
netfactor_results=run_netfactor(path_to_expression,"network_100",path_to_network,
                   path_to_pheno,path_to_biomarker)
```
The data.frame ```netfactor_results``` contains the results of the netfactor algorithm. It has 5 columns with gene names as row names of the data.frame.
The top-7 rows for the asthma data is as follows:
<center>

|         |Lasso Weight  | Genes_Regulated  | Cumulative_Genes_Covered|Viper_FDR|Fisher_Regulon_FDR
| ------------- |:-------------:| -----:|---:|---:|---:|
| PPARG     | 1.072180| 16 |16 |0.0151|0.0061|
| ETV4     | 1.0150     |   15 | 24|0.0328|0.0000679|
| GTF2A2| 1.0130      |    11|30|0.00793|0.3951|
| EGR1|   1.0060    |    3|  33 | 0.568  | 0.2856  |
|SPI1   |  1.0020 |  8 |  38 | 0.606  |  0.3633 |
|CEBPB   | 1.0007  |  7 | 41 |  0.89 |  0.0620 |
|XBP1   |  1.0004 |  11 |  52 |  0.99 |  0.0002 |
</center>
Above the first column has the weights coming from our lasso based convex optimization, second column has the number of biomarker genes regulated by each TF, third column has the number of genes cumulatively regulated by it and all the TFs preceding it, fourth column has the FDR associated with viper inferred activity and the fifth column has the FDR associated with the statistical significance of the overlap between the set of genes regulated by each TF in the network and the biomarker gene set.

# Citation
Ahsen, M.E., Chun, Y., Grishin, A. et al. NeTFactor, a framework for identifying transcriptional regulators of gene expression-based biomarkers. Sci Rep 9, 12970 (2019) doi:[10.1038/s41598-019-49498-y](https://www.nature.com/articles/s41598-019-49498-y)
