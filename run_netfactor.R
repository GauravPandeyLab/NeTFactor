run_netfactor <-function(dummy_exp,n_name,dummy_net,d_pheno,
                          gene_set,plotx,dolasso)
{


if(length(intersect(colnames(dummy_net),c("Correlation","FDR")))<1)
{
  ncol_net=dim(dummy_net)[2]
  dummy_net=cbind(dummy_net,0)
  colnames(dummy_net)[ncol_net+1]="Correlation"
  dummy_net=cbind(dummy_net,0)
  colnames(dummy_net)[ncol_net+2]="FDR"
  for(i in 1:length(dummy_net[,1]))
  {
    dummy=cor.test(dummy_exp[dummy_net[i,1],],dummy_exp[dummy_net[i,2],],method="spearman")
    dummy_net[i,ncol_net+1]=as.numeric(dummy$estimate)
    dummy_net[i,ncol_net+2]=as.numeric(dummy$p.value)
  }
  dummy_net[,ncol_net+2]=p.adjust(dummy_net[,ncol_net+2],method = "BH")
  #assign(n_name,dummy_net,envir = .GlobalEnv)
  rm(ncol_net,dummy)
}


all_genes_network=union(dummy_net[,1],dummy_net[,2])
all_genes_data=rownames(dummy_exp)
d_genes=intersect(all_genes_data,all_genes_network)

## Assign

d_enrich=overlap_regulon_gene_set(dummy_net,
rownames(dummy_exp),gene_set[,1],n_name,plotx)

# if(exists(paste0("Regulon_enrichment",n_name))==FALSE)
# {
# d_enrich=overlap_regulon_gene_set(dummy_net,
#  rownames(dummy_exp),gene_set[,1],n_name,plotx)
# }
#
# else
# {
#   d_enrich=get(paste0("Regulon_enrichment",n_name))
# }

#assign(x = paste0("Regulon_enrichment",n_name),value=d_enrich,
#       envir = .GlobalEnv)

##Assign
require(mixtools)

if(exists(paste0("regul_",n_name))==FALSE)
{
d_regul=aracne2reguloneren(dummy_net[,1:3],
                           dummy_exp,format="3col",gene=FALSE)
assign(x = paste0("regul_",n_name),value=d_regul,envir = .GlobalEnv)
}

else
{
d_regul=get(paste0("regul_",n_name))
}
#assign(x = paste0("regul_",n_name),value=d_regul,envir = .GlobalEnv)

## Assign

if(exists(paste0("viper_",n_name))==FALSE)
{
d_viper=run_viper(dummy_exp[d_genes,],
                  dummy_exp[d_genes,],
                  d_pheno,dummy_net,
                  "eren",20,1,d_regul)
assign(x = paste0("viper_",n_name),value=d_viper,envir = .GlobalEnv)
}
else
{
d_viper=get(paste0("viper_",n_name))
}

#assign(x = paste0("viper_",n_name),value=d_viper,envir = .GlobalEnv)


xx=summary(d_viper[[1]],length(d_viper[[1]]$es$nes))
xy=xx$FDR
names(xy)=xx$Regulon

set_status=xy
gset_status[names(gset_status)]=-5
gset_status[intersect(names(gset_status),asthma_gene_panel[,1])]=3


if(plotx==TRUE)
{
plot_viper_eren_v1(d_viper[[1]],22,pval_name="FDR",
                   pval_n=xy,geneset_info = gset_status,
                   gset_name=paste0("Known Asthma",n_name),
                   geneset_info1 = gset_info1)

}

#View(d_enrich)



if(dolasso=="TRUE")
{
xx=summary(d_viper[[1]],length(names(d_viper[[1]]$es$nes)))
tfs=as.character(xx$Regulon)
weights1=matrix(0,length(xx$FDR),1)
rownames(weights1)=as.character(xx$Regulon)
weights1[,1]=1/(abs(xx$NES+0.01)) #v0 1/(abs(xx$NES+0.01)) #v1 1-xx$FDR
weights2=matrix(0,length(xx$FDR),1)
rownames(weights2)=as.character(xx$Regulon)
weights2[,1]= qnorm((10^(-d_enrich[as.character(xx$Regulon),1]) + 0.99999)/2)/sqrt(2)
#v0 qnorm((10^(-d_enrich[as.character(xx$Regulon),1]) + 0.99999)/2)/sqrt(2)
#v1 1-10^(-d_enrich[as.character(xx$Regulon),1])


####################
##########
#v1 include
#v1 qnorm((10^(-d_enrich[as.character(xx$Regulon),1]) + 0.99999)/2)/sqrt(2)
#10^{-d_enrich[as.character(xx$Regulon),1]}
#*10^(-d_enrich[as.character(xx$Regulon),1])
weights_comb=matrix(0,length(xx$FDR),1)
rownames(weights_comb)=as.character(xx$Regulon)
#weights_comb[,1]=weights1[,1] * (weights2[,1])

##New Weight
alpha=0.5

##v0
#weights_comb[,1]=(weights1[,1]) * (weights2[,1])

##v1

weights_comb[,1]=as.matrix((weights1[,1]) * (weights2[,1]))

idx=which(is.infinite(weights_comb[,1])==TRUE)

weights_comb[idx,1]=max(weights_comb[-idx,1])

#######

tfs= as.character(xx$Regulon)
idx=which(dummy_net[,1] %in% tfs)
d_network=dummy_net[idx,]

#print(length(tfs))
## Assign
dummy_coverage_weighted=gene_set_coverage_network_weighted(gene_set[,1],
rownames(dummy_exp),d_network,paste0("Asthma Panel, Network= ",n_name),
weights_comb,tfs,plotx)
}
#assign(x = paste0("Regulon_enrichment",n_name),value=d_enrich,
#       envir = .GlobalEnv)

#Assign
df=matrix(0,length(xx$FDR),6)
rownames(df)=as.character(xx$Regulon)
colnames(df)=1:6
if(dolasso==TRUE)
{
df[,1:3]=dummy_coverage_weighted[[1]]
colnames(df)[1:3]=colnames(dummy_coverage_weighted[[1]])
}
df[,4]=-log10(xx$FDR)
colnames(df)[4]="Viper_FDR"
df[,5]=d_enrich[as.character(xx$Regulon),1]
colnames(df)[5]="Fisher_Regulon_FDR"
idx=which(df[,5]==0)
df[idx,5]=min(df[-idx,5])
df[,6]=df[,4]*df[,5]
colnames(df)[6]="FDR_Mult"
#assign(x = paste0("Enrichment_total_",n_name),value=df,envir = .GlobalEnv)
return(df)
}

gene_set_coverage_network_weighted <- function(gene_set,all_genes,
  network,figure_name,weights,tfs,plotx)
{
  require(CVXfromR)
  setup.dir <- "C:/Users/mehmeteren/Desktop/cvx-w64/cvx"



  #tfs=unique(network[,1])
  genes_network=union(tfs,unique(network[,2]))
  all_genes=intersect(all_genes,genes_network)
  genes_cluster=gene_set#intersect(gene_set,all_genes)
  genes_set_network=intersect(gene_set,all_genes)

  d_adj_matrix=matrix(0,length(genes_set_network),length(tfs))
  rownames(d_adj_matrix)=genes_set_network
  colnames(d_adj_matrix)=tfs

  kk=1
  View(weights)
  for(i in genes_set_network)
  {
    d_neigh=network[which(network[,2]==i),1]
    d_adj_matrix[i,intersect(d_neigh,tfs)]= (weights[intersect(d_neigh,tfs),1]) #v0  (weights[intersect(d_neigh,tfs),1]) # v1 1
    kk=kk+1
    }

  lasso=list()
  percent_covered=vector()

  idx=which(rowSums(d_adj_matrix)!=0)



  adj_matrix=d_adj_matrix[idx,]
  #print(dim(adj_matrix))

  ### CVX part
  #norm(w'*b, 1)
  p=length(tfs)

  View(weights)
  # V1
  # cvxcode <- paste("variables b(p)",
  #                  "minimize(alpha*(nsett/p)*norm(w'*b,1)+(1-alpha)*norm(1-adj_matrix*b,1))",
  #                  "subject to",
  #                  "adj_matrix*b>=1",
  #                  "0<=b",
  #                  #"(0.8*nsett)<e_p*adj_matrix*b",
  #                  sep=";")
  # lasso=CallCVX(cvxcode, const.vars=list(adj_matrix=adj_matrix,p=p,
  #                                        e_p=matrix(1,1,length(genes_set_network)),
  #                                        w=(1/(weights[colnames(adj_matrix),]+0.0001)),
  #                                        alpha=0.5,
  #                                        nsett=length(genes_set_network)),
  #               opt.var.names="b", setup.dir=setup.dir)




  ### v0
  # # p=length(tfs)
  cvxcode <- paste("variables b(p)",
                   "minimize(norm(w'*b, 1))","subject to","1<adj_matrix*b",
                   "0<=b",
                   sep=";")
  lasso=CallCVX(cvxcode, const.vars=list(adj_matrix=adj_matrix,
                                         p=p,w=weights),
                opt.var.names="b", setup.dir=setup.dir)
  #
  ### v2
  # p=length(tfs)
   # cvxcode <- paste("variables b(p)",
   #                 "minimize(norm(b, 2))","subject to",#"1<adj_matrix*b",
   #                "0<=b",
   #                "adj_matrix*b>=1",
   #               "(0.8*nsett)<e_p*adj_matrix*b",
   #                  sep=";")
   # lasso=CallCVX(cvxcode, const.vars=list(adj_matrix=adj_matrix,
   #                                        p=p,w=weights,
   #                                        e_p=matrix(1,1,length(genes_set_network)),
   #                                        nsett=length(genes_set_network)),
   #               opt.var.names="b", setup.dir=setup.dir)

  #names(lasso)=figure_name
  #View(as.matrix(lasso$b))


  ##Get the coverage with lasso solution
  d_order=order(abs(lasso$b),decreasing=TRUE)
  n_neightbours=vector()
  for(kl in 1:length(tfs))
  {
    tfss=colnames(d_adj_matrix)[d_order[1:kl]]
    neighbours=network[which(network[,1] %in% tfss  ),2]
    neighbours11=network[which(network[,1] %in% colnames(d_adj_matrix)[d_order[kl]] ),2]
    percent_covered[kl]=100*(length(intersect(neighbours,genes_cluster)))/length(genes_cluster)
    n_neightbours[kl]=length(intersect(neighbours11,genes_cluster))
    rm(tfss,neighbours,neighbours11)
  }



  lasso_coverage_matrix=matrix(0,length(tfs),3)
  rownames(lasso_coverage_matrix)=colnames(d_adj_matrix)
  colnames(lasso_coverage_matrix)=c("Lasso Weight","Total_Percent_Explained",
                                    "Coverage")


  lasso_coverage_matrix[colnames(d_adj_matrix)[d_order],1]=abs(lasso$b)[d_order]
  lasso_coverage_matrix[colnames(d_adj_matrix)[d_order],2]=percent_covered
  lasso_coverage_matrix[colnames(d_adj_matrix)[d_order],3]=n_neightbours

  idx=order(n_neightbours,decreasing=TRUE)

  # if(plotx==TRUE)
  # {
  #
  #   require(gridExtra)
  #
  #
  # # plot1=qplot(1:length(tfs),n_neightbours[idx])+
  # #   geom_line()+
  # #   xlab("TFs")+
  # #   ylab("Percent Explained of the gene set")+
  # #   ggtitle(paste0("TF Coverage (Gene Set= ",figure_name," )"))
  #
  #
  # plot1=qplot(1:length(tfs),percent_covered,geom="line")+
  #   xlab("TFs")+
  #   ylab("Percent Explained of the gene set")+
  #   ggtitle(paste0("Percent Covered  Lasso (Gene Set= ",figure_name," )"))+
  #   theme_classic()
  #
  #
  # }

  enrichment_matrix=matrix(0,length(tfs),1)
  rownames(enrichment_matrix)=tfs
  enrichment_matrix[tfs[d_order],1]=n_neightbours

  idx=order(enrichment_matrix[,1],decreasing=TRUE)

  # d_order=order(abs(lasso$b),decreasing=TRUE)
  n_neightbours=vector()
  percent_covered=vector()
  for(kl in 1:length(tfs))
  {
    tfss=tfs[idx[1:kl]]
    neighbours=network[which(network[,1] %in% tfss  ),2]
    neighbours11=network[which(network[,1] %in% colnames(d_adj_matrix)[idx[kl]] ),2]
    percent_covered[kl]=100*(length(intersect(neighbours,genes_cluster)))/length(genes_cluster)
    n_neightbours[kl]=length(intersect(neighbours11,genes_cluster))
    rm(tfss,neighbours,neighbours11)
  }

  # if(plotx==TRUE)
  # {
  #
  #   plot2=qplot(1:length(tfs),percent_covered,geom="line")+
  #     xlab("TFs")+
  #     theme_classic()+
  #     ylab("Percent Explained of the gene set")+
  #     ggtitle(paste0("Percent Covered Ordering (Gene Set= ",figure_name," )"))
  #
  # #plot(1:length(tfs),percent_covered,type="l",xlab="TFs",ylab="Percent Explained of the gene set")
  # #title(paste0("Percent Covered Ordering (Gene Set= ",figure_name," )"))
  # }




  order_coverage_matrix=matrix(0,length(tfs),2)
  rownames(order_coverage_matrix)=tfs
  colnames(order_coverage_matrix)=c("Total_Percent_Explained",
                                    "Coverage")


  order_coverage_matrix[tfs[idx],2]=n_neightbours
  order_coverage_matrix[tfs[idx],1]=percent_covered

  aa=list(lasso_coverage_matrix,order_coverage_matrix)
  names(aa)=c("lasso","order")
  xxy=data.frame(lasso_coverage_matrix,order_coverage_matrix)

  idx1=order(xxy[,2],decreasing = FALSE)
  idx2=order(xxy[,4],decreasing = FALSE)

  coverage_df=data.frame(tfs=1:length(idx1),lasso=0.01*length(genes_cluster)*as.numeric(xxy[idx1,2]),
                         order=0.01*length(genes_cluster)*as.numeric(xxy[idx2,4]),
                         stringsAsFactors = FALSE)
  View(coverage_df)

  pp=ggplot()+
    geom_line(data=coverage_df,aes(tfs,lasso,color="red"),size=1.5)+
    geom_line(data=coverage_df,aes(tfs,order,color="blue"),size=1.5)+
    xlab("TFs ordered by coverage weights")+
    ylab("Number of asthma panel genes covered")+
    scale_color_manual("Method",values=c("red","blue"),label=c("LASSO","Greedy"))+
    theme_classic()+
    theme(text = element_text(size=20))+
    geom_vline(xintercept = 2,linetype="dashed",aes(size=1.5))+
    geom_hline(yintercept = 24,linetype="dashed",aes(size=1.5))+
    scale_x_continuous(breaks = round(seq(0, 132, by = 20),1))
  plot(pp)

  return(aa)
}

overlap_regulon_gene_set <- function(network,all_genes,gene_set,nname,plotx)
{

all_genes_network=union(unique(network[,1]),unique(network[,2]))
all_genes=intersect(all_genes,all_genes_network)

TFs=unique(network[,1])
nTFs=length(TFs)
Enrichment_matrix=matrix(1,nTFs,5)
colnames(Enrichment_matrix)=c("FDR","Regulon Size","Gene Set Size","Overlap Size","pval")
rownames(Enrichment_matrix)=TFs
gene_set=intersect(gene_set,all_genes)
Enrichment_matrix[,3]=length(gene_set)
for(i in TFs)
{
idx=which(network[,1] %in% i)
regul=network[idx,2]
regul=intersect(regul,all_genes)
Enrichment_matrix[i,2]=length(regul)
Enrichment_matrix[i,4]=length(intersect(regul,gene_set))


###
### For manuscript figure
# if(i=="ETV4")
# {
# print(i)
# }



if(length(intersect(regul,gene_set))>1)
{
Enrichment_matrix[i,]=fisher_test_enrichment(gene_set,regul,all_genes)[1]
}
}


Enrichment_matrix[,5]=Enrichment_matrix[,1]
Enrichment_matrix[,1]=p.adjust(Enrichment_matrix[,1],method="BH")
Enrichment_matrix[,1]=-log10(Enrichment_matrix[,1])
##Convert to FDR
#View(Enrichment_matrix)
## Plot


require(ggplot2)
require(reshape2)

#names(dat) <- paste("X", 1:10)
#dat2 <- melt(dat, id.var = "X1")
nth=0.5
idx1=order(Enrichment_matrix[,1],decreasing = TRUE)
idx1=idx1[1:10]



Enrichment_matrix1=Enrichment_matrix[idx1,]

idx=order(Enrichment_matrix1[,1],decreasing=FALSE)

dat<-data.frame(name=rownames(Enrichment_matrix1)[idx],
                value=Enrichment_matrix1[idx,1])


if(plotx==TRUE)
{
p=ggplot(dat, aes(1,factor(name,levels=rownames(Enrichment_matrix1)[idx]))) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(fill = dat$value, label = round(dat$value, 4),angle=0)) +
  scale_fill_gradient(low = "white", high = "red")+
  xlab("-log10(FDR)")+
  ylab("Samples")+
  ggtitle(paste0("Enrichment of Regulons, Network=",nname))


  plot(p)
}
return(Enrichment_matrix)
}

aracne2reguloneren <- function(afile, eset, gene = FALSE, format=c("adj", "3col"), verbose = TRUE) {
  format <- match.arg(format)
  if (is(eset, "ExpressionSet")) eset <- exprs(eset)
  if (verbose) message("\nLoading the dataset...")
  if (length(eset)==1) {
    tmp <- strsplit(readLines(eset), "\t")
    dset <- t(sapply(tmp[-1], function(x) as.numeric(x[-(1:2)])))
    colnames(dset) <- tmp[[1]][-(1:2)]
    rownames(dset) <- sapply(tmp[-1], function(x) x[1])
    annot <- t(sapply(tmp[-1], function(x) x[1:2]))
  }
  else {
    dset <- eset
    annot <- rownames(eset)
    names(annot) <- rownames(eset)
    rm(eset)
  }
  #Collapsing interactomes
  switch(format,
         adj={
           aracne <- readAracneAdj(afile)
         },
         "3col"={
           tmp <- afile#t(sapply(strsplit(readLines(afile), "\t"), function(x) x[1:3]))
           aracne <- data.frame(tf=tmp[, 1], target=tmp[, 2], mi=as.numeric(tmp[, 3])/max(as.numeric(tmp[, 3])))
         })
  if (gene) {
    if (verbose) message("Collapsing the interactomes to the gene level...")
    tmp <- aracne[order(aracne$mi, decreasing=TRUE), ]
    tmp$tf <- annot[match(tmp$tf, annot[, 1]), 2]
    tmp$target <- annot[match(tmp$target, annot[, 1]), 2]
    aracne <- tmp[!duplicated(paste(tmp$tf, tmp$target, sep="_")), ]
    #Generating the gene centric datasets
    rownames(dset) <- annot[match(rownames(dset), annot[, 1]), 2]
    dset <- filterCV(dset)
  }
  if (verbose) message("Generating the regulon objects...")
  tmp <- aracne[!is.na(aracne$mi), ]
  tmp <- tmp[rowSums(matrix(as.matrix(tmp[, 1:2]) %in% rownames(dset), nrow(tmp), 2))==2, ]
  aracne <- tapply(1:nrow(tmp), tmp$tf, function(pos, tmp) {
    tfmode <- rep(0, length(pos))
    names(tfmode) <- tmp$target[pos]
    list(tfmode=tfmode, likelihood=tmp$mi[pos])
  }, tmp=tmp)
  names(aracne) <- levels(tmp$tf)
  aracne <- TFmode1(aracne, dset)
  rm(dset)
  # removing missing data from the aracne regulon
  aracne <- aracne[names(aracne) != "NA"]
  aracne <- lapply(aracne, function(x) {
    filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood))
    x$tfmode <- x$tfmode[filtro]
    x$likelihood <- x$likelihood[filtro]
    return(x)
  })
  aracne <- aracne[sapply(aracne, function(x) length(names(x$tfmode)))>0]
  regul <- TFscore(aracne, verbose=verbose)
  class(regul) <- "regulon"
  return(regul)
}

run_viper <- function(exp_data_net,exp_data_viper,phenotype,
                    network,net_name,threshold,cores,regul,
                    minsize1=25,minsizeAREA=20)
{
## exp_data row gene columns samples
## phenotype 2nd row is status
## network is in 3D format
#require(viper)
require(mixtools)
require(Biobase)
  require(BiocGenerics)
require(parallel)

source("C:/Users/mehmeteren/OneDrive/Rcode_generic/aracne2reguloneren.R")
source("C:/Users/mehmeteren/OneDrive/Rcode_generic/internal.R")
  source("C:/Users/mehmeteren/OneDrive/Rcode_generic/viper/msviper.R")
  source("C:/Users/mehmeteren/OneDrive/Rcode_generic/viper/viper.R")
  source("C:/Users/mehmeteren/OneDrive/Rcode_generic/viper/ledge.R")
  source("C:/Users/mehmeteren/OneDrive/Rcode_generic/viper/shadow.R")
  source("C:/Users/mehmeteren/OneDrive/Rcode_generic/viper/general.R")
#  source("C:/Users/mehmeteren/OneDrive/Rcode_generic/viper/aracne.R")
#network=network[,1:3]


#regul=aracne2reguloneren(network,exp_data_net,format="3col")
phenoData <- new("AnnotatedDataFrame", data=data.frame(phenotype))
eset <- ExpressionSet(assayData=as.matrix(exp_data_viper), phenoData=phenoData)
name_pheno=colnames(phenotype)[2]

signature <- rowTtest(eset, name_pheno, "A", "N")
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
nullmodel <- ttestNull(eset, name_pheno, "A", "N", per = 1000, repos = TRUE, verbose = FALSE)
print("eren")
mrs <- msviper(ges=signature, regulon=regul,
               nullmodel=nullmodel, verbose = FALSE,
               minsize=minsize1)
mrs <- ledge(mrs)

print(minsize1)
#mrshadow <- shadow(mrs,regulators=threshold, verbose = FALSE) ## regulators signify number of regulators
mrshadow=mrs
print("eren1")
#mrssynergy <- msviperCombinatorial(mrs, regulators = threshold,
#                                   verbose = FALSE,cores = cores)
mrssynergy=mrs
print("eren2")
#mrssynergy <- msviperSynergy(mrssynergy, verbose = FALSE)
print("eren3")
aa=list(mrs,mrshadow,mrssynergy)
#aa=list(mrs,mrshadow)
return(aa)

}

plot_viper_eren_v1 <- function(x, mrs=10, color=c("cornflowerblue","salmon"), pval=NULL,
         bins=500,pval_name="p-value" ,pval_n=NULL,
         cex=0, density=0, smooth=0, sep=.2, hybrid=TRUE,
         include=c("expression", "activity"), gama=2,gset_name="Gset Info",
         geneset_info=NULL,geneset_info1=NULL,...) {
  maobject <- x


  if (is.null(pval_n))
  {
    pval_n<-maobject$es$p.value

  }

  rm(x)
  marg <- par("mai")
  rlist <- maobject$signature
  if (ncol(rlist)>0) rlist <- rowMeans(rlist)
  if (length(mrs)==1 & is.numeric(mrs[1])) {
    mrs <- names(maobject$es$nes)[order(maobject$es$p.value)[1:round(mrs)]]
    mrs <- mrs[order(maobject$es$nes[match(mrs, names(maobject$es$nes))], decreasing=TRUE)]
  }
  groups <- maobject$regulon[match(mrs, names(maobject$regulon))]
  if (is.null(pval)) pval <- pval_n[match(mrs, names(pval_n))]#maobject$es$p.value[match(mrs, names(maobject$es$nes))]
  if (is.data.frame(pval))
  {
    pval <- as.matrix(pval)[match(mrs, names(maobject$es$nes)),1]
  }

  if (is.matrix(rlist)) rlist <- rlist[,1]
  if (min(rlist)<0) rlist <- sort(rlist)
  else rlist <- sort(-rlist)
  groups <- groups[length(groups):1]
  color1 <- color
  layout(matrix(1:2, 1, 2), widths=c(10-length(include), length(include)))
  color <- rgb2hsv(col2rgb(color))
  satval <- color[3,]
  color <- color[1,]
  preg <- as.numeric(any(sapply(groups, function(x) any(x$tfmode<0))))+1
  textsize <- 1
  xlimit <- c(0, length(rlist)*(1+.2*max(nchar(names(groups)))/8))
  if (!is.null(pval)) {
    if (is.matrix(pval)) {
      pval <- pval[nrow(pval):1, ]
      pval <- pval[, ncol(pval):1]
      xlimit <- c(-length(rlist)*(sep*ncol(pval)+.02), length(rlist)*(1+.2*max(nchar(names(groups)))/9))
    }
    else {
      pval <- pval[length(pval):1]
      xlimit <- c(-length(rlist)*.12, length(rlist)*(1+.2*max(nchar(names(groups)))/9))
      textsize=.8
    }
  }
  if (length(include)>0 & include[1] != "") par(mai=c(marg[1:3], .05))
  plot(1, type="n", ylim=c(0, length(groups)), xlim=xlimit, axes=FALSE, ylab="", xlab="", yaxs="i")
  if (cex>0) textsize <- (length(groups)<=20)*cex+(length(groups)>20)*(20/length(groups))*cex
  switch(preg,
         {
           for (i in 1:length(groups)) {
             densi <- rep(0, length(rlist))
             x <- which(names(rlist) %in% names(groups[[i]]$tfmode))
             if (length(x)>0) {
               densi[x] <- 1
               denStep <- round(length(densi)/bins)
               x1 <- x[x<denStep]
               x2 <- x[x>=denStep & x <= (length(rlist)-denStep)]
               x3 <- x[x>(length(rlist)-denStep)]
               densiRes <- sapply(x2, function(i, densi, denStep)
                 sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
               densiRes <- densiRes/max(densiRes)
               temp <- rlist[x]
               if (satval[1]==0) temp <- hsv((temp<0)*color[1] + (temp>0)*color[2], satval, 1-densiRes)
               else temp <- hsv((sign(temp)<0)*color[1] + (sign(temp)>0)*color[2], densiRes, satval)
               for (ii in order(densiRes)) lines(c(x[ii], x[ii]),c(i-1, i), col=temp[ii])
               if (density>0) {
                 denStep <- round(length(densi)/density)
                 xpos <- seq(denStep, length(rlist)-denStep, length=density)
                 densiRes <- sapply(xpos, function(i, densi, denStep) {
                   sum(densi[(i-denStep):(i+denStep)])
                 }, densi=densi, denStep=denStep)
                 densiRes <- densiRes/max(densiRes)
                 if (smooth>0) densiRes <- smooth.spline(xpos, densiRes, spar=smooth)$y
                 lines(xpos, i+densiRes-1)
               }
             }
           }
           text(rep(length(rlist)*1.02, length(groups)), 1:length(groups)-.5, names(groups), adj=0, cex=textsize)
           if (!is.null(pval)) {
             if (is.matrix(pval)) {
               for (i in 1:ncol(pval)) {
                 text(rep(-length(rlist)*(sep*(i-1)+.02), length(groups)), 1:length(groups)-.5, pval[, i], adj=1, cex=.85*textsize)
                 text(-length(rlist)*(sep*(i-1)+.02), length(groups)+.5, colnames(pval)[i], adj=1, cex=1)
               }
             }
             else {
               text(rep(-length(rlist)*.02, length(groups)), 1:length(groups)-.5, signif(pval, 3), adj=1, cex=.85*textsize)
               text(0, length(groups)+.5, ifelse(max(pval)>1, "oddsR", pval_name), adj=1, cex=1.2)
             }
           }
           text(length(rlist)*1.02, length(groups)+.5, "Set", adj=0, cex=1.2)
         },
         {
           for (i in 1:length(groups)) {
             for (ii in 1:2) {
               tset <- groups[[i]]$tfmode
               tset <- tset[(tset<0 & ii==1)|(!(tset<0) & ii==2)]
               tset1 <- names(tset)[abs(tset)>.5]
               tset2 <- names(tset)[abs(tset)<.5]
               if (length(tset)>1) {
                 densi <- rep(0, length(rlist))
                 x <- match(names(tset), names(rlist))
                 tw1 <- rep(1, length(x))
                 if (hybrid) {
                   x <- match(tset1, names(rlist))
                   if (ii==1) {
                     x <- c(x, match(tset2, names(sort(-abs(rlist)*sign(maobject$es$nes[names(maobject$es$nes) == names(groups)[i]])))))
                   }
                   else {
                     x <- c(x, match(tset2, names(sort(abs(rlist)*sign(maobject$es$nes[names(maobject$es$nes) == names(groups)[i]])))))
                   }
                   tw1 <- groups[[i]]$likelihood[match(c(tset1, tset2), names(groups[[i]]$tfmode))]
                   tw1 <- tw1/max(tw1)
                 }
                 densi[x] <- 1
                 denStep <- round(length(densi)/bins)
                 x1 <- x[x<denStep]
                 x2 <- x[x>=denStep & x <= (length(rlist)-denStep)]
                 x3 <- x[x>(length(rlist)-denStep)]
                 densiRes <- sapply(x2, function(i, densi, denStep) {
                   sum(densi[(i-denStep):(i+denStep)])
                 }, densi=densi, denStep=denStep)
                 densiRes <- densiRes*(tw1[x>=denStep & x <= (length(rlist)-denStep)])
                 densiRes <- densiRes/max(densiRes)
                 temp <- rlist[x]
                 temp <- hsv(color[ii], densiRes, satval)
                 for (iii in order(densiRes)) {
                   lines(c(x[iii], x[iii]),c(i-1+(ii-1)/2, i-1+ii/2), col=temp[iii])
                 }
                 if (density>0) {
                   denStep <- round(length(densi)/density)
                   xpos <- seq(denStep, length(rlist)-denStep, length=density)
                   densiRes <- sapply(xpos, function(i, densi, denStep)
                     sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
                   densiRes <- densiRes/max(densiRes)
                   if (smooth>0) densiRes <- smooth.spline(xpos, densiRes, spar=smooth)$y
                   lines(xpos, i-1+densiRes/2+(ii-1)/2)
                 }
               }
             }
           }
           text(rep(length(rlist)*1.02, length(groups)), 1:length(groups)-.5, names(groups), adj=0, cex=textsize)
           if (!is.null(pval)) {
             if (is.matrix(pval)) {
               for (i in 1:ncol(pval)) {
                 text(rep(-length(rlist)*(sep*(i-1)+.02), length(groups)), 1:length(groups)-.5, pval[, i], adj=1, cex=.85*textsize)
                 text(-length(rlist)*(sep*(i-1)+.02), length(groups)+.5, colnames(pval)[i], adj=1, cex=1)
               }
             }
             else {
               text(rep(-length(rlist)*.02, length(groups)), 1:length(groups)-.5, signif(pval, 3), adj=1, cex=.85*textsize)
               axis(3, -.05*length(rlist), ifelse(max(pval)>1, "oddsR", pval_name), adj=1, cex=1.2, tick=FALSE, line=-.5)
             }
           }


           axis(3, length(rlist)*1.05, "Set", adj=0, cex=1.2, line=-.5, tick=FALSE)
         })
  abline(h=0:length(groups))
  lines(c(0, 0), c(0, length(groups)))
  lines(c(length(rlist), length(rlist)), c( 0, length(groups)))
  ss <- maobject$signature
  if (!is.null(dim(ss))) ss <- ss[, 1]
  x <- NULL
  xn <- NULL
  xpos <- NULL
  if("activity" %in% include) {
    x <- maobject$es$nes[match(mrs, names(maobject$es$nes))]/max(abs(maobject$es$nes))
    xn <- "Act"
  }
  if ("expression" %in% include) {
    x <- cbind(x, ss[match(mrs, names(ss))]/max(abs(ss)))
    xn <- c(xn, "Exp")
    xpos <- rank(-abs(ss))[match(mrs, names(ss))]
  }

  if (!is.null(geneset_info)) {
    g_min=min(x[,1])
    g_max=max(x[,1])
    geneset_info[geneset_info<0]=g_min
    geneset_info[geneset_info>0]=g_max
    x <- cbind(x, geneset_info[match(mrs, names(geneset_info))])
    xn <- c(xn, gset_name)

  }

  # if (!is.null(geneset_info)) {
  #   x <- cbind(x, ss[match(mrs, names(ss))]/max(abs(ss)))
  #   xn <- c(xn, "Exp1")
  #   xpos1 <- rank(-abs(ss))[match(mrs, names(ss))]
  #
  # }


  if (!is.null(x)) {
    if(is.null(dim(x))) dim(x) <- c(length(x), 1)

    rownames(x) <- colnames(x) <- NULL
    marg[2] <- .05
    par(mai=marg)

    scmax <- max(abs(x), na.rm=TRUE)
    x <- abs(x/scmax)^gama*sign(x)
    x <- filterRowMatrix(x, nrow(x):1)
    x1 <- x
    x1[is.na(x1)] <- 0
    coli <- hsv(ifelse(x1<0, color[1], color[2]), abs(x1), 1)
    coli[is.na(x)] <- hsv(0, 0, .5)
    nxx=ncol(x)

    dy=list(nxx,x,coli)
    image(1:nxx, 1:nrow(x), t(matrix(1:(nxx*nrow(x)),
                                     nrow(x), nxx)), col=coli, ylab="", xlab="", axes=FALSE, yaxs="i")
    box()

    grid(nxx, nrow(x), col="black", lty=1)
    axis(3, 1:length(xn), xn, las=2, line=-.5, tick=FALSE)
    axis(4, length(xpos):1,xpos , las=1, tick=FALSE, line=-.4, cex.axis=.85*textsize)

    if (!is.null(geneset_info1)) {

      axis(2, length(xpos):1,geneset_info1[match(mrs, names(geneset_info1))] , las=1, tick=FALSE, line=-.4, cex.axis=.85*textsize)
    }

  }
  par(mai=marg)

 return(dy)
}
