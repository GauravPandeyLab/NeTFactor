run_netfactor <-function(path_to_expression,n_name,path_to_network,path_to_pheno,
                         path_to_biomarker,plotx,dolasso)
{
  # Load Required Libraries
  library(viper)
  library(CVXR)
  library(mixtools)
  require(Biobase)
  require(BiocGenerics)
  require(parallel)

  #Load Data
  expression_data=read.table(path_to_expression,header = TRUE,row.names = 1,stringsAsFactors = FALSE)
  network=read.table(path_to_network,header=TRUE,stringsAsFactors = FALSE)
  pheno=read.table(path_to_pheno,row.names = 1,header=TRUE,stringsAsFactors = FALSE)
  gene_set=read.table(path_to_biomarker,header = TRUE,row.names = NULL,stringsAsFactors = FALSE)
  expression_data=as.matrix(expression_data)



  ### Choose genes that is in the intersection of network and
  ## expression data

  filtered_data=filter_data(network,expression_data)
  network=filtered_data[[1]]
  expression_data=filtered_data[[2]]
  rm(filtered_data)




  ## Run fisher enrichment
  fisher_enrich=overlap_regulon_gene_set(network,
                                         rownames(expression_data),gene_set[,1],n_name,plotx)


  # Create a regulon object for viper

  regul=networktoregulon(network[,1:3],expression_data,format="3col",gene=FALSE)


  # Run viper

  viper=run_viper(expression_data,pheno,network,regul)


  ## Create network weights

  viper_results=summary(viper,length(names(viper$es$nes)))
  tfs=as.character(viper_results$Regulon)
  weights_viper=matrix(0,length(viper_results$FDR),1)
  rownames(weights_viper)=as.character(viper_results$Regulon)
  weights_viper[,1]=(1-viper_results$FDR)
  weights_fisher=matrix(0,length(viper_results$FDR),1)
  rownames(weights_fisher)=as.character(viper_results$Regulon)
  weights_fisher[,1]= (1-fisher_enrich[as.character(viper_results$Regulon),1])

  #combine the two weights

  weights_combo=matrix(0,length(viper_results$FDR),1)
  rownames(weights_combo)=as.character(viper_results$Regulon)
  weights_combo[,1]=abs((weights_viper[,1]) * (weights_fisher[,1]))

  ## For numerical stability add small weights to 0s
  idx=which(weights_combo[,1]==0)

  if(length(idx)>0)
  {
    weights_combo[idx,1]= 10^{-8}
  }


  ####### Retain tfs from viper

  tfs=as.character(viper_results$Regulon)
  idx=which(network[,1] %in% tfs)
  d_network=network[idx,]

  ## Run Netfactor
  netfactor_results=gene_set_coverage_network_weighted(gene_set[,1],
                                                       rownames(expression_data),d_network,paste0("Asthma Panel, Network= ",n_name),
                                                       weights_combo,tfs,plotx)

  # Create output object
  output=matrix(0,length(viper_results$FDR),5)
  rownames(output)=as.character(viper_results$Regulon)
  colnames(output)=1:5

  output[,1:3]=netfactor_results
  colnames(output)[1:3]=colnames(netfactor_results)

  output[,4]=viper_results$FDR
  colnames(output)[4]="Viper_FDR"
  output[,5]=10^{-fisher_enrich[as.character(viper_results$Regulon),1]}
  colnames(output)[5]="Fisher_Regulon_FDR"
  idx=which(output[,5]==0)
  return(output)
}

gene_set_coverage_network_weighted <- function(gene_set,all_genes,
                                               network,figure_name,weights,tfs,plotx)
{

  # Get intersection
  genes_network=union(tfs,unique(network[,2]))
  all_genes=intersect(all_genes,genes_network)
  genes_cluster=gene_set#intersect(gene_set,all_genes)
  genes_set_network=intersect(gene_set,all_genes)


  ## Create adj matrix from network file for convex optimization
  adj_matrix=matrix(0,length(genes_set_network),length(tfs))
  rownames(adj_matrix)=genes_set_network
  colnames(adj_matrix)=tfs

  kk=1
  idx_tf=which(network[,1] %in% tfs)
  for(i in genes_set_network)
  {
    idx_neigh=which(network[,2]==i)
    d_neigh=network[intersect(idx_tf,idx_neigh),1]
    adj_matrix[i,d_neigh]=1 #(weights[intersect(d_neigh,tfs),1])*(1-network[idx_neigh,4]) #v0  (weights[intersect(d_neigh,tfs),1]) # v1 1
    kk=kk+1
  }


  ### Remove genes not covered
  idx=which(rowSums(adj_matrix)>=1)
  adj_matrix=adj_matrix[idx,]


  ### Convex Optimization LASSO
  p=length(tfs)
  b=Variable(p,name="b")
  w=as.matrix((1/(weights[colnames(adj_matrix),])))
  objective <- norm1(t(w)%*%b)
  constrnt <-list(1<=adj_matrix %*% b,0 <= b  )
  prob <- Problem(Minimize(objective), constrnt)
  results_lasso <- solve(prob)#,solver="MOSEK")

  ##Get the coverage with lasso solution
  lasso_b=results_lasso$getValue(b)

  d_order=order(abs(lasso_b),decreasing=TRUE)
  n_neighbours=vector()
  percent_covered=vector()
  for(kl in 1:length(tfs))
  {
    tfss=colnames(adj_matrix)[d_order[1:kl]]
    neighbours=network[which(network[,1] %in% tfss  ),2]
    neighbours11=network[which(network[,1] %in% colnames(adj_matrix)[d_order[kl]] ),2]
    percent_covered[kl]=100*(length(intersect(neighbours,genes_cluster)))/length(genes_cluster)
    n_neighbours[kl]=length(intersect(neighbours11,genes_cluster))
    rm(tfss,neighbours,neighbours11)
  }



  lasso_coverage_matrix=matrix(0,length(tfs),3)
  rownames(lasso_coverage_matrix)=colnames(adj_matrix)
  colnames(lasso_coverage_matrix)=c("Lasso Weight","Total_Percent_Explained",
                                    "Coverage")

  rownames(lasso_coverage_matrix)=colnames(adj_matrix)
  colnames(lasso_coverage_matrix)=c("Lasso Weight","Total_Percent_Explained",
                                    "Coverage")


  lasso_coverage_matrix[colnames(adj_matrix)[d_order],1]=abs(lasso_b)[d_order]
  lasso_coverage_matrix[colnames(adj_matrix)[d_order],2]=percent_covered
  lasso_coverage_matrix[colnames(adj_matrix)[d_order],3]=n_neighbours


  return(lasso_coverage_matrix)
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
    if(length(intersect(regul,gene_set))>1)
    {
      Enrichment_matrix[i,]=fisher_test_enrichment(gene_set,regul,all_genes)[1]
    }
  }


  Enrichment_matrix[,5]=Enrichment_matrix[,1]
  Enrichment_matrix[,1]=p.adjust(Enrichment_matrix[,1],method="BH")
  Enrichment_matrix[,1]=Enrichment_matrix[,1]
  return(Enrichment_matrix)
}

networktoregulon <- function(afile, eset, gene = FALSE, format=c("adj", "3col"), verbose = TRUE) {
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
    #dset <- filterCV(dset)
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


fisher_test_enrichment<-function(gene_set,diff_set,all_genes)
{
  gene_set=intersect(gene_set,all_genes)
  diff_set=intersect(diff_set,all_genes)
  n_gene_set=setdiff(all_genes,gene_set)
  n_diff_set=setdiff(all_genes,diff_set)

  aa=matrix(0,2,2)
  aa[1,1]=length(intersect(gene_set,diff_set))
  aa[1,2]=length(intersect(gene_set,n_diff_set))
  aa[2,1]=length(intersect(n_gene_set,diff_set))
  aa[2,2]=length(intersect(n_gene_set,n_diff_set))

  xx=fisher.test(aa,alternative = "greater")#,log.p=TRUE)

  dd=c(xx$p.value,length(diff_set),length(gene_set),length(intersect(diff_set,gene_set)))
  return(dd)
}

TFmode1 <- function (regulon, expset, method = "spearman") {
  regulon <- updateRegulon(regulon)
  regulon <- regulon[names(regulon) %in% rownames(expset)]
  regulon <- lapply(regulon, function(x, genes) {
    filtro <- names(x$tfmode) %in% genes
    x$tfmode <- x$tfmode[filtro]
    if (length(x$likelihood) == length(filtro))
      x$likelihood <- x$likelihood[filtro]
    return(x)
  }, genes = rownames(expset))
  tf <- unique(names(regulon))
  tg <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))

  cmat <- cor(t(expset[rownames(expset) %in% tf, ]), t(expset[rownames(expset) %in% tg, ]), method = "spearman")
  reg <- lapply(1:length(regulon), function(i, regulon, cmat) {
    tfscore <- cmat[which(rownames(cmat) == names(regulon)[i]), match(names(regulon[[i]]$tfmode), colnames(cmat))]
    list(tfmode = tfscore, likelihood = regulon[[i]]$likelihood)
  }, regulon = regulon, cmat = cmat)
  names(reg) <- names(regulon)
  return(reg)
}



updateRegulon <- function(regul) {
  if (is.null(names(regul[[1]]))) {
    tmp <- lapply(regul, function(x) {
      tmp <- rep(0, length(x))
      names(tmp) <- x
      list(tfmode=tmp, likelihood=rep(1, length(tmp)))
    })
    return(tmp)
  }
  if (names(regul[[1]])[1]=="tfmode") return(regul)
  return(lapply(regul, function(x) list(tfmode=x, likelihood=rep(1, length(x)))))
}

TFscore <- function (regul, mu = NULL, sigma = NULL, verbose=TRUE) {
  if (length(mu) == 3 & length(sigma) == 3)
    fit <- list(mu = mu, sigma = sigma)
  else {
    tmp <- unlist(lapply(regul, function(x) x$tfmode), use.names = FALSE)
    fit <- list(mu = c(-0.5, 0, 0.5), sigma = c(0.15, 0.25, 0.15), lambda = c(0.2, 0.4, 0.4), all.loglik = rep(0, 1001))
    count <- 0
    while (length(fit$all.loglik) > 1000 & count < 3) {
      fit <- normalmixEM(tmp, mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda, mean.constr = c(NA, 0, NA), maxit = 1000, verb = FALSE)
      count <- count + 1
    }
  }
  if (verbose) message("mu: ", paste(fit$mu, collapse = ", "), ". sigma: ", paste(fit$sigma, collapse = ", "))
  regul <- lapply(regul, function(x, fit) {
    g2 <- pnorm(x$tfmode, fit$mu[3], fit$sigma[3], lower.tail = TRUE)
    g1 <- pnorm(x$tfmode, fit$mu[1], fit$sigma[1], lower.tail = FALSE)
    g0 <- pnorm(x$tfmode, fit$mu[2], fit$sigma[2], lower.tail = FALSE)
    g00 <- pnorm(x$tfmode, fit$mu[2], fit$sigma[2], lower.tail = TRUE)
    x$tfmode <- g2/(g1 + g0 + g2) * (x$tfmode >= 0) - g1/(g1 + g00 + g2) * (x$tfmode < 0)
    return(x)
  }, fit = fit)
  return(regul)
}


filter_data <- function(network,expression_data)
{
  all_genes_network=union(network[,1],network[,2])
  all_genes_data=rownames(expression_data)
  genes_intersect=intersect(all_genes_data,all_genes_network)
  idx_eset=which(rownames(expression_data) %in% genes_intersect)
  expression_data=expression_data[idx_eset,]
  idx_tfs=which(network[,1] %in% genes_intersect)
  idx_targets=which(network[,2] %in% genes_intersect)
  idx_intersect=intersect(idx_tfs,idx_targets)
  network=network[idx_intersect,]
  filtered_data=list()
  filtered_data[[1]]=network
  filtered_data[[2]]=expression_data
  return(filtered_data)
}

run_viper <- function(exp_data_viper,phenotype,
                      network,regul,
                      minsize1=25,minsizeAREA=20)
{
  phenoData <- new("AnnotatedDataFrame", data=data.frame(phenotype))
  eset <- ExpressionSet(assayData=as.matrix(exp_data_viper), phenoData=phenoData)
  name_pheno=colnames(phenotype)[1]

  signature <- rowTtest(eset, name_pheno, "A", "N")

  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  nullmodel <- ttestNull(eset, name_pheno, "A", "N", per = 1000, repos = TRUE, verbose = FALSE)
  mrs <- msviper(ges=signature, regulon=regul,
                 nullmodel=nullmodel, verbose = FALSE,
                 minsize=minsize1)
  mrs <- ledge(mrs)
  return(mrs)

}
