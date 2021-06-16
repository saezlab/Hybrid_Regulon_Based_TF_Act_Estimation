# ARACNE regulons inference is based on the calculation of the mutual information (MI) between each pair of genes. 
# MI computaion has a high computational complexity, therefore, we recommend to run the script on a cluster.

library(viper)
library(dorothea)
library(WGCNA)
library(igraph)
library(doParallel)
WGCNA::enableWGCNAThreads()

# A fucntion to extract matching rows between two data frames:
check.match= function(dataframe1, columns1, dataframe2, columns2) {
  
  df1_interactions = paste(dataframe1[[columns1[1]]], dataframe1[[columns1[2]]])
  df2_interactions = paste(dataframe2[[columns2[1]]], dataframe2[[columns2[2]]])
  df2_rev_interactions = paste(dataframe2[[columns2[2]]], dataframe2[[columns2[1]]])
  
  match_direction1=match(df1_interactions, df2_interactions,nomatch = 0)
  to_skip=which(match_direction1 !=0)
  match_direction2=match(df1_interactions[-c(to_skip)], df2_rev_interactions, nomatch = 0)

  
  match_direction1[which(match_direction1==0)]=match_direction2
  match_df1=which(match_direction1!=0)
  match_df2=match_direction1[which(match_direction1!=0)]
 
  
  cat(length(match_df1),"matching rows","\n","$df1_rows: matching rows from", substitute(dataframe1),"\n", "$df2_rows: matching rows from", substitute(dataframe2),"\n")
  return(list("df1_rows"= match_df1, "df2_rows"=match_df2))
}


# Load the expression set (already normalized):
norm_counts=as.data.frame(t(readRDS(url("https://zenodo.org/record/4322914/files/dorothea_bench_expr.rds?download=1"))))

# Compute the MI between gene pairs:
mim=build.mim(norm_counts, estimator = "mi.shrink", disc = "equalwidth", nbins = sqrt(nrow(norm_counts))) # Continuous counts must be discretized

# Run ARACNE algorithm (eps=1: keeps all the indirect gene-gene interactions):
aracne_mat=aracne(mim, eps=1)
aracne_mat[is.na(aracne_mat)]=0
sim_matrix=aracne_mat

# Transform the MI values to make them range from 0 to 1:
sim_matrix=sim_matrix/max(sim_matrix)

# Soft thresholds (ST) to test:
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Pick the optimal soft threshold (the lowest ST resulting in a high R^2):
threshold_test = pickSoftThreshold.fromSimilarity(sim_matrix, powerVector = powers, verbose = FALSE)

# Rise the MI to the ST:
adj_mat=adjacency.fromSimilarity(sim_matrix,type = "unsigned",power=6)

# Apply the Topological Overlap Measure to remove spurious interactions:

TOM = TOMsimilarity(adj_mat)

# Compute the distances between the genes and cluster them into modules:
dissTOM=1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = 30)

# Create the Gene-Module dataframe:
gene_mod_df=as.data.frame(cbind(colnames(sim_matrix), dynamicMods))
rownames(gene_mod_df)=NULL
colnames(gene_mod_df)=c("Gene", "Module")

# Extract the gene-gene interactions within each module:
weighted_inter_df=data.frame()
mod_vect=NULL
for (i in 1:length(table(dynamicMods))) {
  genes_ids=gene_mod_df[gene_mod_df$Module==i,1]
  sub_sim_mat=sim_matrix[c(genes_ids),c(genes_ids)]
  pairwise_interactions=get.data.frame(graph.adjacency(sub_sim_mat, mode = 'upper', weighted = TRUE, diag=FALSE))
  mod_id=rep(i, dim(pairwise_interactions)[1])
  mod_vect=c(mod_vect,mod_id)
  weighted_inter_df=rbind(weighted_inter_df, pairwise_interactions)}
  weighted_inter_df=cbind(weighted_inter_df, mod_vect)

# Create the hybrid_dorothea_A:
dorothea_a=as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A"), c(1,3)])
matches_a=check.match(weighted_inter_df,c(1,2),dorothea_a,c(1,2))
hybrid_a=as.data.frame(cbind(dorothea_a[matches_a$df2_rows,], weighted_inter_df[matches_a$df1_rows,3]))
rownames(hybrid_a)=NULL
colnames(hybrid_a)=c("tf", "target", "mi")
write.table(hybrid_a, "hybrid_a.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
hybrid_a=aracne2regulon("hybrid_a.txt", as.matrix(norm_counts), format="3col")


# Create the hybrid_dorothea_ABC:
dorothea_regulons=as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A", "B", "C"), c(1,3)])
matches=check.match(weighted_inter_df,c(1,2),dorothea_regulons,c(1,2))
hybrid_abc=as.data.frame(cbind(dorothea_regulons[matches$df2_rows,], weighted_inter_df[matches$df1_rows,3]))
rownames(hybrid_abc)=NULL
colnames(hybrid_abc)=c("tf", "target", "mi")
write.table(hybrid_abc, "hybrid_abc.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
hybrid_abc=aracne2regulon("hybrid_abc.txt", as.matrix(norm_counts), format="3col")


# Create the hybrid_dorothea_ABCDE:
dorothea_all=as.data.frame(dorothea_hs[, c(1,3)])
matches_all=check.match(weighted_inter_df,c(1,2),dorothea_all,c(1,2))
hybrid_abcde=as.data.frame(cbind(dorothea_all[matches_all$df2_rows,], weighted_inter_df[matches_all$df1_rows,3]))
rownames(hybrid_abcde)=NULL
colnames(hybrid_abcde)=c("tf", "target", "mi")
write.table(hybrid_abcde, "hybrid_abcde.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
hybrid_abcde=aracne2regulon("hybrid_abcde.txt", as.matrix(norm_counts), format="3col")



