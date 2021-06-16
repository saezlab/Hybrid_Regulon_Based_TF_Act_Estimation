library(WGCNA)
library(igraph)
library(dorothea)
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

# Compute the Bicor correlation between gene pairs:
sim_matrix=WGCNA::bicor(norm_counts)

# Soft thresholds (ST) to test:
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Pick the optimal soft threshold (the lowest ST resulting in a high R^2):
threshold_test = pickSoftThreshold.fromSimilarity(sim_matrix, powerVector = powers, verbose = FALSE)

# Rise the Bicor cor-coefficients to the ST:
adj_mat=adjacency.fromSimilarity(sim_matrix,type = "signed",power=12)

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
hybrid_a=as.data.frame(cbind(dorothea_a[matches_a$df2_rows,], weighted_inter_df[matches_a$df1_rows,3], rep(1,length(matches_a$df1_rows))))
rownames(hybrid_a)=NULL
colnames(hybrid_a)=c("tf", "target", "mor", "likelihood")
hybrid_a$likelihood=abs(hybrid_a$mor/max(abs(hybrid_a$mor)))

# Create the hybrid_dorothea_ABC:
dorothea_regulons=as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A", "B", "C"), c(1,3)])
matches=check.match(weighted_inter_df,c(1,2),dorothea_regulons,c(1,2))
hybrid_abc=as.data.frame(cbind(dorothea_regulons[matches$df2_rows,], weighted_inter_df[matches$df1_rows,3], rep(1,length(matches$df1_rows))))
rownames(hybrid_abc)=NULL
colnames(hybrid_abc)=c("tf", "target", "mor", "likelihood")
hybrid_abc$likelihood=abs(hybrid_abc$mor/max(abs(hybrid_abc$mor)))

# Create the hybrid_dorothea_ABCDE:
dorothea_all=as.data.frame(dorothea_hs[, c(1,3)])
matches_all=check.match(weighted_inter_df,c(1,2),dorothea_all,c(1,2))
hybrid_abcde=as.data.frame(cbind(dorothea_all[matches_all$df2_rows,], weighted_inter_df[matches_all$df1_rows,3], rep(1,length(matches_all$df1_rows))))
rownames(hybrid_abcde)=NULL
colnames(hybrid_abcde)=c("tf", "target", "mor", "likelihood")
hybrid_abcde$likelihood=abs(hybrid_abcde$mor/max(abs(hybrid_abcde$mor)))
