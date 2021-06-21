# Hybrid_Regulon_Based_TF_Act_Estimation

The objective of this work is to try to improve the Transcription Factors Activity Estimation (TFAE). TFAE functions (e.g. viper()) requires two inputs: an
expression set and a gene regulatory network (regulons). Previously, we were
using DoRothEA regulons [1] containing TFs regulating their gene targets according to a given mode of regulation (MOR): 1 for upregulation and -1 for
downregulation. Here, we have substituted these discrete MOR by continuous
ones, which were inferred from the expression set using WGCNA and ARACNE.
Seemingly, using these hybrid regulons (DoRothEA regulons with an inferred
MOR) increases the performance of TFAE.
