### The SGLT2 Inhibitor Empagliflozin Reduces Circulating Monocyte Platelet Aggregates: A Pilot Study
#### Methods  
PLINK files were converted to VCF format using the plink[1] command. Lift over to GRCh38 was done with picard[2]. The UCSC GRCh38/hg38 assembly[3] was used as the target reference build. FASTQ files were processed using 10x Genomics Cell Ranger Count v7.1.0 using 10x Genomics Cloud Analysis[4]. A KIR-modified GRCh38 reference[5] was used for read alignment. Souporcell[6] was used for genetic demultiplexing of the Cell Ranger counts and annotation of inter-sample doublets. The Demuxafy[16] Assign_Indiv_by_Geno method was used to calculate the strength of correlation between variants called by genomic sequencing and those called by Souporcell. These variant call correlations were used to map sample IDs to Souporcell genetic clusters and subsequently to each barcode.
<br>
<br>
The Solo algorithm[7] as implemented by the scvi-tools[8] python library was used to identify intra-sample doublets. This was done on each Cell Ranger count output separately. Barcodes called doublet or unassigned by Souporcell were dropped for Solo. Additionally, cells with counts from less than 200 genes and genes not expressed in at least 10 cells were removed. Highly variable genes (HVG, n = 5000) were calculated and a single-cell Variational Inference (scVI) model was trained on the HVG. The trained model was used to facilitate the Solo scoring and predictions. Predicted doublets were annotated.
<br>
<br>
Using the R programming language v4.3.1[9] and Seurat v5[10], a Seurat object was constructed containing counts from all counts matrices with barcode-linked metadata. Barcodes assigned 'doublet' or 'unassigned' by Souporcell or barcodes assigned 'doublet' by Solo were dropped. The raw counts matrix, object metadata, and gene names were saved to file. To make use of the many scverse tools[11], these core elements were used to construct an anndata[12] object in Python v3.9.18[13]. 
<br>
<br>
Cells with counts from less than 200 genes or more than 10% counts from mitochondrial (MT) genes were removed. Genes not expressed in at least 10 cells were removed. HVG were selected and a scVI model with hyperparameter ray tuning was trained on the data for integration of the counts data. The selection of HVG included calculation of the 5000 most HVG followed by dropping TCR, IG, and XY genes from the set. TRAV1-2, a MAIT marker, was then added back to the set. Batch key was set to sequencing lane with categorical covariate set to sample ID. Two additional continuous model covariates were used: percentage of counts from MT genes and total counts. A nearest-neighbor search (k = 30, kNN) was done on the model latent space. The nearest neighbor graph was used for Leiden clustering[14] and Uniform Manifold Approximation and Projection (UMAP)[15] coordinate calculation. Myeloid cell clusters were identified using top gene markers. This cell group contained various flavors of Monocytes and Dendritic cells. Platelet cells were similarly identified. 
<br>
<br>
To better annotated the myeloid cell group, it was re-integrated and re-clustered in the same way as was done before on all cells. Cells were labeled one of: cMono (classical monocyte, CD14+ monocyte), nMono (non-classical monocyte, CD16+ monocyte), cDC1 (conventional type 1 dendritic cell), cDC2 (conventional type 2 dendritic cell), or pDC (plasmacytoid dendritic cell). Additionally, a cluster showing evidence of both platelet and classical monocyte genes was given the label MPA (monocyte-platelet aggregate). Care was taken to determine whether the MPA group should be treated as persistent technical doublets or as a biologically relevant type.
<br>
<br>
To start, Solo-assigned doublets and Souporcell-assigned doublets were reintroduced into the complete, cleaned dataset. This dataset was integrated and clustered as before. The myeloid and platelet cell clusters were annotated. These two groups were then re-integrated and re-clustered. Any previously annotated doublets that expressed primarily platelet or myeloid markers were included in this integration and clustering. Several metrics were then calculated per cell cluster: % singlet, % Solo doublet, and % Souporcell doublet. Souporcell doublets define inter-sample, purely technical doublets. Within our dataset, Solo doublets define intra-sample doublets and are classified through doublet simulation. [note: results show rates of Souporcell doublets in MPA == rates of souporcell doublets in cMono; results show rates of Solo doublets in MPA >> rates of Solo doublets in cMono which is an expected outcome based on how Solo identifies doublets]. Rates for these three types were compared across cell clusters. Cells in clusters with much higher Souporcell doublet rates than others were annotated as "true" doublets. [note: MPA did not cluster with these groups].
<br>
<br>
Platelet- and cMono-labeled cells were then taken from the previous labeling that did not include Solo or Souporcell doublets. For each study participant, for each study day the Platelet and cMono counts matrices were subset and a cMono-Platelet doublet was simulated by randomly taking the counts from one Platelet and one cMono and summing them. This would ensure that simulated doublets (SimMPA) would only include counts summed from cells from the same study participant and same study day. Each Platelet in the dataset was only sampled once when summing counts. For each patient-study timepoint there were more cMono than Platelet. Subsequently, the final number of simulated doublets was equal to the number of Platelets in the dataset. This number was similar to the number of annotated MPAs. The final counts for the three cell types were: MPA (n = 3565), Platelet (n = 2843), SimMPA (n = 2843). SimMPAs were assigned barcode labels with a trailing '-dbl' tag. Similarly, metadata entries for these simulated doublets were mapped based on study participant and study day from which the counts were summed. Care was taken to count-normalize and log1p transform the simulated counts whenever comparisons were made to account for the inflated library sizes.
<br>
<br>
MAST methods go here

> REF  
[1] for citing plink: https://zzz.bwh.harvard.edu/plink/cite.shtml  
[2] see 'How should I cite Picard' here: https://broadinstitute.github.io/picard/faq.html  
[3] build obtained here: https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38; this may be useful: https://genome.ucsc.edu/goldenPath/credits.html#human_credits  
[4] Zheng, G. X. Y. et al. (2017). Massively parallel digital transcriptional profiling of single cells. Nature Communications 8: 1-12, doi:10.1038/ncomms14049  
[5] Alves, E. et al. Underrepresentation of activating KIR gene expression in single‐cell RNA‐seq data is due to KIR gene misassignment. Eur. J. Immunol. 54, e2350590 (2024)  
[6] Souporcell preprint: https://doi.org/10.1101/699637  
[7] Nicholas J. Bernstein, , Nicole L. Fong, Irene Lam, Margaret A. Roy, David G. Hendrickson, and David R. Kelley (2020), Solo: doublet identification in single-cell RNA-Seq via semi-supervised deep learning, Cell Systems.  
[8] see "Reference" section under the README here: https://github.com/scverse/scvi-tools  
[9] obtained through the R function citation(): R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.  
[10] for citing Seurat see: https://cran.r-project.org/web/packages/Seurat/citation.html  
[11] scverse tools suite: https://www.nature.com/articles/s41587-023-01733-8  
[12] for citing anndata see doi: 10.21105/joss.04371  
[13] note it's not stated how to cite python and many(most?) suggest you don't since it's so common  
[14] V. A. Traag, L. Waltman, and N. J. van Eck. From louvain to leiden: guaranteeing well-connected communities. Scientific Reports, mar 2019. URL: https://doi.org/10.1038/s41598-019-41695-z  
[15] umap technical source: https://arxiv.org/abs/1802.03426; umap less technical source: https://www.nature.com/articles/s43586-024-00363-x. Either work.  
[16] demuxafy: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03224-8
