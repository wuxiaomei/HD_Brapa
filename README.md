# HD_Brapa
Compound Heat-Drought Stress in Brassica rapa

# Codes for Brassica rapa HD analysis (codes written by Xiaomei Wu, Donghui Hu, and Rachelle R.Q. Lee)

## 1. WGCNA analysis of individual stress-responsive gene groups
1. Identify WGCNA modules (Supplemental Figure S8)
1. Plot gene expression profiles (Figures 4A-C, 7C and Supplemental Figure S9)
1. GO enrichment analysis of each module (Figure 4D, Supplemental Figure S12)

## 2. Reconstruct the Gene regulatory network (GRN), keyGRN and TF-TF network
1. GRN and TF-TF GRN (Supplemental Figure S13)
1. keyGRN

## 3. Network analysis (https://github.com/rlrq/HDstress)
1. Network dispersion analysis for GRN, keyGRN, and TF-TF network (Figure 5 and Supplemental Figure S14)
1. Master regulatory nodes analysis (Figure 6 and Supplemental Figure S15)

## 4. Modified GO enrichment analysis of BLAST-mappable and -unmappable genes to A. thaliana orthologues (https://github.com/rlrq/HDstress)
1. Supplemental Figures S5: Gene counts of GO terms mapped via orthologous families with different normalisation methods against counts of GO terms mapped via A. thaliana orthologue identified by BLAST.
1. Supplemental Figure S6: GO enrichment analysis of all B. rapa genes with different mapping methods and universes.
1. Supplemental Figure S7: GO enrichment analysis of treatment DEGs (HD-responsible genes) with different mapping methods and universes.

## 5. Identify ancestral GO terms (https://github.com/rlrq/HDstress)
1. Figure 3B and 4D, Supplemental Figure S3 and S11: GO hierarchy containing GO terms of interest by mapping higher level GO terms to lower level terms

## 6. Detect putative cis-regulatory elements (pCRE)
1. Supplemental Figures S16 and S17: pCREs enriched in the target genes of a TF
1. Data
    1. install MEMEsuit and download motif database (https://meme-suite.org/meme/meme-software/)
1. Scripts
    1. predict_pcre.sh: steps of the pipeline
    1. predict_pcre.pl: main perl script to detect k-mer motifs and annotate the motifs using MEMEsuit
    1. predict_pcre.fisher.R: R script to do fisher.test and adjust p-values

## 7. Other plots
1. Figure 1D and E: PCA analysis of 48 samples
1. Figure 3A: Significant GO terms for the four HD-responsive gene groups under treatment comparisons and accession comparisons
1. Supplemental Figure 4: Mapped and unmapped genes by BLAST to A. thaliana genes
