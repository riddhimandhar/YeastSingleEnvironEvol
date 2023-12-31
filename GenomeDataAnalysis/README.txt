Prorgams used for processing of NGS data (Illumina platform 2 x 150 bp) 
===================================================================

These programs are run in the following order. 


1. mapping_chromosomes.pl 
    - map chromosome names to NCBI ids 

2. process_samfiles.pl 
   - process generated sam files 
   - output - Alignments and Mutation list 

3. PE_merge.pl 
   - merge mapping information and mutation information from paired-end (PE) reads  

4. calc_coverage.pl 
   - calculate sequencing coverage across the genome  

5. find_mutations.pl
   - call mutations after all filtering 

6. mut_anc_evol.pl
   - compare mutations between ancestral and evolved lines 

7. mutation_annotation.pl
   - classify mutations into coding, noncoding, synonymous etc. 

8. statistics_mutations.pl

   - identify common and unique mutations among the evolved lines 
   - identify common and unique genes that are hit by mutations in the evolved lines 

9. check_gene_properties.pl 

   - check genes that contain mutations in our experiment for expression changes under 
     different stressors (data from Gasch et al., MBoC 2000, Causton et al., MBoC 2001) 

10. knockout_phenotype_mapping.pl 

   - check fitness of knockouts of genes that contain mutations in our experiment 
     (data from Warringer et al., PNAS 2003; Tucker and Fields, Comp. Func. Genomics 2004; 
      Brown et al., Mol. Syst. Biol. 2006; St Onge et al., Nat. Genet. 2007;
      Hillenmeyer et al., Science 2008; Baryshnikova et al., Nat. Met. 2010)


11. CNV_analysis_P1.pl 
    CNV_analysis_P2_cov_v2.pl 

    - Two step process for normalization and identification of copy number variations (CNVs)
      in the evolved lines compared to the ancestral lines 




