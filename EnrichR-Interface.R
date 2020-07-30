
require(enrichR)

EnrichR.dbs = c("ARCHS4_Tissues",
                        "GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018",
                        "KEGG_2016","ChEA_2016",
                        "Human_Gene_Atlas",
                        "GeneSigDB", 
                        "Cancer_Cell_Line_Encyclopedia", "Disease_Perturbations_from_GEO_down",
                        "Disease_Perturbations_from_GEO_up", "Human_Phenotype_Ontology"
)

# Jenia's version
EnrichR.dbs=listEnrichrDbs()
EnrichR.dbs=EnrichR.dbs$libraryName
EnrichR.removedbs=c("Rare_Diseases_GeneRIF_ARCHS4_Predictions","Rare_Diseases_AutoRIF_ARCHS4_Predictions","NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions",
            "NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions","Ligand_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_up",
            "Ligand_Perturbations_from_GEO_up","MCF7_Perturbations_from_GEO_down","MCF7_Perturbations_from_GEO_up","Microbe_Perturbations_from_GEO_down",
            "Microbe_Perturbations_from_GEO_up","LINCS_L1000_Ligand_Perturbations_down" ,"LINCS_L1000_Ligand_Perturbations_up","LINCS_L1000_Kinase_Perturbations_down",
            "LINCS_L1000_Kinase_Perturbations_up","NIH_Funded_PIs_2017_Human_GeneRIF","NIH_Funded_PIs_2017_Human_AutoRIF","Rare_Diseases_GeneRIF_Gene_Lists",
            "Rare_Diseases_AutoRIF_Gene_Lists","DisGeNET" ,"Drug_Perturbations_from_GEO_2014","SILAC_Phosphoproteomics","Old_CMAP_up","Old_CMAP_down","HMDB_Metabolites",
            "Allen_Brain_Atlas_up","Allen_Brain_Atlas_down","Achilles_fitness_decrease","Achilles_fitness_increase","DrugMatrix","Jensen_DISEASES",
            "HMS_LINCS_KinomeScan", "Genes_Associated_with_NIH_Grants", 
            "Data_Acquisition_Method_Most_Popular_Genes", "Enrichr_Libraries_Most_Popular_Genes", "ENCODE_Histone_Modifications_2015",
            "Epigenomics_Roadmap_HM_ChIP-seq", "ESCAPE", "ChEA_2013", "ChEA_2015", "BioCarta_2013", "BioCarta_2015", "Disease_Signatures_from_GEO_down_2014",
            "Disease_Signatures_from_GEO_up_2014", "ENCODE_Histone_Modifications_2013", "ENCODE_TF_ChIP-seq_2014", "Enrichr_Submissions_TF-Gene_Coocurrence",
            "GO_Biological_Process_2013", "GO_Biological_Process_2015", "GO_Biological_Process_2017", "GO_Biological_Process_2017b", "GO_Cellular_Component_2013",
            "GO_Cellular_Component_2015", "GO_Cellular_Component_2017", "GO_Cellular_Component_2017b", "GO_Molecular_Function_2013", "GO_Molecular_Function_2015",
            "GO_Molecular_Function_2017", "GO_Molecular_Function_2017b", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down", "KEGG_2013", "KEGG_2015", "KEGG_2016",
            "MGI_Mammalian_Phenotype_2017", "MGI_Mammalian_Phenotype_Level_4_2019", "MSigDB_Computational", "NCI-Nature_2015", "Pfam_Domains_2019", "InterPro_Domains_2019",
            "Reactome_2013", "Reactome_2015", "WikiPathways_2013", "WikiPathways_2015", "WikiPathways_2016",  "TargetScan_microRNA_2017",
            "Panther_2015", "MGI_Mammalian_Phenotype_2013", "HumanCyc_2015", "NCI-60_Cancer_Cell_Lines", "KEA_2013", "SysMyo_Muscle_Gene_Sets" 
            )
EnrichR.dbs=setdiff(EnrichR.dbs,EnrichR.removedbs)

#
# Assumes a sorted list (according to significance or fold change)
#
CollectEnrichRAnnotation = function( Diff.list.sorted, 
                                     sorted=TRUE, 
                                     DBs = EnrichR.dbs,
                                     Q.value.threshold = 0.1) {
  list.temp = list()
  n.list =  c(100,500,1000,2000)
  if( !sorted )
    n.list = c(length(Diff.list.sorted))
  for( n in n.list ) 
    if( n == 100 || length(Diff.list.sorted) >= n ) {
      n = min(n, length(Diff.list.sorted))
      cat("...", n, "\n")
      enriched <- enrichr(Diff.list.sorted[1:n], DBs)
      en.list = lapply(names(enriched), function(db) {
        if(nrow(enriched[[db]]) > 0) {
          x = enriched[[db]]
          x$DB = db
          x$QuerySize = n
          x = x[,c("DB", "Term", "QuerySize", "Overlap", "P.value", "Adjusted.P.value","Combined.Score", "Genes")]
        }
      })
      list.temp[[as.character(n)]] = en.list
    }
  cat("\n")
  
  if( length(list.temp) == 1 ) {
    list.temp = list.temp[[1]]
  } else
    list.temp = do.call(c,list.temp)
  enriched = do.call(rbind,list.temp)
  
  if( nrow(enriched) > 0 ) {
    enriched = enriched[order(enriched$Combined.Score, decreasing = TRUE),]
    enriched = enriched[enriched$Adjusted.P.value < Q.value.threshold, ]
  }
  if( nrow(enriched) > 0 ) {
    I = duplicated(enriched[,c("Term", "DB")])
    if(any(I))
      enriched = enriched[!I,]
    enriched$Term = gsub(",", ";", enriched$Term)
    
  }
  return(enriched)
}

