library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

options(future.globals.maxSize = 10000000000)
options(Seurat.object.assay.version = "v5")


astro_crb <- readRDS("reference_astro_annotated.Rds")
astro_crb <- JoinLayers(astro_crb)

astro_crb@meta.data[['x']] <- NULL
astro_crb@meta.data[['y']] <- NULL
astro_crb@meta.data[['cluster']] <- NULL
astro_crb@meta.data[['nCount_SCT']] <- NULL
astro_crb@meta.data[['nFeature_SCT']] <- NULL
astro_crb@meta.data[['tissue']] <- NULL
astro_crb@meta.data[['dissociation']] <- NULL
astro_crb@meta.data[['chemistry']] <- NULL
astro_crb@meta.data[['lab']] <- NULL
astro_crb@meta.data[['cell_barcode.x']] <- NULL
astro_crb@meta.data[['library_label.x']] <- NULL
astro_crb@meta.data[['anatomical_division_label']] <- NULL
astro_crb@meta.data[['cell_label']] <- NULL
astro_crb@meta.data[['cell_barcode.y']] <- NULL
astro_crb@meta.data[['barcoded_cell_sample_label']] <- NULL
astro_crb@meta.data[['library_label.y']] <- NULL
astro_crb@meta.data[['feature_matrix_label']] <- NULL
astro_crb@meta.data[['entity']] <- NULL
astro_crb@meta.data[['brain_section_label']] <- NULL
astro_crb@meta.data[['library_method']] <- NULL
astro_crb@meta.data[['region_of_interest_acronym']] <- NULL
astro_crb@meta.data[['donor_label']] <- NULL
astro_crb@meta.data[['donor_genotype']] <- NULL
astro_crb@meta.data[['donor_sex']] <- NULL
astro_crb@meta.data[['dataset_label']] <- NULL
astro_crb@meta.data[['cluster_alias']] <- NULL
astro_crb@meta.data[['neurotransmitter']] <- NULL
astro_crb@meta.data[['class']] <- NULL
astro_crb@meta.data[['subclass']] <- NULL
astro_crb@meta.data[['supertype']] <- NULL
astro_crb@meta.data[['neurotransmitter_color']] <- NULL
astro_crb@meta.data[['class_color']] <- NULL
astro_crb@meta.data[['subclass_color']] <- NULL
astro_crb@meta.data[['supertype_color']] <- NULL
astro_crb@meta.data[['cluster_color']] <- NULL
astro_crb@meta.data[['region_of_interest_order']] <- NULL
astro_crb@meta.data[['region_of_interest_color']] <- NULL
astro_crb@meta.data[['genotype']] <- NULL
astro_crb@meta.data[['classes']] <- NULL
astro_crb@meta.data[['split']] <- NULL
astro_crb@meta.data[['Vascular1']] <- NULL
astro_crb@meta.data[['endothelial1']] <- NULL
astro_crb@meta.data[['Astro1']] <- NULL
astro_crb@meta.data[['Purkinje1']] <- NULL
astro_crb@meta.data[['Neurons_Dcn_glu1']] <- NULL
astro_crb@meta.data[['Oligo1']] <- NULL
astro_crb@meta.data[['Granule_cells_mature1']] <- NULL
astro_crb@meta.data[['Golgi1']] <- NULL
astro_crb@meta.data[['PLI1']] <- NULL
astro_crb@meta.data[['MLI1']] <- NULL
astro_crb@meta.data[['Gaba_intermediate_prog1']] <- NULL
astro_crb@meta.data[['pan_gaba1']] <- NULL
astro_crb@meta.data[['pan_glu1']] <- NULL
astro_crb@meta.data[['pan_neuron1']] <- NULL
astro_crb@meta.data[['broad_class']] <- NULL
astro_crb@meta.data[['cell_type']] <- NULL
astro_crb@meta.data[['integrated_clusters_astro_K0.3']] <- NULL
astro_crb@meta.data[['seurat_clusters']] <- NULL
astro_crb@meta.data[['MARKERS1']] <- NULL

astro_crb$New_annot <- NA
astro_crb@meta.data[astro_crb@meta.data$integrated_clusters_astro %in% c("1","3","4","6"),]$New_annot <- "BG"
astro_crb@meta.data[astro_crb@meta.data$integrated_clusters_astro %in% c("0","2","10"),]$New_annot <- "CNA"
astro_crb@meta.data[astro_crb@meta.data$integrated_clusters_astro %in% c("7","9"),]$New_annot <- "GLA"
astro_crb@meta.data[astro_crb@meta.data$integrated_clusters_astro == "5",]$New_annot <- "WMA"
astro_crb@meta.data[astro_crb@meta.data$integrated_clusters_astro == "8",]$New_annot <- "Astro_glu"




#astro_crb <- FindClusters(astro_crb, resolution = 0.3, cluster.name = "cca_clusters")

#adult_crb <- RunPCA(adult_crb) #only to run if the PCA hasn't been run before

sce <- as.SingleCellExperiment(astro_crb, assay = "RNA")

reducedDims(sce) <- list(PCA=astro_crb@reductions$integrated.cca.astro@cell.embeddings)
#reducedDims(sce) <- list(PCA=adult_crb@reductions$pca@cell.embeddings)



writeH5AD(sce, file="astrocyte_subtypes_from11-30annot.h5ad")

