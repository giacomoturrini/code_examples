#load broad data cell metadata (you need to remove the first line otherwise the coordinates are saved as a character vector) and count matrix
crnl10_sp <- read.csv("well10_spatial.csv")
crnl10_sp <- crnl10_sp[-1,]
crnl10_sp <- write.csv(crnl10_sp, "well10_corrected_spatial.csv", row.names = FALSE)
crnl10_sp <- read.csv("well10_corrected_spatial.csv")
crnl10raw <- read.csv("well10raw_expression_pd.csv")

# subset the cells based on molecular tissue region or on coordinates
#cerebellar_cells <- sgtl3_sp[which(sgtl3_sp$Main_molecular_tissue_region %in% c("CB_1","CB_2","FbTrt")),]
xmax = 32000
xmin = 7000
ymax = 23000
ymin = 300
cerebellar_cells <- crnl10_sp[which(crnl10_sp$X <= xmax & crnl10_sp$X >= xmin & crnl10_sp$Y <= ymax & crnl10_sp$Y >= ymin) ,]

cerebellar_counts <- crnl10raw[,c("GENE",cerebellar_cells$NAME)]

write.csv(cerebellar_cells, "well10_cerebellum_cells_metadata.csv", row.names = FALSE)
write.csv(cerebellar_counts, "well10_cerebellum_counts.csv", row.names = FALSE)