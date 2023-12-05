library(sf)
library(RImageJROI)
library(pracma)
library(dplyr)
library(parallel)
library(pbapply)

Sys.setenv("R_CSTACK_LIMIT" = "10000000")

## from github https://github.com/davidcsterratt/RImageJROI
#source("~/Desktop/ImageJ/RImageJROI-master/R/read.ijroi.R")


roi<-list.files("/home/giacomo/Downloads/P60_SAMPLE1/RoiSet/",full.names = TRUE)

cl <- makeCluster(10L)
rois_list<-pblapply(roi,function(roi){
  library(RImageJROI)
  
  IJ<-read.ijroi(roi)
  read.ijroi(roi)$coords
  
  },cl=cl)



spots_data <- read.csv("/home/giacomo/Downloads/P60_SAMPLE1/decoding_PRMC_ALBP60_sample1_106genes_PRMC_filt_min_04.csv")
spots_data <- spots_data[, c("xc", "yc","target")]




duplicate_rows <- spots_data[duplicated(spots_data), ]

# Print the duplicate rows (if any)
if (nrow(duplicate_rows) > 0) {
  print("Duplicate rows found:")
  print(duplicate_rows)
} else {
  print("No duplicate rows found.")
}

#cb_roi <- read.ijroi("/home/giacomo/Downloads/P60_SAMPLE1/crops/real_crops/crop_CNA_sample1.roi")
#cb_roi <- cb_roi$coords
#cb_roi <- rbind(cb_roi, cb_roi[1,])

#spots_data <- spots_data[(inpolygon(spots_data[,1],spots_data[,2],cb_roi[,1],cb_roi[,2], boundary=TRUE)),]


spot_count<-pblapply(rois_list,count_spots_inside_roi_by_gene <- function(roi, spots) {
  library(sf)
  library(pracma)
  library(dplyr)
  
  # Open the manually obtained cb ROI
  #cb_roi <- read.ijroi("/home/giacomo/Downloads/P60_SAMPLE1/crops/real_crops/crop_CNA_sample1.roi")
  #cb_roi <- cb_roi$coords
  #cb_roi <- rbind(cb_roi, cb_roi[1,])
  
  
  # Define pixel siye and length (micrometers) to which expand the ROI
    pixel_size = 0.1625
    length_um= 2
    length= length_um/pixel_size
    # Ensure the ROI is closed by adding the starting point to the end
    roi <- rbind(roi, roi[1,])
    polygon_nuc <- st_polygon(list(as.matrix(roi)))
    centroid_nuc <- st_centroid(polygon_nuc)
    #if(inpolygon(centroid_nuc[1],centroid_nuc[2],cb_roi[,1],cb_roi[,2])){
    polygon <- st_buffer(polygon_nuc, length)
    roi_buff <- polygon[[1]]
    centroid <- st_centroid(polygon)
    
    centroid<-paste(round(as.numeric(centroid))[1],round(as.numeric(centroid))[2],sep="_")
    points <- st_multipoint(as.matrix(spots[,c("xc", "yc")]))
    
    
    
    spots$in_polygon<-inpolygon(spots[,1], spots[,2], roi_buff[,1], roi_buff[,2], boundary=TRUE)
    
    
    spots$in_polygon<-as.numeric(spots$in_polygon)
    count=spots%>%group_by(target)%>%summarise(count=sum(in_polygon))
    
    count<-data.frame(count)
    rownames(count)<-count$target
    count$target<-NULL
    
    colnames(count)<-centroid
    
    return(count)
    #}
  }
  ,spots=spots_data,cl=cl)

spot_count_inCB_only <- spot_count[!sapply(spot_count, is.null)]


final_count<-do.call("cbind",spot_count_inCB_only)

final_count
count_spots <- sapply(final_count,sum)
barplot(table(count_spots)/length(count_spots)*100)

stopCluster(cl = cl)
write.csv(final_count, file= "",row.names = TRUE)

#create a dataframe with the coordinates of the centroid of each nucleus
name_of_cells <- colnames(final_count)
coordinates <- strsplit(name_of_cells, split = "_")
coord <- data.frame(row.names=name_of_cells)
coord$x <- rep(NA, length(coordinates))
coord$y <- rep(NA, length(coordinates))
for(i in seq(1, length(coordinates))){
  coord[i,]$x <- coordinates[[i]][1]
  coord[i,]$y <- coordinates[[i]][2]
}

write.csv(coord, file = "Data/Spatial/P0_7/coordinates.csv", row.names = TRUE)

