
#--read and process annotation data
#2_3
setwd("./pathlogical_annotation/")
spot_2_3 <- read.csv(pathlogical_annotation[1])
spot_2_3[, 1] <- sub("1", "2_3", spot_2_3[, 1])
spot_2_3_specifc <- read.csv(pathlogical_annotation[2])
spot_2_3_specifc[, 1] <- sub("1", "2_3", spot_2_3_specifc[, 1])
table(sample6.combined$patientID)
length(intersect(spot_2_3$Barcode, colnames(sample6.combined)))#4348 correct
#2_8
spot_2_8 <- read.csv(pathlogical_annotation[3])
spot_2_8[, 1] <- sub("1", "2_8", spot_2_8[, 1])
spot_2_8_specific <- read.csv(pathlogical_annotation[4])
spot_2_8_specific[, 1] <- sub("1", "2_8", spot_2_8_specific[, 1])
length(intersect(spot_2_8_specific$Barcode, colnames(sample6.combined)))#3445 correct
#T4857
spot_T4857 <- read.csv(pathlogical_annotation[5])
spot_T4857[, 1] <- sub("1", "T4857", spot_T4857[, 1])
spot_T4857_specific <- read.csv(pathlogical_annotation[6])
spot_T4857_specific[, 1] <-
  sub("1", "T4857", spot_T4857_specific[, 1])
length(intersect(spot_T4857$Barcode, colnames(sample6.combined)))#4832 correct
spot_annotation <- rbind(spot_2_3, spot_2_8, spot_T4857)
#change blank into level4
spot_annotation[, 2][which(spot_annotation[, 2] == "")] <-
  "non_AT8_level4"
spot_specifc <-
  rbind(spot_2_3_specifc, spot_2_8_specific, spot_T4857_specific)
spot_specifc <- spot_specifc[-which(spot_specifc$AT8 == ""), ]
