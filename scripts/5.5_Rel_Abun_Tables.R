### Script to obtain info for manuscript - 23/07/2025

library(dplyr)
library(tidyr)
library(purrr)
library(openxlsx)   # for Excel export

# Main datasets
BACT<-melted[which(melted$Barcode=="BACT"),];dim(BACT)
FUNG<-melted[which(melted$Barcode=="FUNG"),];dim(FUNG)
POLL<-melted[which(melted$Barcode=="PLAN"),];dim(POLL)

# Bacterial datasets
BACT_MTL<-BACT[which(BACT$city=="Montreal"),];dim(BACT_MTL)
BACT_MTL_spring<-BACT_MTL[which(BACT_MTL$time=="Spring"),];dim(BACT_MTL_spring)
BACT_MTL_summer<-BACT_MTL[which(BACT_MTL$time=="Summer"),];dim(BACT_MTL_summer)
BACT_MTL_fall<-BACT_MTL[which(BACT_MTL$time=="Fall"),];dim(BACT_MTL_fall)

BACT_QBC<-BACT[which(BACT$city=="Quebec"),];dim(BACT_QBC)
BACT_QBC_spring<-BACT_QBC[which(BACT_QBC$time=="Spring"),];dim(BACT_QBC_spring)
BACT_QBC_summer<-BACT_QBC[which(BACT_QBC$time=="Summer"),];dim(BACT_QBC_summer)
BACT_QBC_fall<-BACT_QBC[which(BACT_QBC$time=="Fall"),];dim(BACT_QBC_fall)

BACT_SHB<-BACT[which(BACT$city=="Sherbrooke"),];dim(BACT_SHB)
BACT_SHB_spring<-BACT_SHB[which(BACT_SHB$time=="Spring"),];dim(BACT_SHB_spring)
BACT_SHB_summer<-BACT_SHB[which(BACT_SHB$time=="Summer"),];dim(BACT_SHB_summer)
BACT_SHB_fall<-BACT_SHB[which(BACT_SHB$time=="Fall"),];dim(BACT_SHB_fall)

# Fungal datasets
FUNG_MTL<-FUNG[which(FUNG$city=="Montreal"),];dim(FUNG_MTL)
FUNG_MTL_spring<-FUNG_MTL[which(FUNG_MTL$time=="Spring"),];dim(FUNG_MTL_spring)
FUNG_MTL_summer<-FUNG_MTL[which(FUNG_MTL$time=="Summer"),];dim(FUNG_MTL_summer)
FUNG_MTL_fall<-FUNG_MTL[which(FUNG_MTL$time=="Fall"),];dim(FUNG_MTL_fall)

FUNG_QBC<-FUNG[which(FUNG$city=="Quebec"),];dim(FUNG_QBC)
FUNG_QBC_spring<-FUNG_QBC[which(FUNG_QBC$time=="Spring"),];dim(FUNG_QBC_spring)
FUNG_QBC_summer<-FUNG_QBC[which(FUNG_QBC$time=="Summer"),];dim(FUNG_QBC_summer)
FUNG_QBC_fall<-FUNG_QBC[which(FUNG_QBC$time=="Fall"),];dim(FUNG_QBC_fall)

FUNG_SHB<-FUNG[which(FUNG$city=="Sherbrooke"),];dim(FUNG_SHB)
FUNG_SHB_spring<-FUNG_SHB[which(FUNG_SHB$time=="Spring"),];dim(FUNG_SHB_spring)
FUNG_SHB_summer<-FUNG_SHB[which(FUNG_SHB$time=="Summer"),];dim(FUNG_SHB_summer)
FUNG_SHB_fall<-FUNG_SHB[which(FUNG_SHB$time=="Fall"),];dim(FUNG_SHB_fall)

# Plant datasets
POLL_MTL<-POLL[which(POLL$city=="Montreal"),];dim(POLL_MTL)
POLL_MTL_spring<-POLL_MTL[which(POLL_MTL$time=="Spring"),];dim(POLL_MTL_spring)
POLL_MTL_summer<-POLL_MTL[which(POLL_MTL$time=="Summer"),];dim(POLL_MTL_summer)
POLL_MTL_fall<-POLL_MTL[which(POLL_MTL$time=="Fall"),];dim(POLL_MTL_fall)

POLL_QBC<-POLL[which(POLL$city=="Quebec"),];dim(POLL_QBC)
POLL_QBC_spring<-POLL_QBC[which(POLL_QBC$time=="Spring"),];dim(POLL_QBC_spring)
POLL_QBC_summer<-POLL_QBC[which(POLL_QBC$time=="Summer"),];dim(POLL_QBC_summer)
POLL_QBC_fall<-POLL_QBC[which(POLL_QBC$time=="Fall"),];dim(POLL_QBC_fall)

POLL_SHB<-POLL[which(POLL$city=="Sherbrooke"),];dim(POLL_SHB)
POLL_SHB_spring<-POLL_SHB[which(POLL_SHB$time=="Spring"),];dim(POLL_SHB_spring)
POLL_SHB_summer<-POLL_SHB[which(POLL_SHB$time=="Summer"),];dim(POLL_SHB_summer)
POLL_SHB_fall<-POLL_SHB[which(POLL_SHB$time=="Fall"),];dim(POLL_SHB_fall)

#Bacterial datasets for relative abundance table
BACTt<-t<-aggregate(BACT$relAb,list(Family=BACT$Family),sum)
BACTt$x<-BACTt$x/(sum(BACT$relAb))*100
BACT_MTLt<-t<-aggregate(BACT_MTL$relAb,list(Family=BACT_MTL$Family),sum)
BACT_MTLt$x<-BACT_MTLt$x/(sum(BACT_MTL$relAb))*100
BACT_MTL_springt<-t<-aggregate(BACT_MTL_spring$relAb,list(Family=BACT_MTL_spring$Family),sum)
BACT_MTL_springt$x<-BACT_MTL_springt$x/(sum(BACT_MTL_spring$relAb))*100
BACT_MTL_summert<-t<-aggregate(BACT_MTL_summer$relAb,list(Family=BACT_MTL_summer$Family),sum)
BACT_MTL_summert$x<-BACT_MTL_summert$x/(sum(BACT_MTL_summer$relAb))*100
BACT_MTL_fallt<-t<-aggregate(BACT_MTL_fall$relAb,list(Family=BACT_MTL_fall$Family),sum)
BACT_MTL_fallt$x<-BACT_MTL_fallt$x/(sum(BACT_MTL_fall$relAb))*100

BACT_QBCt<-t<-aggregate(BACT_QBC$relAb,list(Family=BACT_QBC$Family),sum)
BACT_QBCt$x<-BACT_QBCt$x/(sum(BACT_QBC$relAb))*100
BACT_QBC_springt<-t<-aggregate(BACT_QBC_spring$relAb,list(Family=BACT_QBC_spring$Family),sum)
BACT_QBC_springt$x<-BACT_QBC_springt$x/(sum(BACT_QBC_spring$relAb))*100
BACT_QBC_fallt<-t<-aggregate(BACT_QBC_fall$relAb,list(Family=BACT_QBC_fall$Family),sum)
BACT_QBC_fallt$x<-BACT_QBC_fallt$x/(sum(BACT_QBC_fall$relAb))*100

BACT_SHBt<-t<-aggregate(BACT_SHB$relAb,list(Family=BACT_SHB$Family),sum)
BACT_SHBt$x<-BACT_SHBt$x/(sum(BACT_SHB$relAb))*100
BACT_SHB_springt<-t<-aggregate(BACT_SHB_spring$relAb,list(Family=BACT_SHB_spring$Family),sum)
BACT_SHB_springt$x<-BACT_SHB_springt$x/(sum(BACT_SHB_spring$relAb))*100
BACT_SHB_summert<-t<-aggregate(BACT_SHB_summer$relAb,list(Family=BACT_SHB_summer$Family),sum)
BACT_SHB_summert$x<-BACT_SHB_summert$x/(sum(BACT_SHB_summer$relAb))*100
BACT_SHB_fallt<-t<-aggregate(BACT_SHB_fall$relAb,list(Family=BACT_SHB_fall$Family),sum)
BACT_SHB_fallt$x<-BACT_SHB_fallt$x/(sum(BACT_SHB_fall$relAb))*100

#Fungal datasets
FUNGt<-t<-aggregate(FUNG$relAb,list(Family=FUNG$Family),sum)
FUNGt$x<-FUNGt$x/(sum(FUNG$relAb))*100
FUNG_MTLt<-t<-aggregate(FUNG_MTL$relAb,list(Family=FUNG_MTL$Family),sum)
FUNG_MTLt$x<-FUNG_MTLt$x/(sum(FUNG_MTL$relAb))*100
FUNG_MTL_springt<-t<-aggregate(FUNG_MTL_spring$relAb,list(Family=FUNG_MTL_spring$Family),sum)
FUNG_MTL_springt$x<-FUNG_MTL_springt$x/(sum(FUNG_MTL_spring$relAb))*100
FUNG_MTL_summert<-t<-aggregate(FUNG_MTL_summer$relAb,list(Family=FUNG_MTL_summer$Family),sum)
FUNG_MTL_summert$x<-FUNG_MTL_summert$x/(sum(FUNG_MTL_summer$relAb))*100
FUNG_MTL_fallt<-t<-aggregate(FUNG_MTL_fall$relAb,list(Family=FUNG_MTL_fall$Family),sum)
FUNG_MTL_fallt$x<-FUNG_MTL_fallt$x/(sum(FUNG_MTL_fall$relAb))*100

FUNG_QBCt<-t<-aggregate(FUNG_QBC$relAb,list(Family=FUNG_QBC$Family),sum)
FUNG_QBCt$x<-FUNG_QBCt$x/(sum(FUNG_QBC$relAb))*100
FUNG_QBC_springt<-t<-aggregate(FUNG_QBC_spring$relAb,list(Family=FUNG_QBC_spring$Family),sum)
FUNG_QBC_springt$x<-FUNG_QBC_springt$x/(sum(FUNG_QBC_spring$relAb))*100
FUNG_QBC_fallt<-t<-aggregate(FUNG_QBC_fall$relAb,list(Family=FUNG_QBC_fall$Family),sum)
FUNG_QBC_fallt$x<-FUNG_QBC_fallt$x/(sum(FUNG_QBC_fall$relAb))*100

FUNG_SHBt<-t<-aggregate(FUNG_SHB$relAb,list(Family=FUNG_SHB$Family),sum)
FUNG_SHBt$x<-FUNG_SHBt$x/(sum(FUNG_SHB$relAb))*100
FUNG_SHB_springt<-t<-aggregate(FUNG_SHB_spring$relAb,list(Family=FUNG_SHB_spring$Family),sum)
FUNG_SHB_springt$x<-FUNG_SHB_springt$x/(sum(FUNG_SHB_spring$relAb))*100
FUNG_SHB_summert<-t<-aggregate(FUNG_SHB_summer$relAb,list(Family=FUNG_SHB_summer$Family),sum)
FUNG_SHB_summert$x<-FUNG_SHB_summert$x/(sum(FUNG_SHB_summer$relAb))*100
FUNG_SHB_fallt<-t<-aggregate(FUNG_SHB_fall$relAb,list(Family=FUNG_SHB_fall$Family),sum)
FUNG_SHB_fallt$x<-FUNG_SHB_fallt$x/(sum(FUNG_SHB_fall$relAb))*100

#Plant datasets
POLLt<-t<-aggregate(POLL$relAb,list(Family=POLL$Family),sum)
POLLt$x<-POLLt$x/(sum(POLL$relAb))*100
POLL_MTLt<-t<-aggregate(POLL_MTL$relAb,list(Family=POLL_MTL$Family),sum)
POLL_MTLt$x<-POLL_MTLt$x/(sum(POLL_MTL$relAb))*100
POLL_MTL_springt<-t<-aggregate(POLL_MTL_spring$relAb,list(Family=POLL_MTL_spring$Family),sum)
POLL_MTL_springt$x<-POLL_MTL_springt$x/(sum(POLL_MTL_spring$relAb))*100
POLL_MTL_summert<-t<-aggregate(POLL_MTL_summer$relAb,list(Family=POLL_MTL_summer$Family),sum)
POLL_MTL_summert$x<-POLL_MTL_summert$x/(sum(POLL_MTL_summer$relAb))*100
POLL_MTL_fallt<-t<-aggregate(POLL_MTL_fall$relAb,list(Family=POLL_MTL_fall$Family),sum)
POLL_MTL_fallt$x<-POLL_MTL_fallt$x/(sum(POLL_MTL_fall$relAb))*100

POLL_QBCt<-t<-aggregate(POLL_QBC$relAb,list(Family=POLL_QBC$Family),sum)
POLL_QBCt$x<-POLL_QBCt$x/(sum(POLL_QBC$relAb))*100
POLL_QBC_springt<-t<-aggregate(POLL_QBC_spring$relAb,list(Family=POLL_QBC_spring$Family),sum)
POLL_QBC_springt$x<-POLL_QBC_springt$x/(sum(POLL_QBC_spring$relAb))*100
POLL_QBC_fallt<-t<-aggregate(POLL_QBC_fall$relAb,list(Family=POLL_QBC_fall$Family),sum)
POLL_QBC_fallt$x<-POLL_QBC_fallt$x/(sum(POLL_QBC_fall$relAb))*100

POLL_SHBt<-t<-aggregate(POLL_SHB$relAb,list(Family=POLL_SHB$Family),sum)
POLL_SHBt$x<-POLL_SHBt$x/(sum(POLL_SHB$relAb))*100
POLL_SHB_springt<-t<-aggregate(POLL_SHB_spring$relAb,list(Family=POLL_SHB_spring$Family),sum)
POLL_SHB_springt$x<-POLL_SHB_springt$x/(sum(POLL_SHB_spring$relAb))*100
POLL_SHB_summert<-t<-aggregate(POLL_SHB_summer$relAb,list(Family=POLL_SHB_summer$Family),sum)
POLL_SHB_summert$x<-POLL_SHB_summert$x/(sum(POLL_SHB_summer$relAb))*100
POLL_SHB_fallt<-t<-aggregate(POLL_SHB_fall$relAb,list(Family=POLL_SHB_fall$Family),sum)
POLL_SHB_fallt$x<-POLL_SHB_fallt$x/(sum(POLL_SHB_fall$relAb))*100

#t<-aggregate(data_id$relAb,list(Family=data_id$Family),sum)
#t$x<-t$x/(sum(data_id$relAb))*100
#sum(t$x)
#t[order(t$x,decreasing = T)[1:30],]

#Phylum
BACTp<-melted[which(melted$Barcode=="BACT"),];dim(BACTp)
FUNGp<-melted[which(melted$Barcode=="FUNG"),];dim(FUNGp)
POLLp<-melted[which(melted$Barcode=="PLAN"),];dim(POLLp)
sum(BACTp$relAb)#39

# BACTERIA TABLE FOR MANUSCRIP
# Step 1: Get top 15 families from global dataset
top_families <- BACTt %>%
  group_by(Family) %>%
  summarise(Global = sum(x, na.rm = TRUE)) %>%
  arrange(desc(Global)) %>%
  slice_head(n = 15)

# Step 2: Define datasets to include
datasets <- list(
  BAC = BACTt,
  MTL = BACT_MTLt,
  MTL_spring = BACT_MTL_springt,
  MTL_summer = BACT_MTL_summert,
  MTL_fall = BACT_MTL_fallt,
  QBC = BACT_QBCt,
  QBC_spring = BACT_QBC_springt,
  QBC_fall = BACT_QBC_fallt,
  SHB = BACT_SHBt,
  SHB_spring = BACT_SHB_springt,
  SHB_summer = BACT_SHB_summert,
  SHB_fall = BACT_SHB_fallt
)

# Step 3: Extract relative abundance for top families from each dataset
abundance_table <- map_dfr(names(datasets), function(name) {
  datasets[[name]] %>%
    group_by(Family) %>%
    summarise(!!name := sum(x, na.rm = TRUE)) %>%
    filter(Family %in% top_families$Family)
}, .id = "source")
pas<-15
abundance_table$MTL[1:15]<-abundance_table$MTL[(pas+1):(2*pas)]
abundance_table$MTL_spring[1:15]<-abundance_table$MTL_spring[(2*pas+1):(3*pas)]
abundance_table$MTL_summer[1:15]<-abundance_table$MTL_summer[(3*pas+1):(4*pas)]
abundance_table$MTL_fall[1:15]<-abundance_table$MTL_fall[(4*pas+1):(5*pas)]
abundance_table$QBC[1:15]<-abundance_table$QBC[(5*pas+1):(6*pas)]
abundance_table$QBC_spring[1:15]<-abundance_table$QBC_spring[(6*pas+1):(7*pas)]
abundance_table$QBC_fall[1:15]<-abundance_table$QBC_fall[(7*pas+1):(8*pas)]
abundance_table$SHB[1:15]<-abundance_table$SHB[(8*pas+1):(9*pas)]
abundance_table$SHB_spring[1:15]<-abundance_table$SHB_spring[(9*pas+1):(10*pas)]
abundance_table$SHB_summer[1:15]<-abundance_table$SHB_summer[(10*pas+1):(11*pas)]
abundance_table$SHB_fall[1:15]<-abundance_table$SHB_fall[(11*pas+1):(12*pas)]
abundance_table<-abundance_table[1:15,-1]
abundance_table<-abundance_table[-15,]
total<-c("Total",apply(abundance_table[,-1],2,sum));total
View(abundance_table)
write.xlsx(abundance_table, "Top15_BACT_Families_RelativeAbundance.xlsx")

# FUNGI TABLE FOR MANUSCRIP

# Step 1: Get top 15 families from global dataset
top_families <- FUNGt %>%
  group_by(Family) %>%
  summarise(Global = sum(x, na.rm = TRUE)) %>%
  arrange(desc(Global)) %>%
  slice_head(n = 15)

# Step 2: Define datasets to include
datasets <- list(
  FUN = FUNGt,
  MTL = FUNG_MTLt,
  MTL_spring = FUNG_MTL_springt,
  MTL_summer = FUNG_MTL_summert,
  MTL_fall = FUNG_MTL_fallt,
  QBC = FUNG_QBCt,
  QBC_spring = FUNG_QBC_springt,
  QBC_fall = FUNG_QBC_fallt,
  SHB = FUNG_SHBt,
  SHB_spring = FUNG_SHB_springt,
  SHB_summer = FUNG_SHB_summert,
  SHB_fall = FUNG_SHB_fallt
)

# Step 3: Extract relative abundance for top families from each dataset
abundance_table <- map_dfr(names(datasets), function(name) {
  datasets[[name]] %>%
    group_by(Family) %>%
    summarise(!!name := sum(x, na.rm = TRUE)) %>%
    filter(Family %in% top_families$Family)
}, .id = "source")
pas<-15
abundance_table$MTL[1:15]<-abundance_table$MTL[(pas+1):(2*pas)]
abundance_table$MTL_spring[1:15]<-abundance_table$MTL_spring[(2*pas+1):(3*pas)]
abundance_table$MTL_summer[1:15]<-abundance_table$MTL_summer[(3*pas+1):(4*pas)]
abundance_table$MTL_fall[1:15]<-abundance_table$MTL_fall[(4*pas+1):(5*pas)]
abundance_table$QBC[1:15]<-abundance_table$QBC[(5*pas+1):(6*pas)]
abundance_table$QBC_spring[1:15]<-abundance_table$QBC_spring[(6*pas+1):(7*pas)]
abundance_table$QBC_fall[1:15]<-abundance_table$QBC_fall[(7*pas+1):(8*pas)]
abundance_table$SHB[1:15]<-abundance_table$SHB[(8*pas+1):(9*pas)]
abundance_table$SHB_spring[1:15]<-abundance_table$SHB_spring[(9*pas+1):(10*pas)]
abundance_table$SHB_summer[1:15]<-abundance_table$SHB_summer[(10*pas+1):(11*pas)]
abundance_table$SHB_fall[1:15]<-abundance_table$SHB_fall[(11*pas+1):(12*pas)]
abundance_table<-abundance_table[1:15,-1]
abundance_table<-abundance_table[-15,]
total<-c("Total",apply(abundance_table[,-1],2,sum));total
View(abundance_table)
write.xlsx(abundance_table, "Top15_FUNG_Families_RelativeAbundance.xlsx")

# Step 1: Get top 15 families from global dataset
top_families <- POLLt %>%
  group_by(Family) %>%
  summarise(Global = sum(x, na.rm = TRUE)) %>%
  arrange(desc(Global)) %>%
  slice_head(n = 8)

# Step 2: Define datasets to include
datasets <- list(
  POL = POLLt,
  MTL = POLL_MTLt,
  MTL_spring = POLL_MTL_springt,
  MTL_summer = POLL_MTL_summert,
  MTL_fall = POLL_MTL_fallt,
  QBC = POLL_QBCt,
  QBC_spring = POLL_QBC_springt,
  QBC_fall = POLL_QBC_fallt,
  SHB = POLL_SHBt,
  SHB_spring = POLL_SHB_springt,
  SHB_summer = POLL_SHB_summert,
  SHB_fall = POLL_SHB_fallt
)

# Step 3: Extract relative abundance for top families from each dataset
abundance_table <- map_dfr(names(datasets), function(name) {
  datasets[[name]] %>%
    group_by(Family) %>%
    summarise(!!name := sum(x, na.rm = TRUE)) %>%
    filter(Family %in% top_families$Family)
}, .id = "source")
pas<-8
abundance_table$MTL[1:pas]<-abundance_table$MTL[(pas+1):(2*pas)]
abundance_table$MTL_spring[1:pas]<-abundance_table$MTL_spring[(2*pas+1):(3*pas)]
abundance_table$MTL_summer[1:pas]<-abundance_table$MTL_summer[(3*pas+1):(4*pas)]
abundance_table$MTL_fall[1:pas]<-abundance_table$MTL_fall[(4*pas+1):(5*pas)]
abundance_table$QBC[1:pas]<-abundance_table$QBC[(5*pas+1):(6*pas)]
abundance_table$QBC_spring[1:pas]<-abundance_table$QBC_spring[(6*pas+1):(7*pas)]
abundance_table$QBC_fall[1:pas]<-abundance_table$QBC_fall[(7*pas+1):(8*pas)]
abundance_table$SHB[1:pas]<-abundance_table$SHB[(8*pas+1):(9*pas)]
abundance_table$SHB_spring[1:pas]<-abundance_table$SHB_spring[(9*pas+1):(10*pas)]
abundance_table$SHB_summer[1:pas]<-abundance_table$SHB_summer[(10*pas+1):(11*pas)]
abundance_table$SHB_fall[1:pas]<-abundance_table$SHB_fall[(11*pas+1):(12*pas)]
abundance_table<-abundance_table[1:pas,-1]
abundance_table<-abundance_table[-pas,]
total<-c("Total",apply(abundance_table[,-1],2,sum));total
View(abundance_table)
write.xlsx(abundance_table, "Top8_POLL_Families_RelativeAbundance.xlsx")

