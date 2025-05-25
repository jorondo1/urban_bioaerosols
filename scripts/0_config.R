taxRanks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_16S")

# Variable configuration

kingdoms <- c('BACT' = "Bacteria",
              'FUNG' = "Fungi",
              'PLAN' = "Pollen")

cities <- c('Montreal' = 'Montreal', 
            'Quebec' = 'Quebec', 
            'Sherbrooke' = 'Sherbrooke')

periods <- c('Spring' , 'Summer', 'Fall')

barcodes <- c(
  BACT = '16S',
  FUNG = 'ITS',
  PLAN = 'trnL'
  
)

# Colour variables

period_colours <- c('Spring' = 'springgreen3',
                    'Summer' = 'skyblue3', 
                    'Fall' = 'orange3')
nvdi_colours <- c("#E1BE6A", "#40B0A6", "#E66100")


# Microbial families for community plots:
palettes <- list() 
palettes$BACT <- c(
  Others = "#C6C2C2",
  Unclassified = "#E6E5C1",
  Beijerinckiaceae = "#FFFFFF",
  `67-14` = "hotpink1",
  Streptomycetaceae = "darkmagenta",
  Acetobacteraceae = "brown2",
  Oxalobacteraceae = "sienna2",
  Sphingomonadaceae = "#FFA500",
  Paracoccaceae = "tomato",
  Pseudonocardiaceae = "#FFD700",
  Kineosporiaceae = "lightslateblue",
  Intrasporangiaceae = "mediumvioletred",
  Hymenobacteraceae = "#D700D7",
  Deinococcaceae = "#FF6B6B",
  Microbacteriaceae = "#B2006B", 
  Geodermatophilaceae = "#8A008A",
  Micrococcaceae = "violet",
  Nocardioidaceae = "#6A008A")

palettes$FUNG <- c(
  Others="#C6C2C2",
  Phanerochaetaceae="#FFFFFF",
  Mycosphaerellaceae="paleturquoise1",
  Rickenellaceae="royalblue4",
  Omphalotaceae="lightsteelblue",
  Cortinariaceae="seagreen2",
  Discinellaceae="darkgreen",
  Gloeophyllaceae="palegreen4",
  Atheliaceae="seagreen",
  Strophariaceae="darkslategrey",
  Cerrenaceae="cyan1",
  Incrustoporiaceae="#00CED1",
  Ganodermataceae="#40E0D0",
  Peniophoraceae="#20B2AA",
  Irpicaceae="#008B8B",
  Pleosporaceae="#5F9EA0",
  Fomitopsidaceae="#4682B4",
  Polyporaceae="#1E90FF",
  Cladosporiaceae="deepskyblue4")

palettes$PLAN <- c(
  Others="#C6C2C2",
  Unclassified = "#E6E5C1",
  Poaceae="#32CD32",
  Cupressaceae="#00FF00",
  Betulaceae="#ADFF2F",
  Equisetaceae="#74C476",
  Oleaceae="#238B45",
  Sapindaceae="#41AB5D",
  Asteraceae="#006D2C",
  Pinaceae="#00441B")

#period_colour = c('Sampling time_period 1' = 'springgreen4', 
#                  'Sampling time_period 2' = 'skyblue3', 
#                  'Sampling time_period 3' = 'orange3')
