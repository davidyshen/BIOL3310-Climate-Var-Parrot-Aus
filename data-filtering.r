# Impact of climate variability of Australian parrot range size - Data filtering
# David Shen (2020)
# --------------------------------------------------------------

# Using auk package to filter data
library(auk)

#Import EBD dataset
datadir <- "D:/eBird"

#Synonyms: 
#Neopsephotus bourkii - Neophema bourkii
#Cacatua leadbeateri - Lophochroa leadbeateri
#Cacatua roseicapilla - Eolophus roseicapilla
#Calyptorhynchus funereus - Zanda funerea
#Zanda latirostris - Calyptorhynchus latirostris
#Zanda baudinii - Calyptorhynchus baudinii
#Neopsephotus bourkii - Neophema bourkii
#Psephotus varius - Psephotellus varius
#Psephotus dissimilis - Psephotellus dissimilis
#Psephotus chrysopterygius - Psephotellus chrysopterygius

#Define species for filtering
speciesName <- c("Polytelis swainsonii",
                 "Polytelis anthopeplus",
                 "Polytelis alexandrae",
                 "Alisterus scapularis",
                 "Pezoporus wallicus",
                 "Neophema bourkii",
                 "Neophema chrysostoma",
                 "Neophema elegans",
                 "Neophema petrophila",
                 "Neophema pulchella",
                 "Neophema splendida",
                 "Lathamus discolor",
                 "Barnardius zonarius",
                 "Platycercus caledonicus",
                 "Platycercus elegans",
                 "Platycercus venustus",
                 "Platycercus adscitus",
                 "Platycercus icterotis",
                 "Northiella haematogaster",
                 "Psephotus haematonotus",
                 "Psephotus varius",
                 "Psephotus dissimilis",
                 "Psephotus chrysopterygius",
                 "Purpureicephalus spurius",
                 "Melopsittacus undulatus",
                 "Glossopsitta concinna",
                 "Glossopsitta pusilla",
                 "Glossopsitta porphyrocephala",
                 "Psitteuteles versicolor",
                 "Trichoglossus chlorolepidotus",
                 "Calyptorhynchus banksii",
                 "Calyptorhynchus lathami",
                 "Calyptorhynchus funereus",
                 "Calyptorhynchus latirostris",
                 "Calyptorhynchus baudinii",
                 "Callocephalon fimbriatum",
                 "Lophochroa leadbeateri",
                 "Cacatua tenuirostris",
                 "Cacatua pastinator",
                 "Nymphicus hollandicus")

#Set paths to data
f_ebd <- paste(datadir, "/ebd_AU/ebd_AU_relMar-2020.txt", sep="")
f_smp <- paste(datadir, "/ebd_sampling/ebd_sampling_relMar-2020.txt", sep="")

#Define filters
#Only filtering by species because data was only from Australia
filters <- auk_ebd(f_ebd, file_sampling = f_smp) %>% 
  auk_species(species = speciesName) %>%
  auk_complete() 

#Create filtered dataset
ebd_filtered <- auk_filter(filters,
                           file = "ebd-filtered.txt",
                           filter_sampling = F,
                           overwrite = T)

filtered_ebd <- paste(datadir, "/ebd_proj/ebd-filtered.txt", sep = "")

#Open species filtered data and select columns for use
data <- read.delim(filtered_ebd)
data <- data[, c("COMMON.NAME", "SCIENTIFIC.NAME", "LATITUDE", "LONGITUDE", "OBSERVATION.DATE")]

#Save filtered and cleaned data
write.csv(data, "filtered-ebd.csv")
