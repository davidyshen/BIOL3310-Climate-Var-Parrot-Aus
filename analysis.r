# Impact of climate variability of Australian parrot range size - Data analysis
# David Shen (2020)
# --------------------------------------------------------------

# Load packages
library(rgdal)
library(rasterVis)
library(sf)
library(ggplot2)
library(ape)
library(nlme)
library(geiger)
library(phytools)

# Set coordinate reference
crdref <- CRS('+proj=longlat +datum=WGS84')
# Load filtered Australia data, exclude first column
data <- read.csv("filtered-ebd.csv")[, -1]
# Replace space with underscore so it works with trees
data$SCIENTIFIC.NAME <- gsub(" ", "_", data$SCIENTIFIC.NAME)
# Vector of species, excluding non-species observations and Trichoglossus rubritorquis - no phylogeny
speciesList <-
  levels(as.factor(data$SCIENTIFIC.NAME))[-c(25, 41:53, 55)]
# Load Australia shape file
aus <- readOGR("AusStateSHP/States Map.shp")
# Phylogenetic tree of species
tree <- read.nexus("birdtree.nex")
# Consensus tree
conTree <- consensus.edges(tree)

# Load climate data and crop to only Australia
climVar <-
  raster::getData(name = "worldclim", res = 2.5, var = "bio")
# Temperature variability
tempVar <- climVar$bio4
tempVar <- crop(tempVar, aus)
# Precipitation variability
rainVar <- climVar$bio15
rainVar <- crop(rainVar, aus)

# Create plots of climate variabilites
tempVarMap <- as.data.frame(rasterToPoints(tempVar))
colnames(tempVarMap) <- c("Longitude", "Latitude", "TS")
TSmap <-
  ggplot(data = tempVarMap, aes(x = Longitude, y = Latitude)) +
  geom_raster(aes(fill = TS)) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_bw()

rainVarMap <- as.data.frame(rasterToPoints(rainVar))
colnames(rainVarMap) <- c("Longitude", "Latitude", "PS")
PSmap <-
  ggplot(data = rainVarMap, aes(x = Longitude, y = Latitude)) +
  geom_raster(aes(fill = PS)) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_bw()

climVarMap <-
  as.data.frame(cbind(
    rainVarMap$Longitude,
    rainVarMap$Latitude,
    sqrt(rainVarMap$PS * tempVarMap$TS)
  ))
colnames(climVarMap) <- c("Longitude", "Latitude", "CS")
CSmap <-
  ggplot(data = climVarMap, aes(x = Longitude, y = Latitude)) +
  geom_raster(aes(fill = CS)) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_bw()

# Identify cells that are land
non <- is.na(tempVar@data@values)
land <- !non
# Vector of cell coords which are land
landCells <- which(land, arr.ind = T)

# Load range size data
response <- read.csv("rangeSize.csv")
response$species <- as.character(response$species)
# Create tables for response variables
responseTV <- response
responsePV <- response
responseCV <- response

# Calcualte slope of species response to climate variability
# Run 1000 times to reduce sampling bias
runs <- 1000
{
  startTime <- Sys.time()
  for (j in 1:runs)
  {
    cat("Run", j, "\n")
    finTable <- c()
    for (i in 1:length(speciesList))
    {
      # Get longitude and latitude points for each species
      temp <- subset(data, SCIENTIFIC.NAME == speciesList[i])
      pts <- cbind(temp$LONGITUDE, temp$LATITUDE)
      # Rasterise points to get occurence cells
      ras <- rasterize(pts, tempVar, proj4string = crdref)
      # Convert to presence absence cells (binary)
      ras@data@values[!is.na(ras@data@values)] <- T
      # Number of occupied cells
      number <- length(na.omit(ras@data@values))
      # Coords of cells with occurences
      occur <- which(ras@data@values == T, arr.ind = T)
      # Coords of cells that are land (landCells) AND NOT occurences (occur)
      nonOccur <- landCells[!landCells %in% occur]
      # Random selection coords of non-occupied cells, size same as number of occupied cells
      nonOc <- sample(nonOccur, size = number)
      # Temp var corresponding to each occupied cell
      ocTV <- tempVar@data@values[occur]
      # Temp var corresponding to each non-occupied sampled cell
      nonTV <- tempVar@data@values[nonOc]
      # Precip var corresponding to each occupied cell
      ocPV <- rainVar@data@values[occur]
      # Precip var corresponding to each non-occupied sampled cell
      nonPV <- rainVar@data@values[nonOc]
      # Presence, absence, species columns
      presAbs <- c(rep(1, times = number), rep(0, times = number))
      sp <- rep(speciesList[i], times = number * 2)
      TV <- c(ocTV, nonTV)
      PV <- c(ocPV, nonPV)
      # Table of cell occupation | species | temp var | rain var
      table <- cbind(presAbs, sp, TV, PV)
      finTable <- rbind(finTable, table)
      cat("Completed", speciesList[i], i, "\n")
    }
    # Compile table of species occupation, species name, temp var, rain var
    finTable <- as.data.frame(finTable)
    finTable$TV <- as.numeric(finTable$TV)
    finTable$PV <- as.numeric(finTable$PV)
    # Calculate climate variation temp x precip
    finTable$CV <- finTable$TV * finTable$PV
    cat("Running model", j, "\n")
    # Run climate variation response model
    modelCV <-
      glm(presAbs ~ sqrt(CV) * sp, data = finTable, family = "binomial")
    coefCV <- as.data.frame(coef(modelCV))
    responseCV[, 2 + j] <- c(coefCV[2, 1], coefCV[42:80, 1])
    intCV <- coefCV[1, 1]
    responseCV[1, 2 + j] <- responseCV[1, 2 + j] + intCV
    int2CV <- responseCV[1, 2 + j]
    responseCV[2:nrow(responseCV), 2 + j] <-
      responseCV[2:nrow(responseCV), 2 + j] + int2CV
    # Temperature variation response model
    modelTV <-
      glm(presAbs ~ TV * sp, data = finTable, family = "binomial")
    coefTV <- as.data.frame(coef(modelTV))
    responseTV[, 2 + j] <- c(coefTV[2, 1], coefTV[42:80, 1])
    intTV <- coefTV[1, 1]
    responseTV[1, 2 + j] <- responseTV[1, 2 + j] + intTV
    int2TV <- responseTV[1, 2 + j]
    responseTV[2:nrow(responseTV), 2 + j] <-
      responseTV[2:nrow(responseTV), 2 + j] + int2TV
    # Precipitation variation response model
    modelPV <-
      glm(presAbs ~ PV * sp, data = finTable, family = "binomial")
    coefPV <- as.data.frame(coef(modelPV))
    responsePV[, 2 + j] <- c(coefPV[2, 1], coefPV[42:80, 1])
    intPV <- coefPV[1, 1]
    responsePV[1, 2 + j] <- responsePV[1, 2 + j] + intPV
    int2PV <- responsePV[1, 2 + j]
    responsePV[2:nrow(responsePV), 2 + j] <-
      responsePV[2:nrow(responsePV), 2 + j] + int2PV
    cat("Completed run",
        j,
        "of",
        runs,
        "\n",
        "******************** \n")
  }
  endTime <- Sys.time()
  cat("Runtime",
      difftime(endTime, startTime, units = "mins"),
      "mins \n")
}

# Calcualte the mean slope of response
responseTV$avResponse <-
  apply(responseTV[, 3:ncol(responseTV)], 1, mean)
responseCV$avResponse <-
  apply(responseCV[, 3:ncol(responseCV)], 1, mean)
responsePV$avResponse <-
  apply(responsePV[, 3:ncol(responsePV)], 1, mean)

# Save slope response to csv
# write.csv(responseTV, "responseTV.csv")
# write.csv(responsePV, "responsePV.csv")
# write.csv(responseCV, "responseCV.csv")

# Load response csv if opening from saved
# responseTV <- read.csv("responseTV.csv")[,-1]
# responsePV <- read.csv("responsePV.csv")[,-1]
# responseCV <- read.csv("responseCV.csv")[,-1]

#Label row names their species to make PGLS easier
row.names(responseCV) <- responseCV[, 1]
row.names(responseTV) <- responseTV[, 1]
row.names(responsePV) <- responsePV[, 1]

# Function to predict response from models
predVal <- function(model, input)
{
  dat <-
    data.frame(response = seq(min(input$avResponse), max(input$avResponse), length = 1000))
  dat$area <- coef(model)[2] * dat$response + coef(model)[1]
  return(dat)
}

# Function to predict average response from model from 100 trees
averagePhy <- function(trees, data)
{
  p <- c()
  slope <- c()
  intercept <- c()
  tvalue <- c()
  for (i in 1:100)
  {
    model <-
      gls(sqrt(area) ~ avResponse,
          data = data,
          correlation = corBrownian(phy = trees[[i]]))
    p[i] <- coef(summary(model))[2, 4]
    slope[i] <- coef(summary(model))[2, 1]
    intercept[i] <- coef(summary(model))[1, 1]
    tvalue[i] <- coef(summary(model))[2, 3]
  }
  return(data.frame(
    pvalue = p,
    slope = slope,
    intercept = intercept,
    tvalue = tvalue
  ))
}

### Species response to temperature variation --------------------------------------------------------------
# Average tree
tempMod3 <- averagePhy(tree, responseTV)
# Calculating predicted values
tempPred3 <-
  data.frame(response = seq(
    min(responseTV$avResponse),
    max(responseTV$avResponse),
    length = 1000
  ))
tempPred3$area <-
  mean(tempMod3$slope) * tempPred3$response + mean(tempMod3$intercept)
# Plotting data
ggplot(responseTV, aes(x = avResponse, y = sqrt(area))) +
  geom_point() +
  geom_line(data = tempPred1,
            aes(x = response, y = area, colour = "No_phylogeny"),
            size = 1) +
  geom_line(data = tempPred2,
            aes(x = response, y = area, colour = "Consensus"),
            size = 1) +
  geom_line(data = tempPred3,
            aes(x = response, y = area, colour = "Tree_average"),
            size = 1) +
  scale_colour_manual(values = c(
    No_phylogeny = "black",
    Consensus = "blue",
    Tree_average = "darkgreen"
  )) +
  geom_label_repel(
    aes(label = species),
    data = subset(responseTV, avResponse < 2.56),
    nudge_y = 700,
    size = 4
  ) +
  geom_label_repel(
    aes(label = species),
    data = subset(responseTV, sqrt(area) > 2500),
    nudge_x = -0.04,
    force = 90,
    size = 4
  ) +
  theme_bw() +
  labs(x = "Response to temperature seasonality", y = expression(paste("Range size ", sqrt("km"))))

## Excluding outliers ----------------------------------
responseTV34 <- responseTV[-c(34, 35), ]
conTree34.2 <-
  drop.tip(conTree,
           c("Psephotus_chrysopterygius", "Psephotus_dissimilis"))
tree34.2 <-
  lapply(tree,
         drop.tip,
         tip = c("Psephotus_chrysopterygius", "Psephotus_dissimilis"))
class(conTree34.2) <- "multiPhylo"
# Average tree
tempMod34.3 <- averagePhy(tree34.2, responseTV34)
# Calculating predicted values
tempPred34.3 <-
  data.frame(response = seq(
    min(responseTV34$avResponse),
    max(responseTV34$avResponse),
    length = 1000
  ))
tempPred34.3$area <-
  mean(tempMod34.3$slope) * tempPred34.3$response + mean(tempMod34.3$intercept)
# Plotting data
tempPlot <-
  ggplot(responseTV34, aes(x = avResponse, y = sqrt(area))) +
  geom_point() +
  geom_line(data = tempPred34.3, aes(x = response, y = area), size = 1) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(x = "Response to temperature seasonality",
       y = expression(paste("Range size ", "sqrt(", km ^ 2, ")")),
       title = "a.")

### Species response to precipitation variation ------------------------------------------------------------
# Average tree
rainMod3 <- averagePhy(tree, responsePV)
# Calculating predicted values
rainPred3 <-
  data.frame(response = seq(
    min(responsePV$avResponse),
    max(responsePV$avResponse),
    length = 1000
  ))
rainPred3$area <-
  mean(rainMod3$slope) * rainPred3$response + mean(rainMod3$intercept)
# Plotting data
ggplot(responsePV, aes(x = avResponse, y = sqrt(area))) +
  geom_point() +
  geom_line(data = rainPred1,
            aes(x = response, y = area, colour = "No_phylogeny"),
            size = 1) +
  geom_line(data = rainPred2,
            aes(x = response, y = area, colour = "Consensus"),
            size = 1) +
  geom_line(data = rainPred3,
            aes(x = response, y = area, colour = "Tree_average"),
            size = 1) +
  scale_colour_manual(values = c(
    No_phylogeny = "black",
    Consensus = "blue",
    Tree_average = "darkgreen"
  )) +
  theme_bw() +
  labs(x = "Response to precipitation seasonality", y = expression(paste("Range size ", sqrt("km"))))

## Excluding outliers ----------------------------------
responsePV34 <- responsePV[-c(34, 35), ]
conTree34.2 <-
  drop.tip(conTree,
           c("Psephotus_chrysopterygius", "Psephotus_dissimilis"))
tree34.2 <-
  lapply(tree,
         drop.tip,
         tip = c("Psephotus_chrysopterygius", "Psephotus_dissimilis"))
class(tree34.2) <- "multiPhylo"
# Average tree
rainMod34.3 <- averagePhy(tree34.2, responsePV34)
# Calculating predicted values
rainPred34.1 <- predVal(rainMod34.1, responsePV34)
rainPred34.2 <- predVal(rainMod34.2, responsePV34)
rainPred34.3 <-
  data.frame(response = seq(
    min(responsePV34$avResponse),
    max(responsePV34$avResponse),
    length = 1000
  ))
rainPred34.3$area <-
  mean(rainMod34.3$slope) * rainPred34.3$response + mean(rainMod34.3$intercept)
# Plotting data
rainPlot <-
  ggplot(responsePV34, aes(x = avResponse, y = sqrt(area))) +
  geom_point() +
  geom_line(data = rainPred34.3, aes(x = response, y = area), size = 1) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(x = "Response to precipitation seasonality",
       y = expression(paste("Range size ", "sqrt(", km ^ 2, ")")),
       title = "b.")

### Species response to combined variation -----------------------------------------------------------------
# Average tree
climMod3 <- averagePhy(tree, responseCV)
# Calculating predicted values
climPred1 <- predVal(climMod1, responseCV)
climPred2 <- predVal(climMod2, responseCV)
climPred3 <-
  data.frame(response = seq(
    min(responseCV$avResponse),
    max(responseCV$avResponse),
    length = 1000
  ))
climPred3$area <-
  mean(climMod3$slope) * climPred3$response + mean(climMod3$intercept)
# Plotting data
ggplot(responseCV, aes(x = avResponse, y = sqrt(area))) +
  geom_point() +
  geom_line(data = climPred3, aes(x = response, y = area), size = 1) +
  theme_bw() +
  labs(x = "Response to climate seasonality", y = expression(paste("Range size ", sqrt("km"))))
## Excluding P. chrysop ----------------------------------------
responseCV34 <- responseCV[-c(34, 35), ]
conTree34.2 <-
  drop.tip(conTree,
           c("Psephotus_chrysopterygius", "Psephotus_dissimilis"))
tree34.2 <-
  lapply(tree,
         drop.tip,
         tip = c("Psephotus_chrysopterygius", "Psephotus_dissimilis"))
class(tree34.2) <- "multiPhylo"
# Average tree
climMod34.3 <- averagePhy(tree34.2, responseCV34)
# Calculating predicted values
climPred34.1 <- predVal(climMod34.1, responseCV34)
climPred34.2 <- predVal(climMod34.2, responseCV34)
climPred34.3 <-
  data.frame(response = seq(
    min(responseCV34$avResponse),
    max(responseCV34$avResponse),
    length = 1000
  ))
climPred34.3$area <-
  mean(climMod34.3$slope) * climPred34.3$response + mean(climMod34.3$intercept)
# Plotting data
climPlot <-
  ggplot(responseCV34, aes(x = avResponse, y = sqrt(area))) +
  geom_point() +
  geom_line(data = climPred34.3, aes(x = response, y = area), size = 1) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(x = "Response to climate seasonality",
       y = expression(paste("Range size ", "sqrt(", km ^ 2, ")")),
       title = "c.")

## Histgram of p values and slopes --------------------------------------------------------------------------
tempHist <- ggplot(aes(x = pvalue), data = tempMod34.3) +
  geom_histogram(
    bins = 20,
    linetype = 1,
    fill = "#FFFFFF",
    colour = "black"
  ) +
  theme_bw() +
  labs(title = "a.")
rainHist <- ggplot(aes(x = pvalue), data = rainMod34.3) +
  geom_histogram(
    bins = 20,
    linetype = 1,
    fill = "#FFFFFF",
    colour = "black"
  ) +
  theme_bw() +
  labs(title = "b.")
climHist <- ggplot(aes(x = pvalue), data = climMod34.3) +
  geom_histogram(
    bins = 20,
    linetype = 1,
    fill = "#FFFFFF",
    colour = "black"
  ) +
  theme_bw() +
  labs(title = "c.")

tempslopeHist <- ggplot(aes(x = slope), data = tempMod34.3) +
  geom_histogram(
    bins = 20,
    linetype = 1,
    fill = "#FFFFFF",
    colour = "black"
  ) +
  theme_bw() +
  labs(title = "a.")
rainslopeHist <- ggplot(aes(x = slope), data = rainMod34.3) +
  geom_histogram(
    bins = 20,
    linetype = 1,
    fill = "#FFFFFF",
    colour = "black"
  ) +
  theme_bw() +
  labs(title = "b.")
climslopeHist <- ggplot(aes(x = slope), data = climMod34.3) +
  geom_histogram(
    bins = 20,
    linetype = 1,
    fill = "#FFFFFF",
    colour = "black"
  ) +
  theme_bw() +
  labs(title = "c.")


#### T-test of slope and p-value mean -----------------------------------------------------------------
t.test(tempMod34.3$slope)
t.test(rainMod34.3$slope)
t.test(climMod34.3$slope)
