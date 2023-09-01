#sessionInfo()
#R version 4.2.2 (2022-10-31 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows Server x64 (build 17763)

#_______________________________________________________________________________________

#    --------------------------------- 1. Packages ---------------------------------

#_______________________________________________________________________________________

library(plyr)
library(raster)
library(rgdal)
library(rgeos)
library(ggmap)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggpubr)
library(spatialEco)
library(sf)
library(dplyr)
library(blockCV)
library(spocc)
library(ENMeval)
library(magrittr)
library(MaxentVariableSelection)
library(spThin)
library(data.table)
library(tidyr)
library(biomod2)
library(dismo)
library(readr)
library(forcats)
library(biscale)
library(cowplot)
library(ggfx)
library(parallel)
library(doParallel)

#_______________________________________________________________________________________

#    --------------------------------- 2. Data ---------------------------------

#_______________________________________________________________________________________

species.eu <- c("Andrena.bucephala", "Andrena.chrysosceles", "Andrena.cineraria", 
                "Andrena.dorsata", "Andrena.flavipes", "Andrena.fulva", "Andrena.gravida", 
                "Andrena.haemorrhoa", "Andrena.helvola", "Andrena.minutula", 
                "Andrena.nigroaenea", "Andrena.propinqua", "Andrena.scotica", 
                "Andrena.subopaca", "Andrena.varians", "Anthophora.plumipes", 
                "Bombus.hortorum", "Bombus.hypnorum", "Bombus.jonellus", "Bombus.lapidarius", 
                "Bombus.pascuorum", "Bombus.pratorum", "Colletes.cunicularius", 
                "Eucera.nigrilabris", "Lasioglossum.calceatum", "Lasioglossum.fulvicorne", 
                "Lasioglossum.laticeps", "Lasioglossum.malachurum", "Lasioglossum.marginatum", 
                "Lasioglossum.minutissimum", "Lasioglossum.morio", "Lasioglossum.pauxillum", 
                "Seladonia.tumulora")

tempcol <- colorRampPalette(c("chartreuse4", "yellowgreen", "gold2", 
                              "orange", "firebrick2", "tomato4"))
pal <- colorRampPalette(c("red", "green"))

colramp <- colorRampPalette(c(pal(7)))

#world outline for making maps
world <-
  ne_countries(scale = "medium", returnclass = "sf") #get world map
world.eu <- st_transform(world,
                         "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")

europe <- st_as_sf(readOGR(dsn = './Maps',
                            layer = 'EU_apple_diss'))
#_______________________________________________________________________________________

#    ------------------------------- 3. Model Performance -------------------------

#_______________________________________________________________________________________
enmevalperf <- list.files(
  ".",
  pattern = "mod_perform_enmeval_null_run\\d+\\.csv$",
  full.names = T
)

#create final model performance table -- rbind all files
perf.tab <- rbindlist(lapply(enmevalperf, read.csv), idcol = "origin")
perf.tab[, origin := factor(origin, labels = basename(enmevalperf))]
perf.tab$Run <- as.numeric(gsub("[^\\d]+", "", perf.tab$origin, perl=TRUE))

#only select best models
combined.perf <- data.frame(perf.tab %>% filter(!is.na(Null.AUC.95)))
 
#get list of species and model runs where better than null
goodmodels <- combined.perf[combined.perf$auc.val.avg > combined.perf$Null.AUC.95,]
levels(factor(goodmodels$Species))
setdiff(levels(factor(combined.perf$Species)),levels(factor(goodmodels$Species)))

#list of species and runs
selmod <- paste0(goodmodels$Species,"_run",goodmodels$Run)

#visualise model performance
combined.perf.mean <- goodmodels %>% group_by(Species) %>% 
  summarise_at(names(goodmodels[10]), list(mean = mean, sd = sd)) %>% 
  filter(Species %in% species.eu)

g1 <- ggplot(combined.perf.mean,aes(x=fct_reorder(Species, desc(mean)),y=mean))+ geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) 

g1

#_______________________________________________________________________________________

#    ------------------ 4. Variable Contribution Niche Space----------------------

#_______________________________________________________________________________________
#classify the actual niche of "apple pollinators"

var.contrib <- list.files(
  ".",
  pattern = "var_importance_enmeval_run\\d+\\.csv$",
  full.names = T
)

#convert to single table
var.contrib.tab <- rbindlist(lapply(var.contrib, read.csv), idcol = "origin")
var.contrib.tab[, origin := factor(origin, labels = basename(var.contrib))]
var.contrib.tab$Species <-  sub("\\_.*", "",var.contrib.tab$origin)
var.contrib.tab$Run <- as.numeric(gsub("[^\\d]+", "", var.contrib.tab$origin, perl=TRUE))
var.contrib.tab$SppRun <- paste0(var.contrib.tab$Species,"_run",var.contrib.tab$Run)

#select only good models
var.combi <- data.frame(var.contrib.tab[var.contrib.tab$SppRun %in% selmod,]) %>% 
  filter(Species %in% species.eu) %>% 
  mutate(variable = factor(variable, levels= c("annualPET", 
                                                 "PETseasonality", 
                                                 "bio12", 
                                                 "bio13", 
                                                 "bio14", 
                                                 "bio15", 
                                                 "climaticMoistureIndex", 
                                                 "aridityIndexThornthwaite", 
                                                 "bio1",
                                                 "bio2",
                                                 "bio4", 
                                                 "bio5", 
                                                 "bio6", 
                                                 "embergerQ", 
                                                 "continentality", 
                                                 "thermicityIndex", 
                                                 "growingDegDays0", 
                                                 "growingDegDays5", 
                                                 "LC100_20_25km_EU", 
                                                 "LC100_30_25km_EU", 
                                                 "LC100_40_25km_EU", 
                                                 "LC100_50_25km_EU", 
                                                 "LC100_110_25km_EU", 
                                                 "LC100_120_25km_EU"))) %>% 
  group_by(Species,variable) %>% 
  summarise(permutation.importance=sum(permutation.importance))

#make a new column percentage of total
var.combi1 = var.combi %>%
group_by(Species) %>%
  mutate(freq = round(permutation.importance/sum(permutation.importance), 3)) %>% 
  ungroup() %>%
  arrange(variable, -freq) %>%
  mutate(Species = fct_inorder(Species))

#visualise importance
g2 <- ggplot(var.combi1,aes(x=Species,y=permutation.importance,
                     fill=variable))+
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank())
g2



#_______________________________________________________________________________________

#    ---------------------------------- 5. Range Change------------------------------

#_______________________________________________________________________________________
change <- list.files(
  "./RangeChange",
  pattern = "rangechange_stats_run\\d+\\.csv$",
  full.names = T
)

#Make single table
change.tab <- rbindlist(lapply(change, read.csv), idcol = "origin")
change.tab[, origin := factor(origin, labels = basename(change))]
change.tab$Run <- as.numeric(gsub("[^\\d]+", "", change.tab$origin, perl=TRUE))
change.tab$SppRun <- paste0(change.tab$Binomial,"_run",change.tab$Run)

#select only good models
change.combi <- data.frame(change.tab[change.tab$SppRun %in% selmod,]) %>% 
  filter(Binomial %in% species.eu) %>% 
  filter(!grepl("_41.60", Period))

#new period column
change.full <- change.combi %>% 
          separate(Period, c("Bin", "Pred", "Scenario","RunNum"),
                   sep = "_",remove = F,fill = "left") %>%
          separate(Scenario, c("RCP", "Year1", "Year2"),
                   sep = "\\.",remove = F,fill = "left")

#visualise changes
change.mean <- change.full %>% group_by(Binomial,Scenario,dispersal,Year2) %>% 
  summarise_at(names(change.full[8]), list(mean = mean, sd = sd)) %>% 
  filter(Binomial %in% species.eu) %>% filter(dispersal=="full"&Year2=="80")

#order by change
sel_order <- 
  change.mean %>% 
  filter(Scenario == "rcp85.61.80") %>% 
  arrange(desc(mean)) %>% 
  mutate(Binomial = factor(Binomial))

#visualise change
ggplot(change.mean,aes(x=factor(Binomial, levels = sel_order$Binomial, ordered = TRUE),
                       y=mean,fill=Scenario))+
  geom_bar(stat="identity",position="dodge")+
  facet_grid(.~dispersal)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.background = element_blank())+
  geom_errorbar(aes(ymin=mean - sd, 
                    ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(-90,70),breaks = seq(-100, 1000, by = 20))

#visualise loss
change.loss <- change.full %>% group_by(Binomial,Scenario,dispersal,Year2) %>% 
  summarise_at(names(change.full[6]), list(mean = mean, sd = sd)) %>% 
  filter(Binomial %in% species.eu) %>% filter(dispersal=="full"&Year2=="80")

change.loss$mean = change.loss$mean*(-1)

ggplot(change.loss,aes(x=reorder(Binomial,-mean,sum),y=mean,fill=Scenario))+
  geom_bar(stat="identity",position="dodge")+
  facet_grid(.~dispersal)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank())+
  geom_errorbar(aes(ymin=mean - sd, 
                    ymax=ifelse(mean + sd > 0, 0, mean + sd)), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-90,0),breaks = seq(-100, 0, by = 20))

#visualise gain
change.gain <- change.full %>% group_by(Binomial,Scenario,dispersal,Year2) %>% 
  summarise_at(names(change.full[7]), list(mean = mean, sd = sd)) %>% 
  filter(Binomial %in% species.eu) %>% filter(dispersal=="full"&Year2=="80")


ggplot(change.gain,aes(x=reorder(Binomial,-mean,sum),y=mean,fill=Scenario))+
  geom_bar(stat="identity",position="dodge")+
  facet_grid(.~dispersal)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank())+
  geom_errorbar(aes(ymin=ifelse(mean - sd < 0, 0, mean - sd), 
                    ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,70),breaks = seq(0, 100, by = 20))

#_______________________________________________________________________________________

#    --------------------------------- 6. Country Comparison ----------------------

#_______________________________________________________________________________________

#do countrys load all here
poly.c <-
  st_as_sf(readOGR("./Maps", "European_Countries_Final"))
poly.apple <-
  st_as_sf(readOGR("./Maps", "CountryLevelFAOApples"))


#get list of all change maps and extract relevant data
ras.changes <- list.files("./RangeChange",
                          pattern = ".tif$",
                          full.names = T)

ras.changes2 <- gsub("_binary_","_",ras.changes)


ras.changes3 <- data.frame(ras.changes2) %>%   separate(ras.changes2, c("Drive", "Project", "WP","Folder1","Folder2","Data"),
                                       sep = "/",remove = F,fill = "left") %>% 
                separate(Data, c("Species", "Type", "Dispersal","Pred","ScenarioPeriod","Run","Other"),
                  sep = "_",remove = F,fill = "right")

ras.changes3$Run <- sub(".tif","",ras.changes3$Run)

ras.changes3$selmod <- paste0(ras.changes3$Species,"_",ras.changes3$Run)

ras.changes4 <- ras.changes3[ras.changes3$selmod %in% selmod,]
  
ras.changes5 <- ras.changes4$ras.changes2 

ras.changes5.1 <- ras.changes5[!grepl("full", ras.changes5)]
ras.changes5.2 <- ras.changes5[grepl("full", ras.changes5)]
ras.changes5.3 <- gsub("_prediction","_binary_prediction",ras.changes5.1 )

ras.changes6 <- c(ras.changes5.2,ras.changes5.3)

ras.changes7 <- ras.changes6[grepl("full", ras.changes6)]

ras.changes8 <- ras.changes7[!grepl("_41.60.tif", ras.changes7)]

ras.changes9 <- ras.changes8[grepl("61.80_", ras.changes8)]

#duplciate species list
sp.list = species.eu

# run in parallel the per country statistics
library(parallel)
library(foreach)
library(doParallel)
#setup parallel back-end to use many processors
cores=detectCores(logical = FALSE)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
foreach(i=seq_along(sp.list)) %dopar% {

  library(plyr)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(parallel)
  library(sf)

#for (i in 1:length(sp.list)) {
  TAB = data.frame()#table to fill

  spp <- sp.list[i]

  done.list.run <- list.files(
    "./RangeChange/",
    pattern = paste0(spp,"_regional_range_change_EU.csv"),
    full.names = T
  )

  if(length(done.list.run)>0){print("SKIP")
  }else{
    rasterOptions(tmpdir=file.path("G:/Temp4"))

  spp.changes <- ras.changes9[grepl(spp, ras.changes9)]

  table = data.frame()#table to fill

  for (j in 1:length(spp.changes)) {
    ras.cha <- raster(spp.changes[j])
    plot(ras.cha)
    name_ras <- names(ras.cha)
    period <- sub('.*\\prediction_', '', name_ras)
    period2 <- strsplit(period,"_")[[1]][1]
    hypo <- sub('.*\\_rangechange_', '', name_ras)
    hypo2 <- sub("_.*", "", hypo)
    run <- strsplit(hypo,"_")[[1]][4]

    AT = data.frame()#table to fill

    for (r in 1:length(poly.apple$ADMIN)){

    country <- poly.apple$ADMIN[r]
    reg <- poly.apple[poly.apple$ADMIN == country, ]#subset by site

    reg_sp <- as(reg, 'Spatial')

    #create outside buffer
    b_large <- gBuffer(reg_sp, width = 50000, quadsegs = 100)

    #crop with large buffer
    ras2 <- crop(ras.cha, extent(b_large))#crop to large buffer
    ras3 <-
      raster::mask(ras2, b_large)#extract raster from within large buffer

    #turn polygon into raster, so we do not lose edges
    SpP_ras <-
      rasterize(reg, ras3, getCover = TRUE)#rasterisize the buffer
    SpP_ras[SpP_ras == 0] <- NA #remove everything outside of buffer


    #get data
    ras4 = stack(raster::mask(ras3, SpP_ras))#extract data from within buffer

    #simplify for clipping
    reg.clip <- rgeos::gSimplify(reg_sp, tol = 10)

    #convert raster to polgyon to get exact area
    ras_vec <-
      rasterToPolygons(
        ras4,
        fun = NULL,
        n = 16,
        na.rm = TRUE,
        digits = 12,
        dissolve = FALSE
      )

    clip <-
      raster::intersect(ras_vec, reg.clip)#get only vector within the buffer

    #plot(clip)

    #get percentage cover by area
    clip@data$Area_m2 <-
      raster::area(clip) #make new column of area of each polygon
    totalarea = raster::area(reg.clip) #total area of buffer

    #calculate proportion of area
    clip_df <- data.frame(clip)
    areatab <-
      clip_df %>% group_by(get(name_ras)) %>% summarise(Area_m2 = sum(Area_m2)) #sum by land use type
    areatab$prop = areatab$Area_m2 / totalarea

    #make table of area - addin necessary columns
    at <- data.frame(areatab)
    at$Period <- c(period2)
    at$Binomial <- spp
    at$Country <- c(country)
    at$Dispersal <- c(hypo2)
    at$Run <- c(run)

    AT = rbind(AT, at)#join
    print(paste(country, "Complete"))
  }
  do.call(file.remove, list(list.files("G:/Temp/TEMP3", full.names = TRUE)))#removes temp files
  table = rbind(table, AT)#join
  print(paste(period, "Complete"))
  }

#Create final table
TAB <- rbind(TAB, table)#join

TAB$Period.x <- sub('.*\\prediction_', '', TAB$Period)

change.df <- data.frame(Change=c("Lost","Stable_Present","Stable_Absent","Gained"),
                        Code=c(-2,-1,0,1))

TAB2 <- merge(TAB,change.df,by.x="get.name_ras.",by.y="Code")

write.csv(TAB2,paste0("./RangeChange/",spp,"_regional_range_change_EU.csv"),row.names=F)

print(paste(spp, "Complete"))

  }
}
 stopCluster(cl)

#join all tables together
change.country <- list.files(
  "./RangeChange",
  pattern = "_regional_range_change_EU.csv$",
  full.names = T
)

#final  table
change.c.tab <- unique(rbindlist(lapply(change.country, read.csv), idcol = "origin"))

#mean proportional species range change per country 2080 per scenario
change.c <- change.c.tab %>% filter(Period %like% "61.80") %>% 
  filter(Change=="Lost"|Change=="Gained") %>% 
  dplyr::select(prop,Period,Binomial,Country, Change, Run) %>% 
  unique() %>% 
  pivot_wider(names_from = Change,values_from = prop, values_fill = 0)

change.c$RangeChange <- change.c$Gained-change.c$Lost

change.c.prop <- change.c %>% group_by(Country,Period) %>% 
  summarise_at(names(change.c[7]), list(mean = mean, sd = sd)) %>% 
  pivot_wider(names_from = Period,values_from = c(mean,sd),names_prefix = "SRC_")

poly.c.change <- merge(poly.apple, change.c.prop, by.x="ADMIN",by.y="Country")

loss.c.prop <- change.c %>% group_by(Country,Period) %>% 
  summarise_at(names(change.c[5]), list(mean = mean, sd = sd)) %>% 
  pivot_wider(names_from = Period,values_from = c(mean,sd),names_prefix = "Loss_")

poly.c.change.l <- merge(poly.c.change, loss.c.prop, by.x="ADMIN",by.y="Country")

#make bivariate maps
#Yield/Loss Map
data <- bi_class(poly.c.change.l, x = Yld.h.., y = mean_Loss_rcp85.61.80, style = "fisher", dim = 3)
map <- ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "DkBlue", dim = 3) +
  labs(
  ) +
  bi_theme()
legend <- bi_legend(pal = "DkBlue",
                    dim = 3,
                    xlab = "Higher Yield (hg/ha)",
                    ylab = "Higher Range Loss (%)",
                    size = 8)
finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.2, .65, 0.2, 0.2)

finalPlot

#Production/Loss Map
data <- bi_class(poly.c.change.l, x = Prdct.., y = mean_Loss_rcp85.61.80, style = "fisher", dim = 3)
map <- ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "DkCyan", dim = 3) +
  labs(
  ) +
  bi_theme()
legend <- bi_legend(pal = "DkCyan",
                    dim = 3,
                    xlab = "Higher Production (t)",
                    ylab = "Higher Range Loss (%)",
                    size = 8)
finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.2, .65, 0.2, 0.2)

finalPlot

#_______________________________________________________________________________________

#    ------------------------------ 7. Range Overlap -----------------------

#_______________________________________________________________________________________

#make consensus apple maps per scenario

#load binary rasters
apple.rasters1 <- list.files("./Binary",
                             pattern = ".tif$",
                             full.names = T)

#select out apple rasters
apple.rasters.binary <- apple.rasters1[grepl("apples", apple.rasters1)]

good.apples <- c("run1","run4","run7","run8","run9","run10")

apple.rasters.binary2 <- (grep(paste(good.apples,collapse="|"), 
      apple.rasters.binary, value=TRUE))

scenario.list.pres <- c("present","rcp26.41.60","rcp45.41.60","rcp85.41.60",
                        "rcp26.61.80","rcp45.61.80","rcp85.61.80")

#create consensus rasters of apple predictions
for(sc in 1:length(scenario.list.pres)){
  sl <- scenario.list.pres[sc]
  arb <- apple.rasters.binary2[grepl(sl, apple.rasters.binary2)]
  arb.agri <- arb[grepl("prediction_agri", arb)]
  arb.agri <- arb.agri[!grepl("limited", arb.agri)]
  arb.lim <- arb[grepl("prediction_limited", arb)]

  arb.agri.sum <- sum(stack(arb.agri))
  arb.agri.sum[arb.agri.sum < 6] = 0
  arb.agri.sum[arb.agri.sum >= 6] = 1
  plot(arb.agri.sum)

  name.agri <- gsub("run1","consensus",arb.agri[1])
  name.agri2 <- gsub("Binary.","Consensus.",name.agri)
  writeRaster(arb.agri.sum,paste0(name.agri2))

  arb.lim.sum <- sum(stack(arb.lim))
  arb.lim.sum[arb.lim.sum < 6] = 0
  arb.lim.sum[arb.lim.sum >= 6] = 1
  plot(arb.lim.sum)

  name.lim <- gsub("run1","consensus",arb.lim[1])
  name.lim2 <- gsub("Binary.","Consensus.",name.lim)
  writeRaster(arb.lim.sum,paste0(name.lim2))
}

#load consensus rasters
apple.rasters2 <- list.files("./Consensus",
                             pattern = ".tif$",
                             full.names = T)
apple.rasters.con <- apple.rasters2[grepl("apples", apple.rasters2)]


# load all binary sepcies rasterts
binary.rasters2 <- list.files("./Binary/All.Binary",
                             pattern = ".tif$",
                             full.names = T)

scenario.list <- c("rcp26.41.60","rcp45.41.60","rcp85.41.60",
                   "rcp26.61.80","rcp45.61.80","rcp85.61.80")

#relevant when multiple dispersals used - not in final model
dispersal.list <- c("_full_")

runlist <- paste0("run",1:10,".tif")

apple.dispersal.list <- c("agri","limited")

#calculate overlap between species and apple distributions

#setup parallel back-end to use many processors
cores=detectCores(logical = FALSE)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
foreach(sp=seq_along(sp.list)) %dopar% {

  library(plyr)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(parallel)
  library(sf)
  library(ggplot2)

  #for (sp in 1:length(sp.list)) {
    spp <- sp.list[sp]

    kap <- binary.rasters2[grepl(spp, binary.rasters2)]

    for (s in 1:length(scenario.list.pres)) {
      scenario <- scenario.list.pres[s]

      kap2 <-  kap[grepl(scenario, kap)]

      for (d in 1:length(dispersal.list)) {
        dispersal <- dispersal.list[d]

        kap3 <-  kap2[grepl(dispersal, kap2)]

        for (run in 1:length(kap3)) {

          kap4 <-  kap3[run]

          ras.bee <- raster(kap4)

          runnum <- strsplit(kap4,"_")[[1]][7]

          for (a in 1:2) {
            apple.dispersal <- apple.dispersal.list[a]

            done.list.good <- list.files(
              ".U/Overlap/",
              pattern = paste0(
                spp,
                "_",
                scenario,
                dispersal,
                apple.dispersal,
                "_",
                gsub(".tif", "", runnum),
                ".csv"
              ),
              full.names = T
            )

            if (length(done.list.good) > 0) {
              print("SKIP")
            } else{

            if (a == 1) {
                apple.ras.d <-
                  apple.rasters.con[!grepl("limited", apple.rasters.con)]
              } else{
                apple.ras.d <-
                apple.rasters.con[grepl("limited", apple.rasters.con)]
            }



          apple.ras.d2 <- apple.ras.d[grepl(scenario,apple.ras.d)]

          ras.apple <- raster(apple.ras.d2)

          ras.bee[ras.bee == 1] = 2 #reclassify to make distinction between apple
          kap.apple <-
            ras.bee + ras.apple #1=apple only, 2=species only, 3=overlap

          # #grids and area of species and apple and overlap
          spp.grids <- sum(na.omit(kap.apple[]) >= 2)
          spp.area <- spp.grids * (res(kap.apple)[1] / 1000) ^ 2 #area in km2
          apple.grids <-
            sum(na.omit(kap.apple[]) == 1 | na.omit(kap.apple[]) == 3)
          apple.area <-
            apple.grids * (res(kap.apple)[1] / 1000) ^ 2 #area in km2

          overlap.grids <- sum(na.omit(kap.apple[]) == 3)
          overlap.area <-
            overlap.grids * (res(kap.apple)[1] / 1000) ^ 2 #area in km2

          #what percentage/area of species range within apple?
          prop.spp.range.in.apple <- overlap.grids / spp.grids

          #what percentage/area of apple range occupied by species?
          prop.apple.range.with.spp <- overlap.grids / apple.grids

          #create a polygon of overlap
          kap.apple.poly <- kap.apple
          kap.apple.poly[kap.apple.poly < 3] = NA

          if(!is.na(cellStats(kap.apple.poly, "mean"))){
            overlap.poly <-  st_as_sf(
              rasterToPolygons(
                kap.apple.poly,
                fun = NULL,
                n = 4,
                na.rm = T,
                digits = 12,
                dissolve = F
              )
            )
          }

          #save a map
          g3 <- ggplot() +
            geom_sf(data = world.eu, colour = 'black') +
            geom_tile(data = data.frame(rasterToPoints(
            kap.apple
            )),
            aes(
              x = x,
              y = y,
              fill = factor(layer)
            )) +
            {if(!is.na(cellStats(kap.apple.poly, "mean")))geom_sf(
              data = overlap.poly,
              fill = NA,
              size = 0.01,
              colour = 'white',
              alpha = 0.5
            )} +
            #geom_sf(data = europe, fill= NA, size = 0.1, colour='black') +
            #geom_sf(data = morocco, fill= NA, size = 0.1, colour='black') +
            coord_sf(
              xlim = c(1878430, 5935415),
              ylim = c(1450620, 5424004),
              expand = FALSE
            ) +
            theme(
              plot.title = element_text(size = 12, hjust = 0.5),
              panel.grid.major = element_line(
                color = gray(.5),
                linetype = "dashed",
                size = 0.5
              ),
              panel.background = element_rect(fill = "aliceblue"),
              plot.caption = element_text(hjust = 0.5, size = rel(0.6)),
              legend.position = c(0.15, 0.75),
              legend.title = element_blank()
            ) +
            labs(
              x = "Longitude",
              y = "Latitude",
              title = paste0(spp, " ", scenario, " ",gsub("_","",dispersal), " ",gsub(".tif","",runnum)," ",apple.dispersal),
              caption = "Source : WorldClim"
            ) +
            scale_fill_manual(
              labels = c("Absent", "Apple Range", "Species Range", "Overlap"),
              values = c("grey", "#7C1124", "#0C2CAA", "#652597")
            )


          #save map here #e.g. pdf
          print(g3)
          

          #build table
          df <- data.frame(
            Binomial = spp,
            Dispersal = gsub("_","",dispersal),
            Scenario = scenario,
            Run = gsub(".tif", "", runnum),
            SpeciesGrids = spp.grids,
            SpeciesArea = spp.area,
            AppleGrids = apple.grids,
            AppleArea = apple.area,
            OverlapGrids = overlap.grids,
            OverlapArea = overlap.area,
            Prop.in.Apple = prop.spp.range.in.apple,
            Prop.of.Apple = prop.apple.range.with.spp,
            AppleDispersal = apple.dispersal
          )
          write.csv(df,paste0("./Overlap/",spp,"_",
                    scenario,dispersal,apple.dispersal,"_",gsub(".tif","",runnum),".csv"),
                    row.names=F)
          print(apple.dispersal)
        }
        print(dispersal)
      }
      print(scenario)

    }
    print(runnum)
  }
  }
  print(spp)
}
stopCluster(cl)


#get all csvs and join them
overlap.csv <- list.files(
  "./Overlap",
  pattern = ".csv$",
  full.names = T
)

#make single table
overlap.tab <- rbindlist(lapply(overlap.csv, read.csv))

#visualise
mean.per.species <- overlap.tab %>% group_by(Scenario,AppleDispersal,Dispersal,Binomial) %>% 
  summarise_at(names(overlap.tab)[11:12], list(mean = mean, sd = sd)) %>% 
  filter(AppleDispersal!="nod"&AppleDispersal!="int"&Dispersal!="nod"&AppleDispersal!="limited"&Dispersal!="int") %>% 
  filter(Scenario=="present"|Scenario=="rcp26.41.60"|Scenario=="rcp26.61.80")

#calculate per species percentage change between present and 2080
mean.pres <- mean.per.species %>% filter(Scenario=="present")
mean.2080 <- mean.per.species %>% filter(Scenario=="rcp26.61.80")

percchange <- data.frame(Binomial=mean.pres$Binomial,
           PercChange=mean.2080$Prop.in.Apple_mean-mean.pres$Prop.in.Apple_mean)
percchange <- percchange %>% mutate(Group = case_when(PercChange <= (-0.10) ~ "Decrease",
                                                      PercChange > (-0.10) ~ "Stable",
                                                      PercChange >= 0.1 ~ "Increase"))

mean.per.species2 <- merge(mean.per.species,percchange, by="Binomial")
library(ggrepel)

mean.per.species2 %>% 
  mutate(label = if_else(Scenario=="rcp26.61.80", as.character(Binomial), NA_character_)) %>% 
  ggplot(aes(y = Prop.in.Apple_mean, x=Scenario, group=Binomial,color=factor(Group))) + 
  geom_point(shape=1) +
  geom_line(aes(color=factor(Group))) +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  geom_text_repel(aes(label = gsub("^.*$", " ", label)), # This will force the correct position of the link's right end.
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE,
  ) +
  geom_text_repel(data = . %>% filter(!is.na(label)),
                  aes(label = paste0("  ", label)),
                  segment.alpha = 0, ## This will 'hide' the link
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  # segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE)+
  scale_x_discrete(expand = expansion(mult = c(0, .5)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1))


mean.per.species <- overlap.tab %>% group_by(Scenario,AppleDispersal,Dispersal,Binomial) %>% 
  summarise_at(names(overlap.tab)[11:12], list(mean = mean, sd = sd)) %>% 
  filter(AppleDispersal!="nod"&AppleDispersal!="int"&Dispersal!="nod"&AppleDispersal!="limited"&Dispersal!="int") %>% 
  filter(Scenario=="present"|Scenario=="rcp45.41.60"|Scenario=="rcp45.61.80")
#calculate per species percentage change between present and 2080
mean.pres <- mean.per.species %>% filter(Scenario=="present")
mean.2080 <- mean.per.species %>% filter(Scenario=="rcp45.61.80")

percchange <- data.frame(Binomial=mean.pres$Binomial,
                         PercChange=mean.2080$Prop.in.Apple_mean-mean.pres$Prop.in.Apple_mean)
percchange <- percchange %>% mutate(Group = case_when(PercChange <= (-0.10) ~ "Decrease",
                                                      PercChange > (-0.10) ~ "Stable",
                                                      PercChange >= 0.1 ~ "Increase"))

mean.per.species2 <- merge(mean.per.species,percchange, by="Binomial")

mean.per.species2 %>% 
  mutate(label = if_else(Scenario=="rcp45.61.80", as.character(Binomial), NA_character_)) %>% 
  ggplot(aes(y = Prop.in.Apple_mean, x=Scenario,  group=Binomial,color=factor(Group))) + 
  geom_point(shape=1) +
  geom_line(aes(color=factor(Group))) +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  geom_text_repel(aes(label = gsub("^.*$", " ", label)), # This will force the correct position of the link's right end.
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE,
  ) +
  geom_text_repel(data = . %>% filter(!is.na(label)),
                  aes(label = paste0("  ", label)),
                  segment.alpha = 0, ## This will 'hide' the link
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  # segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE)+
  scale_x_discrete(expand = expansion(mult = c(0, .5)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1))



mean.per.species <- overlap.tab %>% group_by(Scenario,AppleDispersal,Dispersal,Binomial) %>% 
  summarise_at(names(overlap.tab)[11:12], list(mean = mean, sd = sd)) %>% 
  filter(AppleDispersal!="nod"&AppleDispersal!="int"&Dispersal!="nod"&AppleDispersal!="limited"&Dispersal!="int") %>% 
  filter(Scenario=="present"|Scenario=="rcp85.41.60"|Scenario=="rcp85.61.80")

#calculate per species percentage change between present and 2080
mean.pres <- mean.per.species %>% filter(Scenario=="present")
mean.2080 <- mean.per.species %>% filter(Scenario=="rcp85.61.80")

percchange <- data.frame(Binomial=mean.pres$Binomial,
                         PercChange=mean.2080$Prop.in.Apple_mean-mean.pres$Prop.in.Apple_mean)
percchange <- percchange %>% mutate(Group = case_when(PercChange <= (-0.10) ~ "Decrease",
                                                      PercChange > (-0.10) ~ "Stable",
                                                      PercChange >= 0.1 ~ "Increase"))

mean.per.species2 <- merge(mean.per.species,percchange, by="Binomial")

mean.per.species2 %>% 
  mutate(label = if_else(Scenario=="rcp85.61.80", as.character(Binomial), NA_character_)) %>% 
  ggplot(aes(y = Prop.in.Apple_mean, x=Scenario, group=Binomial,color=factor(Group))) + 
  geom_point(shape=1) +
  geom_line(aes(color=factor(Group))) +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  geom_text_repel(aes(label = gsub("^.*$", " ", label)), # This will force the correct position of the link's right end.
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE,
  ) +
  geom_text_repel(data = . %>% filter(!is.na(label)),
                  aes(label = paste0("  ", label)),
                  segment.alpha = 0, ## This will 'hide' the link
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  # segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE)+
  scale_x_discrete(expand = expansion(mult = c(0, .5)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1))



library(ggrepel)

mean.per.species %>% 
  mutate(label = if_else(Scenario=="rcp85.61.80", as.character(Binomial), NA_character_)) %>% 
  ggplot(aes(y = Prop.of.Apple_mean, x=Scenario, group=Binomial)) + 
  geom_point(shape=1) +
  geom_line() +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  geom_text_repel(aes(label = gsub("^.*$", " ", label)), # This will force the correct position of the link's right end.
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE,
  ) +
  geom_text_repel(data = . %>% filter(!is.na(label)),
                  aes(label = paste0("  ", label)),
                  segment.alpha = 0, ## This will 'hide' the link
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  # segment.color = 'grey',
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 0.15,
                  force = 0.5,
                  hjust = 0,
                  direction="y",
                  na.rm = TRUE)+
  scale_x_discrete(expand = expansion(mult = c(0, .5)))



#_______________________________________________________________________________________

#    --------------------------------- 8. Richness Maps ----------------------

#_______________________________________________________________________________________

#mean change in species richness per cell per country
tempcol <- colorRampPalette(c("grey", "blue", "turquoise3",
                              "chartreuse4", "yellowgreen", "gold2",
                              "orange", "firebrick2", "tomato4"))

kap.rasters.binary <- ras.changes6[!grepl("_41.60.tif", ras.changes6)]
kap.rasters.binary.rich <- kap.rasters.binary[!grepl("apples_", kap.rasters.binary)]

prd.list <- c("rcp26.41.60","rcp45.41.60","rcp85.41.60",
                 "rcp26.61.80","rcp45.61.80","rcp85.61.80")
plot.list <- list()

#make a species loop first where we get a consensus map per scenario per dispersal
for (sp in 1:length(sp.list)) {

  spp <- sp.list[sp]
  kap.ras.spp <-
    kap.rasters.binary.rich[grepl(spp, kap.rasters.binary.rich)]

  for (d in 1:length(dispersal.list)) {
    dispersal <- dispersal.list[d]
    kap.ras.d <-
      kap.ras.spp[grepl(dispersal, kap.ras.spp)]

    for (s in 1:length(scenario.list)) {
      scenario <- scenario.list[s]

      kap.ras.s <-
        kap.ras.d[grepl(scenario, kap.ras.d)]

      kap.ras.spp2 <- stack(kap.ras.s)

      if (s == 1) {
        #do present as well
        kap.pres <- kap.ras.spp2
        kap.pres[kap.pres == 1] = 0
        kap.pres[kap.pres == -1] = 1
        kap.pres[kap.pres == -2] = 1

        d2 <- unstack(kap.pres)

        for(i in seq_along(d2)){
        file <- paste0("/Binary/All.Binary/",names(d2[[i]]),".tif")
        file1 <- gsub("_binary_", "_", file)
        file2 <- gsub("rangechange", "binary", file1)
        file3 <- gsub(scenario, "present", file2)
        writeRaster(d2[[i]], file=file3,overwrite=T)}

        spp.full.sum <- sum(kap.pres)
        spp.full.sum[spp.full.sum < 10] = 0
        spp.full.sum[spp.full.sum >= 10] = 1

        name.full <- gsub("run1", "consensus", kap.ras.s[1])
        name.full1 <- gsub("_binary_", "_", name.full)
        name.full2 <- gsub("RangeChange.", "Consensus.", name.full1)
        name.full3 <- gsub(scenario, "present", name.full2)
        name.full4 <- gsub("_rangechange_", "_binary_", name.full3)
        writeRaster(spp.full.sum, paste0(name.full4),overwrite=T)

        kap.sc <- kap.ras.spp2
        kap.sc[kap.sc == -1] = 1
        kap.sc[kap.sc == -2] = 0

        d3 <- unstack(kap.sc)
        for(i in seq_along(d3)){
          file <- paste0("/Binary/All.Binary/",names(d3[[i]]),".tif")
          file1 <- gsub("_binary_", "_", file)
          file2 <- gsub("rangechange", "binary", file1)
          writeRaster(d3[[i]], file=file2,overwrite=T)}

        spp.full.sum2 <- sum(kap.sc)
        spp.full.sum2[spp.full.sum2 < 10] = 0
        spp.full.sum2[spp.full.sum2 >= 10] = 1

        name.full.1 <- gsub("run1", "consensus", kap.ras.s[1])
        name.full.1.1 <- gsub("_binary_", "_", name.full.1)
        name.full.2 <- gsub("RangeChange.", "Consensus.", name.full.1.1)
        name.full.3 <- gsub("_rangechange_", "_binary_", name.full.2)
        writeRaster(spp.full.sum2, paste0(name.full.3),overwrite=T)
      }
      if (s != 1) {
        kap.sc <- kap.ras.spp2
        kap.sc[kap.sc == -1] = 1
        kap.sc[kap.sc == -2] = 0

        d3 <- unstack(kap.sc)
        for(i in seq_along(d3)){
          file <- paste0("./Binary/All.Binary/",names(d3[[i]]),".tif")
          file1 <- gsub("_binary_", "_", file)
          file2 <- gsub("rangechange", "binary", file1)
          writeRaster(d3[[i]], file=file2,overwrite=T)}

        spp.full.sum2 <- sum(kap.sc)
        spp.full.sum2[spp.full.sum2 < 10] = 0
        spp.full.sum2[spp.full.sum2 >= 10] = 1

        name.full <- gsub("run1", "consensus", kap.ras.s[1])
        name.full1 <- gsub("_binary_", "_", name.full)
        name.full2 <- gsub("RangeChange.", "Consensus.", name.full1)
        name.full3 <- gsub("_rangechange_", "_binary_", name.full2)
        writeRaster(spp.full.sum2, paste0(name.full3),overwrite=T)
      }
      print(paste(spp,scenario,dispersal))
    }
  }
}

#make richness maps
consensus.rasters <- list.files("./Consensus",
                              pattern = ".tif$",
                              full.names = T)

consensus.rasters.bees <- consensus.rasters[!grepl("apples",consensus.rasters)]

for (d in 1:length(dispersal.list)) {
  dispersal <- dispersal.list[d]
  kap.ras.d <-
    consensus.rasters.bees[grepl(dispersal, consensus.rasters.bees)]

  for (s in 1:length(scenario.list.pres)) {
    scenario <- scenario.list.pres[s]

    kap.ras.s <-
      kap.ras.d[grepl(scenario, kap.ras.d)]

    rich.map <- sum(stack(kap.ras.s))

    name <- kap.ras.s[1]
    writeRaster(rich.map,paste0("./Richness.Maps/binary"
                                ,dispersal,
                                "richness_consensus_",
                                scenario,".tif"))
  }
}


#visualize them
richness.rasters <- list.files("Richness.Maps",
                                pattern = ".tif$",
                                full.names = T)
consensus.rasters.apples <- consensus.rasters[grepl("apples",consensus.rasters)]

library(ggfx)
greens <- colorRampPalette(c("#E9F7EF", "#D4EFDF", "#A9DFBF",
  "#A9DFBF", "#7DCEA0", "#52BE80","#27AE60",
                               "#229954", "#1E8449", "#196F3D", "#145A32"))
world <-
  ne_countries(scale = "medium", returnclass = "sf") #get world map
world.eu <- st_transform(world,
                         "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")

world.eu2 <- st_crop(world.eu,extent(2635415,6526354,1385618,5424004))


for (d in 1:length(dispersal.list)) {
  dispersal <- dispersal.list[d]
  kap.ras.d <-
    richness.rasters[grepl(dispersal, richness.rasters)]
  max.val <- max(maxValue(stack( kap.ras.d)))
  plot.list=list()
for (i in 1:length(scenario.list.pres)){
  prd <- scenario.list.pres[i]
  apple.ras.s <- consensus.rasters.apples[grepl(scenario,consensus.rasters.apples)]
  
  ras.binary.eu.sel <- 
      raster(kap.ras.d[grepl(prd, kap.ras.d)])
    names(ras.binary.eu.sel) <- "layer"
    
    pa.eu <- mask(crop(ras.binary.eu.sel,europe),europe)
    extent(europe)
    
    p.p <- ggplot() +
   geom_sf(data = world.eu2,colour='black') +
   geom_raster(data = data.frame(rasterToPoints(pa.eu)), 
                              aes(x = x, y = y, fill = layer)) +
  geom_sf(data = europe2,fill='NA',lwd=0.5,colour='black') +
      # geom_path(aes(x = long, y = lat, group = group), data = outline.af, 
      #           size=1, col="white")+
  coord_sf(xlim = c(2635415,6526354),
           ylim = c(1385618,5424004),
           expand = FALSE) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.major = element_line(
      color = gray(.5),
      linetype = "dashed",
      size = 0.5
    ),
    panel.background = element_rect(fill = "aliceblue"),
    plot.caption = element_text(hjust = 0.5, size = rel(0.6)),
    legend.position = c(0.9, 0.8),
    legend.title = element_blank(),
    legend.background = element_blank()
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title=paste0(prd," ")
  )+
  scale_fill_gradientn(colors = greens(100),
                       limits = c(0, max.val))
p.p

plot.list[[i]] <- p.p

}

print(ggarrange(plot.list[[1]]+ rremove("legend"),
                labels = "a",ncol= 2,heights=c(0.5,1),
          ggarrange(plot.list[[2]],plot.list[[3]],
                    plot.list[[4]],plot.list[[5]],
                    plot.list[[6]],plot.list[[7]],
                    col = 3,common.legend = TRUE,
                    labels = c("b","c","d","e","f","g"))))


}
scenario.list <- c("rcp26.41.60","rcp45.41.60","rcp85.41.60",
                   "rcp26.61.80","rcp45.61.80","rcp85.61.80")

for (d in 1:length(dispersal.list)) {
  dispersal <- dispersal.list[d]
  kap.ras.d <-
    richness.rasters[grepl(dispersal, richness.rasters)]
  max.val <- max(maxValue(stack( kap.ras.d)))
  plot.list=list()
 
  pres.rich <- raster(kap.ras.d[1])
  names(pres.rich) <- "layer"
   for (i in 1:length(scenario.list)){
    prd <- scenario.list[i]
     
    fut.rich <- 
      raster(kap.ras.d[grepl(prd, kap.ras.d)])
    names(fut.rich) <- "layer"
    
    fut.change <- fut.rich-pres.rich
    
    pa.eu.pres <- mask(crop(pres.rich,europe2),europe2)
    pa.eu <- mask(crop(fut.change,europe2),europe2)
    extent(europe2)
   
    #max gain and max loss of species - manual input
    if(dispersal=="_full_"){maxv = 5;minv=-18}
    
    p.p.pres <- ggplot() +
      geom_sf(data = world.eu2,colour='black') +
      geom_raster(data = data.frame(rasterToPoints(pa.eu.pres)), 
                  aes(x = x, y = y, fill = layer)) +
      geom_sf(data = europe2,fill='NA',lwd=0.5,colour='black') +
      # geom_path(aes(x = long, y = lat, group = group), data = outline.af, 
      #           size=1, col="white")+
      coord_sf(xlim = c(2635415,6526354),
               ylim = c(1385618,5424004),
               expand = FALSE) +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        panel.grid.major = element_line(
          color = gray(.5),
          linetype = "dashed",
          size = 0.5
        ),
        panel.background = element_rect(fill = "aliceblue"),
        plot.caption = element_text(hjust = 0.5, size = rel(0.6)),
        legend.position = c(0.9, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank()
      ) +
      labs(
        x = "Longitude",
        y = "Latitude",
        title=paste0(prd," ")
      )+
      scale_fill_gradientn(colors = greens(100),
                           limits = c(0, max.val))
    p.p <- ggplot() +
      geom_sf(data = world.eu2,colour='black') +
      geom_raster(data = data.frame(rasterToPoints(pa.eu)), 
                  aes(x = x, y = y, fill = layer)) +
      geom_sf(data = europe2,fill='NA',lwd=0.5,colour='black') +
      # geom_path(aes(x = long, y = lat, group = group), data = outline.af, 
      #           size=1, col="white")+
      coord_sf(xlim = c(2635415,6526354),
               ylim = c(1385618,5424004),
               expand = FALSE) +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        panel.grid.major = element_line(
          color = gray(.5),
          linetype = "dashed",
          size = 0.5
        ),
        panel.background = element_rect(fill = "aliceblue"),
        plot.caption = element_text(hjust = 0.5, size = rel(0.6)),
        legend.position = c(0.9, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank()
      ) +
      labs(
        x = "Longitude",
        y = "Latitude",
        title=paste0(prd," ")
      )+
      scale_fill_gradientn(colours = c("#641E16","#C0392B","#F9EBEA", "white", "#EBF5FB", "#3498DB","#1B4F72"), 
                           values = scales::rescale(c(minv,(minv+5),(minv+10),-.1,0,.1,maxv-0.8,maxv-0.5,maxv)),
                           limits = c(minv, maxv))
    p.p
    
    if(i==1){
    plot.list[[i]] <- p.p.pres
    plot.list[[i+1]] <- p.p
    }else{
    plot.list[[i+1]] <- p.p 
    }
    
    print(paste(dispersal,prd,min(minValue(stack(pa.eu))),max(maxValue(stack(pa.eu)))))
  }
  

  print(ggarrange(plot.list[[1]]+ rremove("legend"),
                  labels = "a",ncol= 2,heights=c(0.5,1),
                  ggarrange(plot.list[[2]],plot.list[[3]],
                            plot.list[[4]],plot.list[[5]],
                            plot.list[[6]],plot.list[[7]],
                            col = 3,common.legend = TRUE,
                            labels = c("b","c","d","e","f","g"))))

  
}


#_______________________________________________________________________________________

#    --------------------------- 9. Regional Range Overlap  -----------------------

#_______________________________________________________________________________________
#calcualte the overlap between apple and species distribution within each country/region

cores=detectCores(logical = FALSE)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
foreach(sp=seq_along(sp.list)) %dopar% {

  library(plyr)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(parallel)
  library(sf)
  library(ggplot2)

  #for (sp in 1:length(sp.list)) {
  spp <- sp.list[sp]

  kap <- binary.rasters2[grepl(spp, binary.rasters2)]

  table = data.frame()#table to fill

  for (s in 1:length(scenario.list.pres)) {
    scenario <- scenario.list.pres[s]

    kap2 <-  kap[grepl(scenario, kap)]

    for (d in 1:1) {
      dispersal <- dispersal.list[d]

      kap3 <-  kap2[grepl(dispersal, kap2)]

      for (run in 1:length(kap3)) {

        kap4 <-  kap3[run]

        ras.bee <- raster(kap4)

        runnum <- strsplit(kap4,"_")[[1]][7]

        for (a in 1:2) {
          apple.dispersal <- apple.dispersal.list[a]

              if (a == 1) {
              apple.ras.d <-
                apple.rasters.con[!grepl("limited", apple.rasters.con)]
            } else{
              apple.ras.d <-
                apple.rasters.con[grepl("limited", apple.rasters.con)]
            }

            apple.ras.d2 <- apple.ras.d[grepl(scenario,apple.ras.d)]

            ras.apple <- raster(apple.ras.d2)

            ras.bee[ras.bee == 1] = 2 #reclassify to make distinction between apple
            kap.apple <-
              ras.bee + ras.apple #1=apple only, 2=species only, 3=overlap

            AT = data.frame()#table to fill

            for (r in 1:length(poly.apple$ADMIN)){

              country <- poly.apple$ADMIN[r]
              reg <- poly.apple[poly.apple$ADMIN == country, ]#subset by site

              reg_sp <- as(reg, 'Spatial')

              #create outside buffer
              b_large <- gBuffer(reg_sp, width = 50000, quadsegs = 100)

              #crop with large buffer
              ras2 <- crop(kap.apple, extent(b_large))#crop to large buffer
              ras3 <-
                raster::mask(ras2, b_large)#extract raster from within large buffer

              #turn polygon into raster, so we do not lose edges
              SpP_ras <-
                rasterize(reg, ras3, getCover = TRUE)#rasterisize the buffer
              SpP_ras[SpP_ras == 0] <- NA #remove everything outside of buffer

              #get data
              ras4 = stack(raster::mask(ras3, SpP_ras))#extract data from within buffer

              #simplify for clipping
              reg.clip <- rgeos::gSimplify(reg_sp, tol = 10)

              #convert raster to polgyon to get exact area
              ras_vec <-
                rasterToPolygons(
                  ras4,
                  fun = NULL,
                  n = 16,
                  na.rm = TRUE,
                  digits = 12,
                  dissolve = FALSE
                )

              clip <-
                raster::intersect(ras_vec, reg.clip)#get only vector within the buffer

              #plot(clip)

              #get percentage cover by area
              clip@data$Area_m2 <-
                raster::area(clip) #make new column of area of each polygon
              totalarea = raster::area(reg.clip) #total area of buffer

              #calculate proportion of area
              clip_df <- data.frame(clip)
              areatab <-
                clip_df %>% group_by(layer) %>% summarise(Area_m2 = sum(Area_m2)) #sum by land use type
              areatab$prop = areatab$Area_m2 / totalarea

              at <- data.frame(areatab)
              at$Period <- c(scenario)
              at$Binomial <- spp
              at$Country <- c(country)
              at$Dispersal <- c(dispersal)
              at$Run <- c(gsub(".tif", "", runnum))
              at$AppleDispersal = c(apple.dispersal)

              AT = rbind(AT, at)#join
            }


            write.csv(AT,paste0("./Overlap/Regional/regional_overlap_",spp,"_",
                                scenario,dispersal,apple.dispersal,"_",gsub(".tif","",runnum),"_",".csv"),
                      row.names=F)
    }
  }
  print(spp)
    }
  }
}
stopCluster(cl)


#get all csv and join them
overlap.region.csv <- list.files(
  "/Overlap/Regional",
  pattern = ".csv$",
  full.names = T
)

orcsv <- unique(rbindlist(lapply(overlap.region.csv, read.csv), idcol = "origin"))

orcsv2 <- orcsv %>% filter(layer==3) %>% group_by(Period,AppleDispersal,Dispersal,Country) %>% 
  summarise_at(names(orcsv)[4], list(mean = mean, sd = sd)) %>% 
  filter(AppleDispersal!="nod"&AppleDispersal!="int"&Dispersal!="nod"&AppleDispersal!="limited"&Dispersal!="int") %>% 
  filter(Period=="present"|Period=="rcp26.61.80"|Period=="rcp45.61.80"|Period=="rcp85.61.80") %>% 
  ungroup() %>% select(1,4,5) %>% pivot_wider(names_from = Period, values_from = mean)

orcsv2$OverlapLoss26 <- 1-orcsv2$rcp26.61.80/orcsv2$present
orcsv2$OverlapLoss26 <- orcsv2$OverlapLoss26 %>% replace_na(0)
orcsv2$OverlapLoss45 <- 1-orcsv2$rcp45.61.80/orcsv2$present
orcsv2$OverlapLoss45 <- orcsv2$OverlapLoss45 %>% replace_na(0)
orcsv2$OverlapLoss85 <- 1-orcsv2$rcp85.61.80/orcsv2$present
orcsv2$OverlapLoss85 <- orcsv2$OverlapLoss85 %>% replace_na(0)

#missing countries
setdiff(levels(factor(orcsv$Country)),levels(factor(orcsv2$Country)))

#orcsv2 <- rbind(orcsv2,c("Albania",replicate(7, 0)),c("Ireland",replicate(7, 0)))

#remove Macedonia as an outlier
orcsv2.1 <- orcsv2 %>% filter(Country!="Macedonia")

orcsv.poly <- merge(poly.c.change, orcsv2.1, by.x="ADMIN",by.y="Country")

plot(orcsv.poly["OverlapLoss26"])
orcsv.poly$OverlapLoss26 <- as.numeric(orcsv.poly$OverlapLoss26)
orcsv.poly$OverlapLoss45 <- as.numeric(orcsv.poly$OverlapLoss45)
orcsv.poly$OverlapLoss85 <- as.numeric(orcsv.poly$OverlapLoss85)

#make bivariate maps of overlap
#Yield/Loss Map
data <- bi_class(orcsv.poly, x = Yld.h.., y = OverlapLoss85, style = "fisher", dim = 3)
map <- ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "DkBlue", dim = 3) +
  labs(
  ) +
  bi_theme()
legend <- bi_legend(pal = "DkBlue",
                    dim = 3,
                    xlab = "Higher Yield (hg/ha)",
                    ylab = "Higher Overlap Loss (%)",
                    size = 8)
finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.2, .65, 0.2, 0.2)

finalPlot

#Production/Loss Map
data <- bi_class(orcsv.poly, x = Prdct.., y = OverlapLoss85, style = "fisher", dim = 3)
map <- ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "DkCyan", dim = 3) +
  labs(
  ) +
  bi_theme()
legend <- bi_legend(pal = "DkCyan",
                    dim = 3,
                    xlab = "Higher Production (t)",
                    ylab = "Higher Overlap Loss (%)",
                    size = 8)
finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.2, .65, 0.2, 0.2)

finalPlot
