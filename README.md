## Potential for climate change driven spatial mismatches between apple crops and their wild bee pollinators at a continental scale

### Journal: Global Environmental Change

Article DOI: https://doi.org/10.1016/j.gloenvcha.2023.102742

##### Author list: Leon Marshall<sup>1,2</sup>\*, Nicolas Leclercq<sup>1</sup>, Timothy Weekers<sup>1</sup>, Insafe El Abdouni<sup>3</sup>, Luísa G. Carvalheiro<sup>4,5</sup>, Michael Kuhlmann<sup>6</sup>, Denis Michez<sup>3</sup>, Pierre Rasmont<sup>3</sup>, Stuart P.M. Roberts, Guy Smagghe<sup>7</sup>, Peter Vandamme<sup>8</sup>, Thomas Wood<sup>3</sup>, Nicolas J. Vereecken<sup>1</sup>,

1.	Agroecology Lab, Université libre de Bruxelles École Interfacultaire de Bioingénieurs, Avenue F.D. Roosevelt 50, cp 264/2, B-1050 Bruxelles, Belgium.
2.	Naturalis Biodiversity Center, Darwinweg 2, 2333 CR, Leiden, Netherlands.
3.	Laboratory of Zoology, Université de Mons, Bd Dolez 31, 7000 Mons, Belgium.
4.	Depto de Ecologia, Univ. Federal de Goias, Av. Esperança, s/n - Chácaras de Recreio Samambaia, Goiânia - GO, 74690-900, Brazil.
5.	Centre for Ecology, Evolution and Environmental Changes (CE3C), Faculdade de Ciências Universidade de Lisboa, Edifício C2, 5º Piso, Sala 2.5.46, 1749-016 Lisboa, Portugal.
6.	Zoological Museum, University of Kiel, Hegewischstraße 3, D–24105 Kiel, Germany.
7.	Department of plants and crops, Ghent University, Coupure links 653, B-9000, Ghent, Belgium.
8.	Department of biochemistry and microbiology, Ghent University, Technologiepark 927, B-9052, Ghent, Belgium.
\* Corresponding author

Keywords: Pollination, Species distribution modeling, MaxEnt, Food systems, Land cover, Agriculture

### Abstract
Visitation by wild bee species alongside managed pollinators is necessary to ensure consistent yields and fruit quality in apple fields. Wild bee species are vulnerable to several environmental changes. Climate change is expected to lead to broad-scale changes to wild bee distributions that will impact the service they provide as crop pollinators. We modelled selected wild bee species known to be important for apple production in Europe and we quantified the shifts in distribution range for these key apple-pollinating bee species (KABS) under three climate‌ change scenarios (RCP 2.6, 4.5 and 8.5) for 2041–2060 and 2061–2080. We compared species distribution maps (after the expected range changes) to the distribution of areas with suitable habitat for apple orchards and with national apple production statistics to estimate potential pollination service at the landscape scale. Overall, ‌KABS are widespread species found across Europe and while most species have projected range contractions, these contractions are limited (∼10% loss). Only under the worst-case climate change scenario (RCP8.5) do we project range contractions over 50% and only under RCP8.5 is the average loss of overlap between suitable apple conditions and KABS likely to decrease by over 10%. However, range contractions at the southern limit of many species’ ranges mean that the potential impact of climate change on apple pollination is not evenly shared between apple producing countries; France and Italy for example are projected to have high range loss of KABS and loss in potential pollination service. Climate change is not the only threat to apple pollination and future pollination deficits will also depend on local orchard intensification and ecological infrastructure. Key changes to intensive, commercial apple orchards towards a more agroecological approach are needed to maintain a diverse wild bee community and apple production in areas that may become climatically unsuitable in the future.

![IMG_20220127_105650](https://user-images.githubusercontent.com/33490288/188571470-752677ee-0e22-41e0-875b-d815d9c1849d.jpg)
Parque Natural Aguas de Ramón, Santiago, Chile. Photo by Leon Marshall (CC BY 2.0)

### Data and Script Files
 - ModelScript.R - r script file to run all models
 - rasterstack.RData - includes all rasters used in the modeling process - loaded as part of ModelScript.R
         - format: list of 8 elements 
         
                      1. Full training scale of present period
                      
                      2. Present period at projected scale
                      
                      3. Rcp2.6 2041–2060 
                      
                      4. Rcp4.5 2041–2060 
                      
                      5. Rcp8.5 2041–2060 
                      
                      6. Rcp2.6 2061–2080
                      
                      7. Rcp4.5 2061–2080
                      
                      8. Rcp8.5 2061–2080
                      
                           e.g. Rcp2.6.2041.60 <- stack.ras[[3]]
   
         - see suppl. table TableS4_variabletable.csv for full list of all variables.

- Appledata.csv - Occurerence records of apple orchards in Europe -- Source LUCAS 2009, 2012, 2015, 2018 -- Eurostat, 2021. Land Cover/Use Statistics (LUCAS). Database. EUROSTAT.. http://ec.europa.eu/eurostat/web/lucas/data/database.
- rasterstackapples.RData - raster data for modelling apple occurrences - replaces land use data with soil data.
- post_modeling_analysis.R - script for analysis of model output and making figures.
- Maps folder - all maps needed to duplicate figures in publication - see metadata in folder for more information.
- Bee species collection data needed to run the models can be downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.5tb2rbp8n



 
