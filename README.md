# Climate Change Driven Spatial Mismatches - Bees and Apple Crops
Data and scripts for publication: 

Article DOI: [...](...)

Article: Potential for climate change driven spatial mismatches between apple crops and their wild bee pollinators at a continental scale

Journal: Global Environmental Change

Author list: Leon Marshall\*,
Nicolas Leclercq,
Timothy Weekers,
Insafe El Abdouni,
Luísa G. Carvalheiro,
Michael Kuhlmann,
Denis Michez,
Pierre Rasmont,
Stuart P.M. Roberts,
Guy Smagghe,
Peter Vandamme,
Thomas Wood,
Nicolas J. Vereecken

\* Corresponding author

### Files
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

N.B: All species data can be downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.5tb2rbp8n



 
