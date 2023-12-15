########################################
# script create_folders.R
#
# creates folders for storing simulation outputs
#
########################################

create_folders <- function(exp){
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  GROUP = "/scicore/home/penny/malinga/malaria_hotspots/simulations_nov/"
  
  dir.create(paste0(GROUP,exp))
  dir.create(paste0(GROUP,exp,"/codefiles"))
  
  file.copy(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov/codefiles"), 
            paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov/",exp),recursive = TRUE)
  
}

create_folders(exp)