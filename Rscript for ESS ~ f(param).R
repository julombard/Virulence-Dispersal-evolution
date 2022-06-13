#Try for analysis of multiple output files (densities & traits) for a particular range of parameters
#in order to produce figure of evolutionnary trajectories as function of parameter
library(readr)
library(dplyr)
library(plot3D)
library(plotly)

### FUNCTIONS USED TO BUILD MEAN SIMULATION ###
MeanSimValue <- function(lsdf, var){ # return mean value for each time of a given variable
  
  nbsims <- length(lsdf)
  values <- vector('list', nbsims)
  
  for (i in 1: nbsims){
    idcol <- which(colnames(lsdf[[i]])== var)
    values[[i]] <- lsdf[[i]][,idcol]
  }
  as.data.frame(values)
  meanvalues <- rowSums(as.data.frame(values), na.rm = T)/nbsims
  meanvalues
  
} 



BuildMeanSim <- function(lsdf){ # Build the mean dataframe from the whole set of simulations
  namecols<- names(lsdf[[1]])
  MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(lsdf[[1]]$t)))
  colnames(MeanSim) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(lsdf,namecols[i])
    MeanSim[,i] <- meanvar
  }
}


################# Get Files #######
dir <- choose.dir() # Has to go to global folder in which all sims are stored in subfolders

#### Build a list of folder DENSITIES & TRAITS ####
#Read densities files
cwd <- setwd(dir)
listfiles <- list.files(cwd)
listfiles

DensitiesFiles <- list()
TraitsFiles <- list()
DistribFiles <- list()

for (i in 1:length(listfiles)){
  
  currentdir <- paste(dir, listfiles[i], sep = "\\") #Browse subfolders in main folders
  TraitsDir <- paste(currentdir, "Traits", sep = "\\") # Get the subsubfolder with traits files
  DensitiesDir <- paste(currentdir, "Densities", sep = "\\") # Get the subsubfolder with densities files
  DistribDir <- paste(currentdir, "Distributions", sep = "\\") # Get the subfolder for distribution of traits when discrete trait values are used
  
  cwd <- setwd(TraitsDir)
  
  ls_files_Traits <- list.files(TraitsDir) # list files
  list_tibbles_traits <- lapply(ls_files_Traits, read_csv) # Read CSV returns tibble object but i don't know how to deal with thoses guys so...
  list_df_traits <- lapply(list_tibbles_traits, as.data.frame) #... we turn them to good old dataframes
  
  ### This part contruct the mean simulation from the Traits list of file we got earlier ###
  #Can be improved using Xapply mastery
  
  # I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_traits)){
    
    list_df_traits[[i]]<- subset(list_df_traits[[i]][c(-1)]) # suppress first column
    list_df_traits[[i]]$t <- round(list_df_traits[[i]]$Time, 2) #round t values
    
    #Avoid duplicates
    doubles <- which(duplicated(list_df_traits[[i]]$t))
    list_df_traits[[i]] <- list_df_traits[[i]][-doubles,]
  }
  #Get only commons rows between all dataframes 
  #Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_traits)){
    times <- list_df_traits[[i]]$t
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  counts
  #Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_traits)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  #loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_traits)){
    indexrows <- which(list_df_traits[[i]]$t %in% tcommon)
    list_df_traits[[i]] <- subset(list_df_traits[[i]][indexrows,])
  }
  
  # II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_traits[[1]])
  MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_traits[[1]]$t)))
  colnames(MeanSim) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_traits,namecols[i])
    MeanSim[,i] <- meanvar
  }
  

  #And add to the list
  TraitsFiles <- c(TraitsFiles, list(MeanSim))
  
  
  ### Now we repeat the same process for Densities
  
  cwd <- setwd(DensitiesDir)
  
  ls_files_Densities <- list.files(DensitiesDir) #list files
  list_tibbles_densities <- lapply(ls_files_Densities, read_csv) 
  list_df_densities <- lapply(list_tibbles_densities, as.data.frame) 
  
  
  # I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_densities)){
    
    list_df_densities[[i]]<- subset(list_df_densities[[i]][c(-1)]) # suppress first column
    list_df_densities[[i]]$t <- round(list_df_densities[[i]]$t, 2) #round t values
    
    #Avoid duplicates
    doubles <- which(duplicated(list_df_densities[[i]]$t))
    list_df_densities[[i]] <- list_df_densities[[i]][-doubles,]
  }
  #Get only commons rows between all dataframes 
  #Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_densities)){
    times <- list_df_densities[[i]]$t
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  counts
  #Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_densities)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  #loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_densities)){
    indexrows <- which(list_df_densities[[i]]$t %in% tcommon)
    list_df_densities[[i]] <- subset(list_df_densities[[i]][indexrows,])
  }
  
  # II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_densities[[1]])
  MeanSim2 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_densities[[1]]$t)))
  colnames(MeanSim2) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_densities,namecols[i])
    MeanSim2[,i] <- meanvar
  }
  
  
  DensitiesFiles <- c(DensitiesFiles, list(MeanSim2))
  
  
  # Now we re-repeat the same process for distributions
  ### Now we repeat the same process for Densities
  
  cwd <- setwd(DistribDir)
  
  ls_files_Distribtions <- list.files(DistribDir) #list files
  list_tibbles_Distributions <- lapply(ls_files_Distribtions, read_csv) 
  list_df_Distributions <- lapply(list_tibbles_Distributions, as.data.frame) 
  
  
  # I. Shaping dataframes and getting common rows
  
  for (i in 1:length(list_df_Distributions)){
    
    list_df_Distributions[[i]]<- subset(list_df_Distributions[[i]][c(-1)]) # suppress first column
    list_df_Distributions[[i]]$t <- round(list_df_Distributions[[i]]$Time, 2) #round t values
    
    #Avoid duplicates
    doubles <- which(duplicated(list_df_Distributions[[i]]$t))
    list_df_Distributions[[i]] <- list_df_Distributions[[i]][-doubles,]
  }
  #Get only commons rows between all dataframes 
  #Get all t values
  t_values <- c(0)
  for (i in 1: length(list_df_Distributions)){
    times <- list_df_Distributions[[i]]$t
    t_values <- c(t_values, times)
  }
  t_values <- t_values[-1] # Suppress the value given for initialization
  counts <- as.data.frame(table(t_values)) # Count occurrences into a table
  counts
  #Get values that appear nbsim times (common to every df)
  nbsim <- length(list_df_Distributions)
  tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df
  
  #loop to keep only those t in the whole set of dfs
  for (i in 1:length(list_df_Distributions)){
    indexrows <- which(list_df_Distributions[[i]]$t %in% tcommon)
    list_df_Distributions[[i]] <- subset(list_df_Distributions[[i]][indexrows,])
  }
  
  # II. Building the mean simulation ( returns a single DataFrame)
  namecols<- names(list_df_Distributions[[1]])
  MeanSim3 <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df_Distributions[[1]]$t)))
  colnames(MeanSim3) <- namecols
  
  for (i in 1:length(namecols)){
    
    meanvar <- MeanSimValue(list_df_Distributions,namecols[i])
    MeanSim3[,i] <- meanvar
  }
  
  
  DistribFiles <- c(DistribFiles, list(MeanSim3))
  
  
} # GIANT LOOP THAT EXTRACT, COMPILE AND BUILD A MEAN SIM WITH ALL FILES

### DEMOGRAPHY ANALYSIS ###

####### Average Dynamics of each metapopulation #####

for (i in 1:length(DensitiesFiles)){
  
  ActualDf <- DensitiesFiles[[i]]
  DeuxDf <- DensitiesFiles[[i]] # ?
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  
  Metapop_Data <- data.frame(t = ActualDf$t, S = Smean, I=Imean, N = Nmean)
  
  plot(Metapop_Data$t, Metapop_Data$S, type = 'l', col='blue', main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,40000))
  lines(Metapop_Data$t,Metapop_Data$I, col = 'red')
  lines(Metapop_Data$t, Metapop_Data$N, col = 'black')
  legend("topright",legend = c("S","I", "N"), col = c("blue", "red", "black"), lty = 1 )
  
  
  ###### Phase Space #####
  #Select only equilibrium values
  #plot(Metapop_Data$S, Metapop_Data$I, type = 'l', xlab = 'Susceptibles', ylab = 'Infected')
}


### FIGURE THAT GIVES THE TOTAL POPULATION SIZE AS A FUNCTION OF PARAMETER ####
ActualDf <- DensitiesFiles[[1]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean

Metapop_Data <- data.frame(t = ActualDf$t, S = Smean, I=Imean, N = Nmean)
plot(Metapop_Data$t, Metapop_Data$N, type = 'l', col=1, main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,40000))


for (i in 2:length(DensitiesFiles)){
  
  ActualDf <- DensitiesFiles[[i]]
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  
  Metapop_Data <- data.frame(t = ActualDf$t, S = Smean, I=Imean, N = Nmean)
  lines(Metapop_Data$t, Metapop_Data$N, col = i)


}
legend(150,41000, legend = c('e = 0.1','e =0.2','e =0.3','e =0.4', 'e =0.5'), col = c(1,2,3,4,5), lty = 1:2, cex = 0.8)
### FIGURE THAT GIVE THE SUSCEPTIBLE DENSITY AS A FUNCTION OF PARAMETERS ####

ActualDf <- DensitiesFiles[[1]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)

Metapop_Data <- data.frame(t = ActualDf$t, S = Smean)
plot(Metapop_Data$t, Metapop_Data$S, type = 'l', col=1, main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,40000))


for (i in 2:length(DensitiesFiles)){
  
  ActualDf <- DensitiesFiles[[i]]
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  
  Metapop_Data <- data.frame(t = ActualDf$t, S = Smean)
  lines(Metapop_Data$t, Metapop_Data$S, col = i)
  
  
}
legend(100,41000, legend = c('e = 0.1','e =0.2','e =0.3','e =0.4', 'e =0.5'), col = c(1,2,3,4,5), lty = 1:2, cex = 0.8)

### FIGURE THAT GIVE THE INFECTED FRACTION OF POPULATION AS A FUNCTION OF PARAMETER ####

ActualDf <- DensitiesFiles[[1]]

###### Mean Sim TAU
taille_col = (length(ActualDf)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(ActualDf[generate_even_indexes])
Ipops <- subset(ActualDf[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)
Nmean <- Smean + Imean
Ifrac <- Imean / Nmean

Metapop_Data <- data.frame(t = ActualDf$t, S = Smean, I=Ifrac, N = Nmean)
plot(Metapop_Data$t, Metapop_Data$I, type = 'l', col=1, main = "Metapopulation dynamics", xlab = "t", ylab = "Fraction I", ylim = c(0,1))


for (i in 2:length(DensitiesFiles)){
  
  ActualDf <- DensitiesFiles[[i]]
  
  ###### Mean Sim TAU
  taille_col = (length(ActualDf)-1) # Number of populations columns 
  nbsites <- taille_col/2
  generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
  generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations
  
  Spops <- subset(ActualDf[generate_even_indexes])
  Ipops <- subset(ActualDf[generate_odd_indexes])
  
  # We sum by lines to get global densities
  Smean <- rowSums(Spops)
  Imean <- rowSums(Ipops)
  Nmean <- Smean + Imean
  Ifrac <- Imean / Nmean
  
  Metapop_Data <- data.frame(t = ActualDf$t, I=Ifrac)
  lines(Metapop_Data$t, Metapop_Data$I, col = i)
  
  
}
legend(100,1, legend = c('e = 0.1','e =0.2','e =0.3','e =0.4', 'e =0.5'), col = c(1,2,3,4,5), lty = 1:2, cex = 0.8)



### TRAITS VALUES ANALYSIS ###

#Mean Trait dynamics for average trait value

ConvergeValues <- c()

for (i in 1:length(TraitsFiles)){
  
  ActualDf <- TraitsFiles[[i]]
  DeuxDf <- TraitsFiles[[i]] # ?
  
  Nbsites <- length(ActualDf)-1
  SitesValues <- subset(ActualDf[c(-1)])
  MeanValues <- rowSums(SitesValues, na.rm = T)/Nbsites
  taillemax <- length(MeanValues)
  ConvergeValue <- MeanValues[taillemax]
  ConvergeValues <- c(ConvergeValues, ConvergeValue)
  plot(ActualDf$t, MeanValues, type = 'l', xlab = 'time', ylab = 'Mean Trait Value')
  

  
}


#### TRACE LA FIGURE AVEC TOUTES LES DYNAMIQUES DES TRAITS + LEGENDE ####

Nbsites <- length(TraitsFiles[[1]])-1
SitesValues <- subset(TraitsFiles[[1]][c(-1)])
MeanValues <- rowSums(SitesValues, na.rm = T)/Nbsites

plot(TraitsFiles[[1]]$t, MeanValues, type = 'l',xlab = 'time', ylab = 'Mean Trait Value', ylim = c(0,1.2))

for (i in 1:length(TraitsFiles)-1){
  
  ActualDf <- TraitsFiles[[i+1]]
  
  Nbsites <- length(ActualDf)-1
  SitesValues <- subset(ActualDf[c(-1)])
  MeanValues <- rowSums(SitesValues, na.rm = T)/Nbsites
  taillemax <- length(MeanValues)
  ConvergeValue <- MeanValues[taillemax]
  ConvergeValues <- c(ConvergeValues, ConvergeValue)
  lines(ActualDf$t, MeanValues, type = 'l', xlab = 'time', ylab = 'Mean Trait Value', col = i)
  
  
  
}
legend(100,0.8, legend = c('e = 0.1','e =0.2','e =0.3','e =0.4', 'e =0.5'), col = c(1,2,3,4,5), lty = 1:2, cex = 0.8)

#### CALCUL OF CONVERGENCE VALUES BASED ON THE MEAN ON 30 LAST 


ConvergeValues <- c()

for (i in 1:length(TraitsFiles)){
  
  ActualDf <- TraitsFiles[[i]]
  
  Nbsites <- length(ActualDf)-1
  
  LastValues<- subset(ActualDf, ActualDf$t > 120)
  SitesValues <- subset(LastValues[c(-1)])
  
  RowMeanValues <- rowMeans(SitesValues, na.rm = T) # Fait la moyenne des lignes
  RowMeanValues[1]
  MeanValues <- mean(RowMeanValues, na.rm = T)# Fait la moyenne des colonnes 
  MeanValues

  ConvergeValues <- c(ConvergeValues, MeanValues)
  
  
  
  
}
ConvergeValues
dvalues = c(0.1,0.3,0.5, 0.7, 0.8, 0.9, 0.95, 1)
plot(dvalues, ConvergeValues, xlab = 'e', ylab='Stable Value', col = 'red', pch =19)


###### ANALYSIS OF THE TRAIT DISTRIBUTIONS FOR DISCRETE ALPHA VALUES ###########


for ( i in 1:length(DistribFiles)){
  
  ActualDf <- DistribFiles[[1]]

  time <- as.vector(ActualDf$t)
  alpha <- seq(0.01,3,0.1)
  Densities <- as.matrix(subset(ActualDf[,-1]))
  Densities
  colnames(Densities) <- NULL

  fig <- plot_ly(z = ~Densities, type = "heatmap")%>%layout(title = 'Phenotypic Distribution dynamics, d=0.95', plot_bgcolor = "#e5ecf6", xaxis = list(title = 'Alpha'), 
                                         
                                         yaxis = list(title = 'Time'),zaxis = list(title = "Density"), legend = list(title=list(text='<b> What is that ? </b>')))

  plotly_build(fig)

}
list_df_densities[[1]]
View(list_df_densities[[2]])
