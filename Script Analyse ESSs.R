# Script for Evolutionnary model of parasite virulence
library(readr)
library(dplyr)



################# Get Files #######
dir <- choose.dir() # Has to go to global folder in which all sims are stored in subfolders

#Here we get folder for different traits values ()
dir01 <- paste(dir, "0.1", sep = "\\")
dir03 <- paste(dir, "0.3", sep = "\\")
dir05 <- paste(dir, "0.5", sep = "\\")
dir08 <- paste(dir, "0.8", sep = "\\")
dir1 <- paste(dir, "1", sep = "\\")


dirTraits <- paste(dir, "Traits Out", sep = "\\")
dirDensities <- paste(dir, "Densities Out", sep = "\\")

ls_files_Traits <- list.files(dirTraits)
ls_files_Densities <- list.files(dirDensities)

#Read densities files
cwd <- setwd(dirDensities)
listfiles <- list.files(cwd)
listfiles
list_tibbles <- lapply(listfiles, read_csv) # Return tibble ? What is that **** ?
list_df <- lapply(list_tibbles, as.data.frame) # Convert to good old dataframe


################################################################################
##################### GET FILES AND TAU LEAP DATAFRAMES ########################
################################################################################


#Supress first column (because useless) and round 't' column to 3 digits (to perform comparison of same t in DM and TauLeap)
for (i in 1:length(list_df)){
  
  list_df[[i]]<- subset(list_df[[i]][c(-1)]) # supress first column
  list_df[[i]]$t <- round(list_df[[i]]$t, 2) #round t values
  #Avoid duplicates
  doubles <- which(duplicated(list_df[[i]]$t))
  list_df[[i]] <- list_df[[i]][-doubles,]
}
#Get only commons rows between all dataframes (ca commence a devenir un poil fasitidieux mébon)
#Get all t values
t_values <- c(0)
for (i in 1: length(list_df)){
  times <- list_df[[i]]$t
  t_values <- c(t_values, times)
}
t_values <- t_values[-1] # Suppress the value given for initialization
counts <- as.data.frame(table(t_values)) # Count occurences into a table
counts
#Get values that appear nbsim times (common to every df)
nbsim <- length(list_df)
tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df

#loop to keep only those t in the whole set of dfs
for (i in 1:length(list_df)){
  indexrows <- which(list_df[[i]]$t %in% tcommon)
  list_df[[i]] <- subset(list_df[[i]][indexrows,])
}

######### FUNCTIONS FOR BUILDING MEAN SIMULATION ######
MeanSimValue <- function(lsdf, var){ # return mean value for each time of a given variable
  
  nbsims <- length(lsdf)
  values <- vector('list', nbsims)
  
  for (i in 1: nbsims){
    idcol <- which(colnames(lsdf[[i]])== var)
    values[[i]] <- lsdf[[i]][,idcol]
  }
  as.data.frame(values)
  meanvalues <- rowSums(as.data.frame(values))/nbsims
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

#In a function -> Doesnt work, oustide -> perfectly working.... voodoo
namecols<- names(list_df[[1]])
MeanSim <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df[[1]]$t)))
colnames(MeanSim) <- namecols

for (i in 1:length(namecols)){
  
  meanvar <- MeanSimValue(list_df,namecols[i])
  MeanSim[,i] <- meanvar
}


 ##########################

######### BUILD MEAN SIMULATION FOR TAU LEAP ##################################

#In a function -> Doesnt work, oustide -> perfectly working.... voodoo
namecols<- names(list_df[[1]])
MeanSimTau <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df[[1]]$t)))
colnames(MeanSimTau) <- namecols

for (i in 1:length(namecols)){
  
  meanvar <- MeanSimValue(list_df,namecols[i])
  MeanSimTau[,i] <- meanvar
}

####### Average Dynamics of one site #####
Tau_raw <- list_df[[1]] 
plot(Tau_raw$t, Tau_raw$S2,type = 'l', col = 'blue',main = "Dynamique d'un site", xlab = 't', ylab = 'Densité')
lines(Tau_raw$t, Tau_raw$I2, col = 'red')
legend("topright",legend = c("S","I"), col = c("blue", "red"), lty = 1 )

####### Average Dynamics of metapopulation #####
###### Mean Sim TAU
taille_col = (length(MeanSimTau)-1) # Number of populations columns 
nbsites <- taille_col/2
generate_even_indexes = seq(2, taille_col, 2) # Even indexes are S populations
generate_odd_indexes = seq(3, taille_col+1, 2) # Odd indexes are I populations

Spops <- subset(MeanSimTau[generate_even_indexes])
Ipops <- subset(MeanSimTau[generate_odd_indexes])

# We sum by lines to get global densities
Smean <- rowSums(Spops)
Imean <- rowSums(Ipops)

Mean_Data_LEAP <- data.frame(t = MeanSimTau$t, S = Smean, I=Imean)

plot(Mean_Data_LEAP$t, Mean_Data_LEAP$S, type = 'l', col='blue', main = "Metapopulation dynamics", xlab = "t", ylab = "Densité", ylim = c(0,40000))
lines(Mean_Data_LEAP$t, Mean_Data_LEAP$I, col = 'red')
legend("topright",legend = c("S","I"), col = c("blue", "red"), lty = 1 )


###### Phase Space #####
#Select only equilibrium values
plot(Mean_Data_LEAP$S, Mean_Data_LEAP$I, type = 'l', xlab = 'Susceptbles', ylab = 'Infected')

#Selecting only quasi-equilibrium values
Mean_EQ_LEAP <- subset(Mean_Data_LEAP, Mean_Data_LEAP$t >8.5)
plot(Mean_EQ_LEAP$S, Mean_EQ_LEAP$I, type = 'l')

#Variance Measures for predictions

var(Mean_EQ_LEAP$S) # 194.43
Sd_LEAPS <- sqrt(var(Mean_EQ_LEAP$S))
var(Mean_EQ_LEAP$I) # 414.08
Sd_LEAPI <- sqrt(var(Mean_EQ_LEAP$I))
cov(Mean_EQ_LEAP$S, Mean_EQ_LEAP$I)#58.8403



##################### TRAITS VALUE ANALYSIS ####################

#Read densities files
cwd <- setwd(dirTraits)
listfiles <- list.files(cwd)
listfiles
list_tibbles <- lapply(listfiles, read_csv) # Return tibble ? What is that **** ?
list_df <- lapply(list_tibbles, as.data.frame) # Convert to good old dataframe

#Supress first column (because useless) and round 't' column to 3 digits (to perform comparison of same t in DM and TauLeap)
for (i in 1:length(list_df)){
  
  list_df[[i]]<- subset(list_df[[i]][c(-1)]) # supress first column
  list_df[[i]]$t <- round(list_df[[i]]$t, 2) #round t values
  #Avoid duplicates
  doubles <- which(duplicated(list_df[[i]]$t))
  list_df[[i]] <- list_df[[i]][-doubles,]
}

#Get only commons rows between all dataframes (ca commence a devenir un poil fasitidieux mébon)
#Get all t values
t_values <- c(0)
for (i in 1: length(list_df)){
  times <- list_df[[i]]$t
  t_values <- c(t_values, times)
}
t_values <- t_values[-1] # Suppress the value given for initialization
counts <- as.data.frame(table(t_values)) # Count occurences into a table
counts
#Get values that appear nbsim times (common to every df)
nbsim <- length(list_df)
tcommon<- counts$t_values[counts$Freq==nbsim] #Return the value of t that are common to every df

#loop to keep only those t in the whole set of dfs
for (i in 1:length(list_df)){
  indexrows <- which(list_df[[i]]$t %in% tcommon)
  list_df[[i]] <- subset(list_df[[i]][indexrows,])
}

#### Build Mean Simulation 
#In a function -> Doesnt work, oustide -> perfectly working.... voodoo
namecols<- names(list_df[[1]])
MeanSimTraits <- data.frame(matrix(ncol = length(namecols), nrow =length(list_df[[1]]$t)))
colnames(MeanSimTraits) <- namecols


MeanSimValue4Traits <- function(lsdf, var){ # return mean value for each time of a given variable
  
  
  nbsims <- length(lsdf)
  values <- vector('list', nbsims)
  
  for (i in 1: nbsims){
    idcol <- which(colnames(lsdf[[i]])== var)
    values[[i]] <- lsdf[[i]][,idcol]
  }
  vadf <- as.data.frame(values)
  vadf
  meanvalues <- rowSums(as.data.frame(values), na.rm = T)/nbsims
  meanvalues
  
} 


for (i in 1:length(namecols)){
  
  meanvar <- MeanSimValue4Traits(list_df,namecols[i])
  MeanSimTraits[,i] <- meanvar
}

# Evolutionnary dynamics for a random site
plot(MeanSimTraits$t, MeanSimTraits$Site0, type = 'l')
plot(MeanSimTraits$t, MeanSimTraits$Site1, type = 'l')

# Mean trait dynamics
Nbsites <- length(MeanSimTraits)-1
SitesValues <- subset(MeanSimTraits[c(-1)])
rowSums(is.na(SitesValues))
MeanValues <- rowSums(SitesValues, na.rm = T)/Nbsites
length(MeanValues)
ConvergeValue <- MeanValues[10000]
ConvergeValue
plot(MeanSimTraits$t, MeanValues, type = 'l', xlab = 'time', ylab = 'Mean Trait Value')
