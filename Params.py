#####   MODEL PARAMETERS FOR POPULATION DYNAMICS ####
beta0 = 1  #Infectious contact rate
bi = 3.5  #Per capita birth rate
di = 0.5 #Per capita natural death rate
k = 400  #Carrying capacity
omega = (bi -di) / k # Strength of density dependance on natality
gamma = 3  #Parasite Clearance
alpha = 1.5  #Parasite Virulence
rho = 0.9 #Dispersal Cost
epsilon = 0.1 #Extinction rate