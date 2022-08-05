# Stochastic Simulation Algorithm for Metapopulation models
from copy import deepcopy
import classes_demography
import fonctions_demography
import numpy as np
import pandas as pd


#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
#Including New features for metapopulations modelling, designed by Massol F., Lion S. and bibi
#Damn efficient compared to previous try (on lab laptop 100 sites simulated over 16K iterations took 2'30)

#Simulation parameters
seed = 4

np.random.seed(seed) #Set seed for reproducibility
sim_time = 0 # Simulation time (model time, not an iteration number)
vectime = [0] # to keep t variable
tmax = 30 # Ending time
Nexactsteps = 20  # Number of steps to do if/when performing direct method
nbsite = 40 # Number de sites
Taillepop = 100 # Initial local population sizes

#Model parameters
beta0 = 1  #Infectious contact rate
bi = 2  #Per capita birth rate
di = 0.5 #Per capita natural death rate
omega = 1 - (di /bi) # Strength of density dependance on natality
k = 100  #Carrying capacity
d = 0.5  #Dispersal propensity
gamma = 2.5  #Parasite Clearance
alpha = 0.2  #Parasite Virulence
rho = 0.9 #Dispersal Cost
epsilon = 0.1 #Extinction rate

#Define population as class instances
ListSites = fonctions_demography.SetMetapop(nbsite, Taillepop)

#Event definition
#Further expansion idea : build events from a unique model.txt file read by the program, in order to simulate whathever you want
ReproductionS = classes_demography.Event(name='Reproduction S',propensity='bi*(1- omega * (self.S+self.I)/k) * self.S', Schange='1', Ichange='0', order=1)
DeathS = classes_demography.Event(name='Death S',propensity='di*self.S', Schange='-1', Ichange='0', order=3)
DispersalS = classes_demography.Event(name='Dispersal S',propensity='d*self.S', Schange='-1', Ichange='0', order=1)
DispersalI = classes_demography.Event(name='Dispersal I',propensity='d*self.I', Schange='0', Ichange='-1', order=1)
Extinction = classes_demography.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0)
Infection = classes_demography.Event(name='Infection',propensity='(beta0 *alpha / (1+alpha))*self.S*self.I ', Schange='-1', Ichange='1', order=2)
Recovery = classes_demography.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1)
DeathI = classes_demography.Event(name='Death I',propensity='(di + alpha) *self.I', Schange='0', Ichange='-1', order=1)

#Event vector, cause tidying up things is always nice
Events = [ReproductionS,Infection, DispersalI, DispersalS, DeathS, DeathI,Recovery, Extinction]

#Initializing outputs storage
Densities_out = [] # Collect densities outputs
Propensities_out =[] # Collect propensities outputs
IsTrackPropensites = False #Set false/ true to not track/track propensities
for i in range(len(Events)):
    Propensities_out.append([0]) #Store times series of propensities

#Get initial lists and values in outputs
#We want to get one list per Sx(t) and Ix(t) to store them easily in a dataframe at the end
for i in ListSites:
    Densities_out.append([i.effectifS])
    Densities_out.append([i.effectifI])

#Get ready for the simulation
#Further expansion idea : import multiprocessing to launch as many runs as CPUnits in the same time

#Model main Loop
while sim_time < tmax :
    print('We have currently passed', sim_time,'time in the simulation') # Kind of a loading bar
    vectime.append(sim_time) #Update time vector

    #Compute the propensities
    Propensities, Sum_propensities = fonctions_demography.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites
    print("Propensities", Propensities)
    SumS, SumI = fonctions_demography.SumDensities(ListSites) # Get total densities of species

    #print('Les props', Propensities)
    #print('Les sommes',Sum_propensities)

    #Get Critical Reactions (maybe not useful since we use multinomial samples to identifiy sites where reactions occurs)
    #Criticals = fonctions.GetCriticals(Propensities, ListSites, Events)

    #We now can compute vectors mu and sigma using previous shit
    Mus = fonctions_demography.ComputeMuNSigma(Sum_propensities, Events) # As each statechange is 1 , 0, or -1 we have sigma = mu
    #print('les mumu', Mus)

    #Get epsilon_i
    Epsis = fonctions_demography.GetEpsilonI(SumS, SumI)
    #print('Les ei*xi', Epsis)
    #Get Tau prime
    TauPrime = fonctions_demography.GetTauPrime(Epsis, Mus)

    # Now that main intermediary computations are done, let's get to the main algorithm Decision tree
    aox = sum(Sum_propensities)

    if TauPrime < 10/aox : # Take 10/aox 1 is left for ignoring this part
        print('Direct Method performed')
        Tau = fonctions_demography.DoDirectMethodV2(Propensities,Sum_propensities,Nexactsteps, Events,ListSites)
    else:
        print('Lets leap')
        #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
        #We expect that random sample of the place of reactions will be equivalent
        #As critical reactions in critical population will have low ocurrences
        Tau=TauPrime # Divided by two for a error measurement test

        #print('Voici Tau', Tau)

        #So we directly sample the kjs (number of realisation of a given event during tau) from a poisson distribution
        Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
        #print('Poisson', Poisson_means) # it's working so we are glad


        ### UPDATE 9.06.2022 : Here we sample FOR EACH site in a poisson distribuion #####
        Propensities_per_site, Sum_propensities_per_site = fonctions_demography.GetPropensites_Per_sites(ListSites, Events)


        #Now we sample the kjs in poisson law, aka the number of trigger of each event
        triggers = []
        for i in Poisson_means :
            kj = np.random.poisson(i,1)
            triggers.append(kj[0]) # The [0] is due to array structure of kj
        #print('Occurrences', triggers)

        #Creating a vector with sites indexes, used later and put here to not be computed at each subloop iteration
        sites_indexes = []
        for i in range(len(ListSites)):
            sites_indexes.append(i)

        #Now we sample the sites where events will occur from multinomial law
        #And apply the effect of event
        for index,event in enumerate(Events) : # For each event
            #This part define the number of occurences per site
            Noccur = triggers[index] #We get the number of times it should trigger during tau
            props = Propensities[index] # We get the propensity per sites

            SumProp = sum(props) # We get the total propensity
            Probas = [i /SumProp for i in props] # We get probability of occurence in each site
            #print('les probas', Probas) #Good job boy


            if Noccur == 0 : #Case where the event can't happen
                trigger_persite=[0 for i in range(nbsite)]
            else : # Normal cases
                trigger_persite = np.random.multinomial(Noccur, Probas)
            #print('Occurrences per sites', trigger_persite)
            #This part apply the effect of events in site populations
            for index, Site in enumerate(ListSites) :
                if 'Dispersal' in event.name :
                    # Multiply the state change in population by the number of triggers
                    Site.effectifS += trigger_persite[index] * event.Schange
                    Site.effectifI += trigger_persite[index] * event.Ichange
                    nbmigrants = max(abs(trigger_persite[index] * event.Schange), abs(trigger_persite[index] * event.Ichange))
                    #print('Nombre de migrants', nbmigrants)

                    #Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
                    SuccessfulMigrants = 0
                    for i in range(nbmigrants):
                        roll4urlife = np.random.uniform(0,1,1)
                        if roll4urlife > rho : SuccessfulMigrants += 1
                    #Here we distribute successful migrants among neighboring sites
                    #This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
                    #Implemented this way, i should be easy to simulate on any network type
                    for i in range(SuccessfulMigrants) :
                        index_sites = deepcopy(sites_indexes) # working copy of site indexes vector
                        del index_sites[index]#Drop the current site from the list cause you can't emigrate to the place from which you departed

                        Index_destination = np.random.choice(index_sites) #Get index of destination site

                        #Add individual to destination
                        if abs(event.Schange) > 0 : #if S are dispersers
                            ListSites[Index_destination].effectifS += 1
                        elif abs(event.Ichange) > 0 : # if I are dispersers
                            #print('PYCHARM WAS HERE')
                            ListSites[Index_destination].effectifI += 1
                        else : print('ERROR : disperser is neither S nor I and that is very curious !') #This is useless, the error never raises
                else:
                    #Multiply the state change in population by the number of triggers
                    if event.name == 'Extinction' :
                        #When extinction occur we need to retrieve densities values cause they're initialized to zero otherwise in class definition
                        event.Schange=-Site.effectifS #I dont really understand why attribute is not protected but it's good news
                        event.Ichange=-Site.effectifI
                        #print('Event Schange', trigger_persite[index])
                        Site.effectifS += trigger_persite[index]*event.Schange
                        #print('Effectif après extinction', Site.effectifS)
                        Site.effectifI += trigger_persite[index]*event.Ichange
                    else :
                        Site.effectifS += trigger_persite[index] * event.Schange
                        Site.effectifI += trigger_persite[index] * event.Ichange

    #Update time
    sim_time += Tau
    print('Effectifs', SumS, SumI)
    #print('time increment', Tau)
    #print('Le temps passe si vite',sim_time)

    #Update the output tracking
    #1. Densities
    indexlist = 0
    for i in ListSites :
        if i.effectifS < 0 : #Avoid negative population in the "big fat brute" way
            i.effectifS = 0
        Densities_out[indexlist].append(i.effectifS)
        indexlist += 1
        if i.effectifI<0:
            i.effectifI = 0
        Densities_out[indexlist].append(i.effectifI)
        indexlist += 1
    #2. Propensities
    if IsTrackPropensites == True :
        Sum_propensities # Propensities of each event in a list sorted by event
        for index,propensitiy in enumerate(Sum_propensities) :
            Propensities_out[index].append(propensitiy)


    #Not done yet but define a boolean IsTrackPropensity in sim options and track ioi = 1

    # Break the main loop if there are no infected remaining ( this happens essentially at start if the 1st infected dies)
    # if SumI == 0:
    #     print('WARNING : ABORTED SIMULATION, No infected remaining')
    #     break

#Structuring outputs to get a .csv file event if the loop has broken
if IsTrackPropensites == True : #if we track propensities
    #Creating propensities dataframe
    dataprop = pd.DataFrame(columns=['t'])
    for event in Events :
        EventName = event.name
        dataprop[EventName]= []
    #Filling the dataframe
    for index, colname in enumerate(dataprop):
        if index == 0 : dataprop[colname] = vectime
        else : dataprop[colname] = Propensities_out[index-1]
    #Saving into .csv file
    dataprop.to_csv('Propensities_outputs_'+str(seed)+'.csv')


#Creating the time series dataframe
data = pd.DataFrame(columns=['t'])
for i in range(nbsite): # Define a column for each subpop and site index
    colname_s = 'S'+str(i)
    colname_i = 'I'+str(i)
    data[colname_s] = []
    data[colname_i] = []
#Filling the dataframe
for index,colname in enumerate(data):
    if index == 0 : data[colname] = vectime
    else : data[colname] = Densities_out[index-1] # I don't remember why i put -1, but it doesn't work if it's removed
print(data)
#Saving into .csv file
data.to_csv('Metapop_Outputs_TOdemo'+str(seed)+'.csv')
# Strange, csv is saved in model directory on lab computer but in User/Appdata/local/temp on personal computer...


