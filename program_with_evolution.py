# STOCHASTIC SIMULATION ALGORITHM FOR EVOLUTION IN SIR METAPOPULATION
import os
from copy import deepcopy
import NetworkFunctions
import classes
import fonctions
import numpy as np
import pandas as pd
import multiprocessing
import Params

# A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
# Adapted for metapopulations modelling and networks

#POPULATION DYNAMICS PARAMETERS
beta0 = Params.beta0 # Infectious contact rate
bi = Params.bi  # Per capita net growth rate
di = Params.di # Per capita natural death rate
omega = Params.omega # Strength of density dependance on births
k = Params.k  # Carrying capacity
#d = 0.5# Dispersal propensity
gamma = Params.gamma  # Parasite Clearance
alpha = Params.alpha  # Parasite Virulence
rho = Params.rho  # Dispersal Cost
epsilon = Params.epsilon  # Extinction rate

def RunModel(seed, param) :
    #SIMULATION PARAMETERS
    Sp_config = "Island"
    # This part changes parameters (here d) values for multisim runs where numerous are tested
    d = param
    classes.d = param
    fonctions.d = param

    print('Seed', seed)
    np.random.seed(seed) #Set seed for reproducibility
    nb_iterations = 0 # Store the number of interations, used to define times that are saved (later)
    sim_time = 0 # Simulation time (model time, not an iteration number)
    vectime = [0] # to keep track of t variable
    tmax = 8000 # Ending time
    Nexactsteps = 20  # Number of steps to do if/when performing direct method (USELESS IF nbsite > 20~30)
    nbsite = 80 # Number of sites
    n = 7 #Number of rows  for lattices configurations
    p= 7 # Number of columns for lattices configurations
    Taillepop = Params.k # Initial local population sizes
    Evoltrait = classes.EvolvingTrait('alpha', True) # Define the trait that will evolve in the simulation (only alpha availables for now)

    #ENABLE OUTPUT STORAGE
    dico_densities_df = {} # Track population densities
    dico_distrib_df = {} # Trait traits values in the whole metapop
    dico_traits_df = {} # Track mean trait value per site

    #CREATE THE METAPOPULATION
    ListSites = fonctions.SetMetapop(nbsite, Taillepop)

    #SPATIAL CONFIGURATIONS OPTIONS AND ENABLING
    if Sp_config == "Square lattice" :
        Hastable_adjacency = NetworkFunctions.Build_Square_Lattice(n,p,ListSites)
        nb_neighbors = 4
    if Sp_config == "Accordion":
        nb_neighbors = 4
        Hastable_adjacency = NetworkFunctions.Build_Accordion_Neighboring(ListSites, nb_neighbors)
    if Sp_config == "Island" : pass # The programm was initially designed for island (complete graph) pop so we just go with the flow
    if Sp_config == "Hexagonal lattice" :
        Hastable_adjacency = NetworkFunctions.Build_Hexagonal_Lattice(n,p,ListSites) # n the number or row desired (int), p the number of columns (int) , ListSite a list of site objects

     #Gives adjacency lists
    #print(Hastable_adjacency)

    #EVENT DEFINITIONS
    ReproductionS = classes.Event(name='Reproduction S',propensity='(bi - omega * (self.S+self.I) ) * self.S', Schange='1', Ichange='0', order=1,EvolvingTrait=Evoltrait)
    DeathS = classes.Event(name='Death S',propensity='di*self.S', Schange='-1', Ichange='0', order=3,EvolvingTrait=Evoltrait)
    DispersalS = classes.Event(name='Dispersal S',propensity='d*self.S', Schange='-1', Ichange='0', order=1,EvolvingTrait=Evoltrait)
    DispersalI = classes.Event(name='Dispersal I',propensity='d*self.I', Schange='0', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    Extinction = classes.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0,EvolvingTrait=Evoltrait)
    Infection = classes.Event(name='Infection',propensity='beta0 *(alpha / (1+alpha) )*self.S*self.I', Schange='-1', Ichange='1', order=2,EvolvingTrait=Evoltrait)
    Recovery = classes.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    DeathI = classes.Event(name='Death I',propensity='(di + alpha) *self.I', Schange='0', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    Events = [ReproductionS,Infection, DispersalI, DispersalS, DeathS, DeathI,Recovery, Extinction]

    #Initializing outputs storage
    Propensities_out =[] # Collect propensities outputs
    IsTrackPropensites = False #Set false/ true to not track/track propensities
    for i in range(len(Events)):
        Propensities_out.append([0]) #Store times series of propensities

    # INITIALIZE OUTPUTS FOT T=0
    # TRACK OF DENSITIES (one list of S and I per site) and MEAN TRAIT VALUES PER SITE
    for index, i in enumerate(ListSites):
        dico_densities_df[f"S{index}"]= [i.effectifS]
        dico_densities_df[f"I{index}"]= [i.effectifI]
        if i.effectifI > 0:
            dico_traits_df[f"site{index}"] = [sum(i.traitvalues) / i.effectifI]
        else:
            dico_traits_df[f"site{index}"] = ["NA"]
    #TRACK TRAIT DISTRIBUTION (more exaclty the densitiy of each trait value trough time)
    for i in fonctions.Rounded_alpha_values:  # For each value defined
        count = 0
        for j in ListSites:  # We browse the different sites
            # We count the given value and sum it
            count += j.traitvalues.count(i)
        dico_distrib_df[f"Alpha{i}"] = [count]


    ############################# MODEL MAIN LOOP ########################################
    while sim_time < tmax :
        #DEFINE HOW MUCH YOU LIKE TO SAVE DATAS
        if nb_iterations % 20 == 0: # % x => we save each x times
            vectime.append(sim_time) # Update time vector

        #COMPUTE PROPENSITIES
        Propensities, Sum_propensities = fonctions.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites
        SumS, SumI = fonctions.SumDensities(ListSites) # Get total densities of species

        # Break the main loop if there are no infected remaining ( This happens essentially at start if the 1st infected dies)
        if SumI == 0:
            print('WARNING : ABORTED SIMULATION, No infected remaining')
            break
        if sim_time % 100.0 == 0.0 :
            print(sim_time, 'time')  # Kind of a loading bar but much uglier

        ################################# TAU-LEAP PART #################################

        # Most of this section is about getting Tau or die trying
        # All explicit computations are available in cao et. al. (2006)

        #Get Critical Reactions (NOT USE : NEGATIVE POPULATION ARE AVOIDED IN THE DIRTY WAY)
        #Criticals = fonctions.GetCriticals(Propensities, ListSites, Events)

        # Preliminaries
        Mus = fonctions.ComputeMuNSigma(Sum_propensities, Events) # As each statechange is 1 , 0, or -1 we have sigma = mu
        Epsis = fonctions.GetEpsilonI(SumS, SumI)
        TauPrime = fonctions.GetTauPrime(Epsis, Mus)

        #Main algorithm Decision tree
        aox = sum(Sum_propensities)

        if TauPrime < 10/aox : # Take 10/aox 1 is left for ignoring this part
            print('Direct Method performed')
            Tau = fonctions.DoDirectMethodV2(Propensities,Sum_propensities,Nexactsteps, Events,ListSites)
        else:
            #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
            #We expect that random sample in the place of reactions will be kind equivalent
            #As critical reactions in critical population will have low ocurrences
            #And if not : we let populations having -1 individuals after an event, and set it back to zero after
            Tau=TauPrime

            ################## DETERMINE HOW MANY EVENTS WILL TRIGGER, AND HOW MANY TIMES PER SITE #####################

            #Sample the mean number of realisation of a given event during tau from a poisson distribution
            Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
            #Sample the kjs in poisson law, aka the number of trigger of each event
            triggers = []
            for i in Poisson_means :
                kj = np.random.poisson(i,1)
                triggers.append(kj[0]) # The [0] is due to array structure of kj
            #Creating a vector with sites indexes, used later and put here to not be computed at each subloop iteration
            sites_indexes = []
            for i in range(len(ListSites)):
                sites_indexes.append(i)
            #Now we sample the sites where events will occur from multinomial law
            for index,event in enumerate(Events) : # For each event
                Noccur = triggers[index] #We get the number of times it should trigger during tau
                props = Propensities[index] # We get the propensity per sites
                SumProp = sum(props) # We get the total propensity
                Probas = [float(i /SumProp) for i in props] # We get probability of occurence in each site
                if Noccur == 0 : #Case where the event can't happen
                    trigger_persite=[0 for i in range(nbsite)]
                else : # Normal cases
                    trigger_persite = np.random.multinomial(Noccur, Probas)


                ########## APPLY CHANGES DUE TO EVENTS IN POPULATIONS ##################

                for index, Site in enumerate(ListSites) :

                    if 'Dispersal' in event.name :
                        # Multiply the state change in population by the number of triggers
                        Site.effectifS += trigger_persite[index] * event.Schange
                        Site.effectifI += trigger_persite[index] * event.Ichange
                        nbmigrants = max(abs(trigger_persite[index] * event.Schange), abs(trigger_persite[index] * event.Ichange))

                        # Here we delete the trait value corresponding to dispersing individual
                        # For virulence evolution, hold only for I indviduals
                        dispersers_traitvalues = []
                        dispersers_beta = []
                        if Site.traitvalues and abs(
                                event.Ichange) > 0:  # If there are dispersers AND that those dispersers are infected
                            for i in range(trigger_persite[index]):  # for each disperser
                                if Site.traitvalues:  # In case we remove all trait values during the loop (induces error messages)
                                    disperser = float(np.random.choice(
                                        Site.traitvalues))  # Get the value that is to be depleted and added to the receiving site
                                    Index_Disperser = Site.traitvalues.index(disperser)  # Get his corresponding beta
                                    dispersers_traitvalues.append(disperser)
                                    dispersers_beta.append(Site.betaI[Index_Disperser])
                                    Site.traitvalues.remove(disperser)  # Remove disperser from actual site
                                    Site.betaI.pop(Index_Disperser)  # Remove corresponding beta
                                else:  break

                        #Aplpy dispersal Cost
                        SuccessfulMigrants = 0
                        for i in range(nbmigrants):
                            roll4urlife = np.random.uniform(0,1,1)
                            if roll4urlife > rho : SuccessfulMigrants += 1

                        #DISTRIBUTE SUCCESSFUL MIGRANTS TO NEIGHBORS
                        for i in range(SuccessfulMigrants):
                            if Site.traitvalues and abs(event.Ichange) > 0:
                                # Get trait values of surviving individuals
                                SurvivorTrait = np.random.choice(
                                    dispersers_traitvalues)  # All values are chosen with same probability
                                Index_survivor = dispersers_traitvalues.index(SurvivorTrait) # Get index
                                SurvivorBeta = dispersers_beta[Index_survivor]  # Get corresponding beta
                                dispersers_traitvalues.remove(SurvivorTrait)  # Chosen value is removed
                                dispersers_beta.remove(SurvivorBeta)  # Remove it from the list

                            #MIGRANT DISTRIBUTION FOR COMPLETE GRAPH "ISLAND"
                            if Sp_config == "Island":
                                index_sites = deepcopy(sites_indexes)  # working copy of site indexes vector
                                del index_sites[
                                    index]  # Drop the current site from the list cause you can't emigrate to the place from which you departed
                                Index_destination = np.random.choice(index_sites)  # Get index of destination site

                            #MIGRANT DISTRIBUTION FOR ANY KIND OF THING THAT RETURNED AND ADJACENCY LIST IN THE BEGINING
                            elif Sp_config != "Island" :
                                Current_site = index # Get the current site index
                                Nearest_neighbors = Hastable_adjacency[f"{Current_site}"] # We get in the topology dictionnary neighbors corresponding to actual site
                                Index_destination = np.random.choice(Nearest_neighbors) # And choose one of them at random
                            #Add individual to destination
                            if abs(event.Schange) > 0 : #if S are dispersers
                                ListSites[Index_destination].effectifS += 1
                            elif abs(event.Ichange) > 0 : # if I are dispersers
                                ListSites[Index_destination].effectifI += 1
                                # Add trait value of individual
                                ListSites[Index_destination].traitvalues.append(SurvivorTrait)
                                ListSites[Index_destination].betaI.append(SurvivorBeta)
                            else : print('ERROR : disperser is neither S nor I and that is very curious !') #This is useless, the error never raises
                    else:
                        if event.name == 'Extinction' :
                            #When extinction occur we need to retrieve densities values cause they're initialized to zero otherwise in class definition
                            event.Schange=-Site.effectifS
                            event.Ichange=-Site.effectifI
                            Site.effectifS += trigger_persite[index]*event.Schange
                            Site.effectifI += trigger_persite[index]*event.Ichange

                            if trigger_persite[index] > 0:  # If the extinction has really occured
                                Site.traitvalues = []
                                Site.betaI = []
                        else :
                            if abs(event.Ichange) > 0:  # If the event has an effect on 'evolving population'
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange
                                Site.traitvalues, Site.betaI = fonctions.ChooseTraitValue(Evoltrait, trigger_persite[index],
                                                                            event.Ichange,Site.traitvalues, Site.betaI)
                            else :
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange

        ############# UPDATE THINGS AT THE END OF ITERATION #############
        #Update time
        sim_time += Tau
        #Update the output tracking
        # 1. Densities
        indexlist = 0
        for index, i in enumerate(ListSites):
            if i.effectifS < 0:  # Avoid negative population in the "big fat brute" way
                i.effectifS = 0
            if nb_iterations % 20 == 0 :
                dico_densities_df[f"S{index}"].append(i.effectifS)
            indexlist += 1
            if i.effectifI < 0:
                i.effectifI = 0
            if nb_iterations % 20 == 0:
                dico_densities_df[f"I{index}"].append(i.effectifI)
            indexlist += 1
        #2. Propensities
        if IsTrackPropensites == True :
             # Propensities of each event in a list sorted by event
            for index,propensitiy in enumerate(Sum_propensities) :
                Propensities_out[index].append(propensitiy)
        # 3. Trait Values
        indexlist2 = 0
        list_value = []
        for index, i in enumerate(ListSites):
            if i.effectifI > 0:
                SumTraitValues = sum(i.traitvalues)
                MeanTraitValue = SumTraitValues / i.effectifI
                if nb_iterations % 20 == 0:
                    dico_traits_df[f"site{index}"].append(MeanTraitValue)
                indexlist2 += 1
            else:
                MeanTraitValue = 'NA'
                list_value.append(MeanTraitValue)
                if nb_iterations % 20 == 0:
                    dico_traits_df[f"site{index}"].append(MeanTraitValue)
        # 4 Count the different phenotypes in the metapopulation, in order to follow their distribution over time
        Possible_values = fonctions.Rounded_alpha_values
        for i in Possible_values:  # For each value defined
            count = 0
            if nb_iterations % 20 == 0:
                for j in ListSites:  # We browse the different sites
                    # We count the given value and sum it
                    count += j.traitvalues.count(i)
                dico_distrib_df[f"Alpha{i}"].append(count)
        #Update the number of iterations (for times where datas are saved)
        nb_iterations += 1

    ################ AFTER THE MAIN LOOP IS ACHIEVED ################
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
        dataprop.to_csv('Propensities_outputs_LongRun_Step005'+str(seed)+'.csv')

    #DENSITIES TIME SERIES
    datadensity = pd.DataFrame.from_dict(data=dico_densities_df)
    VectimeDf = pd.DataFrame(data=vectime)
    datadensity.insert(0, "Time", VectimeDf, allow_duplicates=False)
    datadensity.to_csv('Metapop_outputs_1702_' + str(d) + '_' + str(seed) + '.csv')
    #MEAN TRAITS TIME SERIES
    datatrait = pd.DataFrame.from_dict(data=dico_traits_df)
    datatrait.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datatrait.to_csv('Traits_outputs_1702_' + str(d) + '_' + str(seed) + '.csv')
    #JUST THE TRAITS TIME SERIES
    datadistrib = pd.DataFrame.from_dict(data=dico_distrib_df)
    datadistrib.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datadistrib.to_csv('Distribution_outputs_1702_' + str(d) + '_' + str(seed) + '.csv')

################## MULTIPROCESSING PART ################

# Multiprocessing parameters
list_seeds = [1,2,3,4,5,6] # The list of seed you want to test
list_params =[0.1,0.2,0.3,0.4,0.5] # The list of params values you want to test (has to be changed also at the begining)
#,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1
#list_config =["Island"]
nbsims = len(list_seeds)

#Launch a batch of nbsims simulations
# WARNING : MULTISIM MAKES ERROR MESSAGES VANISH
if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count() #Number of CPU, minus 2 by precaution. And to be able to do things meanwhile
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb) #I don't know what that is doing exactly, but it is necessary.
    for j in range(len(list_params)) :
        for i in range(nbsims):
            #Switch the two following lines to enable/disable multisim
            pool.apply_async(RunModel, args=(list_seeds[i],list_params[j])) #Launch multisim
            #RunModel(list_seeds[i],list_params[j]) #Launch without multisim (restores errors messages)
    pool.close() # Mandatory
    pool.join() # Idem