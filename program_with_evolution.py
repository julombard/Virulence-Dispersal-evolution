# Stochastic Simulation Algorithm for Metapopulation models
from copy import deepcopy
import classes
import fonctions
import numpy as np
import pandas as pd
import multiprocessing
import Params

#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
#Including New features for metapopulations modelling, designed by Massol F., Lion S. and bibi
#Damn efficient compared to previous try (on lab laptop 100 sites simulated over 16K iterations took 2'30)
#Parameters extracted from param files
beta0 = Params.beta0 # Infectious contact rate
bi = Params.bi  # Per capita net growth rate
di = Params.di # Per capita natural death rate
omega = Params.omega # Strength of density dependance on births
k = Params.k  # Carrying capacity
d = 1# Dispersal propensity
gamma = Params.gamma  # Parasite Clearance
alpha = Params.alpha  # Parasite Virulence
rho = Params.rho  # Dispersal Cost
epsilon = Params.epsilon  # Extinction rate

def RunModel(seed, param) :
    #Simulation parameters
    # This part changes parameters value for long autonomous runs
    d = param
    classes.d = param
    fonctions.d = param

    # Stocker les sorties dans un dictionnaire
    dico_densities_df = {}
    dico_distrib_df = {}
    dico_traits_df = {}

    # Simulation parameters
    print('Seed', seed)
    np.random.seed(seed) #Set seed for reproducibility
    sim_time = 0 # Simulation time (model time, not an iteration number)
    vectime = [0] # to keep t variable
    tmax = 0.8 # Ending time
    Nexactsteps = 20  # Number of steps to do if/when performing direct method
    nbsite = 100 # Number de sites
    Taillepop = Params.k # Initial local population sizes

    Evoltrait = classes.EvolvingTrait('alpha', True)
    #Define population as class instances
    ListSites = fonctions.SetMetapop(nbsite, Taillepop)

    #Event definition
    #Further expansion idea : build events from a unique model.txt file read by the program, in order to simulate whathever you want
    ReproductionS = classes.Event(name='Reproduction S',propensity='bi*(1- omega * (self.S+self.I)/k) * self.S', Schange='1', Ichange='0', order=1,EvolvingTrait=Evoltrait)
    DeathS = classes.Event(name='Death S',propensity='di*self.S', Schange='-1', Ichange='0', order=3,EvolvingTrait=Evoltrait)
    DispersalS = classes.Event(name='Dispersal S',propensity='d*self.S', Schange='-1', Ichange='0', order=1,EvolvingTrait=Evoltrait)
    DispersalI = classes.Event(name='Dispersal I',propensity='d*self.I', Schange='0', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    Extinction = classes.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0,EvolvingTrait=Evoltrait)
    Infection = classes.Event(name='Infection',propensity='beta0 *(alpha / (1+alpha) )*self.S*self.I', Schange='-1', Ichange='1', order=2,EvolvingTrait=Evoltrait)
    Recovery = classes.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    DeathI = classes.Event(name='Death I',propensity='alpha*self.I', Schange='0', Ichange='-1', order=1,EvolvingTrait=Evoltrait)

    #Event vector, cause tidying up things is always nice
    Events = [ReproductionS,Infection, DispersalI, DispersalS, DeathS, DeathI,Recovery, Extinction]

    #Initializing outputs storage
    Propensities_out =[] # Collect propensities outputs
    IsTrackPropensites = False #Set false/ true to not track/track propensities
    for i in range(len(Events)):
        Propensities_out.append([0]) #Store times series of propensities

    # Get initial lists and values in outputs
    # We want to get one list per Sx(t) and Ix(t) to store them easily in a dataframe at the end
    list_trait_site = []
    for index, i in enumerate(ListSites):
        dico_densities_df[f"S{index}"]= [i.effectifS]
        dico_densities_df[f"I{index}"]= [i.effectifI]
        if i.effectifI > 0:
            dico_traits_df[f"site{index}"] = [sum(i.traitvalues) / i.effectifI]
        else:
            dico_traits_df[f"site{index}"] = ["NA"]

    # Same for individual values for distributions outputs
    for i in fonctions.Rounded_alpha_values:  # For each value defined
        count = 0
        for j in ListSites:  # We browse the different sites
            # We count the given value and sum it
            count += j.traitvalues.count(i)
        dico_distrib_df[f"Alpha{i}"] = [count]


    ############################# Model main Loop ###########
    while sim_time < tmax :
        vectime.append(sim_time) #Update time vector

        #Compute the propensities
        Propensities, Sum_propensities = fonctions.GetPropensites(ListSites, Events) # Get a vector of propensities ordered by event and by sites
        #print("Propensities", Propensities)
        SumS, SumI = fonctions.SumDensities(ListSites) # Get total densities of species

        # Break the main loop if there are no infected remaining ( this happens essentially at start if the 1st infected dies)
        # print('SOMME I', SumI)
        if SumI == 0:
            print('WARNING : ABORTED SIMULATION, No infected remaining')
            break
        print(sim_time, 'time')  # Kind of a loading bar


        #print('Les props', Propensities)
        #print('Les sommes',Sum_propensities)

        #Get Critical Reactions (maybe not useful since we use multinomial samples to identifiy sites where reactions occurs)
        #Criticals = fonctions.GetCriticals(Propensities, ListSites, Events)

        #We now can compute vectors mu and sigma using previous shit
        Mus = fonctions.ComputeMuNSigma(Sum_propensities, Events) # As each statechange is 1 , 0, or -1 we have sigma = mu
        #print('les mumu', Mus)

        #Get epsilon_i
        Epsis = fonctions.GetEpsilonI(SumS, SumI)
        #print('Les ei*xi', Epsis)
        #Get Tau prime
        TauPrime = fonctions.GetTauPrime(Epsis, Mus)

        # Now that main intermediary computations are done, let's get to the main algorithm Decision tree
        aox = sum(Sum_propensities)

        if TauPrime < 10/aox : # Take 10/aox 1 is left for ignoring this part
            print('Direct Method performed')
            Tau = fonctions.DoDirectMethodV2(Propensities,Sum_propensities,Nexactsteps, Events,ListSites)
        else:
            #print('Lets leap')
            #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
            #We expect that random sample of the place of reactions will be equivalent
            #As critical reactions in critical population will have low ocurrences
            Tau=TauPrime # Divided by two for a error measurement test

            #print('Voici Tau', Tau)

            #So we directly sample the kjs (number of realisation of a given event during tau) from a poisson distribution
            Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
            #print('Poisson', Poisson_means) # it's working so we are glad


            ### UPDATE 9.06.2022 : Here we sample FOR EACH site in a poisson distribuion #####
            Propensities_per_site, Sum_propensities_per_site = fonctions.GetPropensites_Per_sites(ListSites, Events)


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
                #print('site propensities',props)
                SumProp = sum(props) # We get the total propensity
                Probas = [float(i /SumProp) for i in props] # We get probability of occurence in each site
                #print('les probas', Probas) #Good job boy

                # if event.name == "Infection" or "Death I":
                #     print(event.name)
                #     print('site propensities', props)

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

                        # Here we delete the trait value corresponding to dispersing individual
                        # For virulence evolution, hold only for I indviduals
                        dispersers_traitvalues = []
                        dispersers_beta = []
                        if Site.traitvalues and abs(
                                event.Ichange) > 0:  # If there are dispersers AND that those dispersers are infected

                            for i in range(trigger_persite[index]):  # for each disperser
                                if Site.traitvalues:  # In case we remove all trait values during the loop (induces error messages)
                                    #print(type(Site.traitvalues))
                                    disperser = float(np.random.choice(
                                        Site.traitvalues))  # Get the value that is to be depleted and added to the receiving site
                                    Index_Disperser = Site.traitvalues.index(disperser)  # Get his corresponding beta
                                    dispersers_traitvalues.append(disperser)
                                    dispersers_beta.append(Site.betaI[Index_Disperser])
                                    Site.traitvalues.remove(disperser)  # Remove disperser from actual site
                                    Site.betaI.pop(Index_Disperser)  # Remove corresponding beta
                                else:  break

                        #Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
                        SuccessfulMigrants = 0
                        for i in range(nbmigrants):
                            roll4urlife = np.random.uniform(0,1,1)
                            if roll4urlife > rho : SuccessfulMigrants += 1
                        #Here we distribute successful migrants among neighboring sites
                        #This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
                        #Implemented this way, i should be easy to simulate on any network type
                        for i in range(SuccessfulMigrants):
                            if Site.traitvalues and abs(event.Ichange) > 0:
                                # Get trait values of surviving individuals
                                #print(dispersers_traitvalues)
                                SurvivorTrait = np.random.choice(
                                    dispersers_traitvalues)  # All values are chosen with same probability
                                Index_survivor = dispersers_traitvalues.index(SurvivorTrait)  # Get index
                                SurvivorBeta = dispersers_beta[Index_survivor]  # Get corresponding beta

                                dispersers_traitvalues.remove(SurvivorTrait)  # Chosen value is removed
                                dispersers_beta.remove(SurvivorBeta)  # Remove it from the list

                            index_sites = deepcopy(sites_indexes)  # working copy of site indexes vector
                            del index_sites[
                                index]  # Drop the current site from the list cause you can't emigrate to the place from which you departed

                            Index_destination = np.random.choice(index_sites)  # Get index of destination site

                            #Add individual to destination
                            if abs(event.Schange) > 0 : #if S are dispersers
                                ListSites[Index_destination].effectifS += 1
                            elif abs(event.Ichange) > 0 : # if I are dispersers
                                #print('PYCHARM WAS HERE')
                                ListSites[Index_destination].effectifI += 1
                                # Add trait value of individual
                                ListSites[Index_destination].traitvalues.append(SurvivorTrait)
                                ListSites[Index_destination].betaI.append(SurvivorBeta)
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

                            if trigger_persite[index] > 0:  # If the extinction has really occured
                                Site.traitvalues = []
                                Site.betaI = []
                        else :

                            if abs(event.Ichange) > 0:  # If the event has an effect on 'evolving population'
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange
                                # print('Avant',len(Site.traitvalues))

                                Site.traitvalues, Site.betaI = fonctions.ChooseTraitValue(Evoltrait, trigger_persite[index],
                                                                            event.Ichange,Site.traitvalues, Site.betaI)
                            else :
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange

        #Update time
        sim_time += Tau
        print('Effectifs', SumS, SumI)
        #print('time increment', Tau)
        #print('Le temps passe si vite',sim_time)

        #Update the output tracking

        # 1. Densities
        indexlist = 0
        for index, i in enumerate(ListSites):
            if i.effectifS < 0:  # Avoid negative population in the "big fat brute" way
                i.effectifS = 0
            dico_densities_df[f"S{index}"].append(i.effectifS)
            indexlist += 1
            if i.effectifI < 0:
                i.effectifI = 0
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
                dico_traits_df[f"site{index}"].append(MeanTraitValue)
                indexlist2 += 1
            else:
                MeanTraitValue = 'NA'
                list_value.append(MeanTraitValue)
                dico_traits_df[f"site{index}"].append(MeanTraitValue)
        # print(Traits_Values_Out)
        # print(SumI, SumS)
        # 4 Count the different phenotypes in the metapopulation, in order to follow their distribution over time
        indexlist3 = 0
        Possible_values = fonctions.Rounded_alpha_values
        list_out_t = []
        for i in Possible_values:  # For each value defined
            count = 0
            for j in ListSites:  # We browse the different sites
                # We count the given value and sum it
                count += j.traitvalues.count(i)
            dico_distrib_df[f"Alpha{i}"].append(count)
            list_out_t.append(count)

    ################ PREPARING THE OUTPUTS #################################################################################
    ################   Structuring outputs to get a .csv file even if the loop has broken ##################################

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



    #Creating the time series dataframe
    datadensity = pd.DataFrame.from_dict(data=dico_densities_df)
    VectimeDf = pd.DataFrame(data=vectime)
    datadensity.insert(0, "Time", VectimeDf, allow_duplicates=False)
    datadensity.to_csv('Metapop_Outputs_LongRunStep005' + str(d) + '_' + str(seed) + '.csv')
    # Creating Mean trait dataframe
    datatrait = pd.DataFrame.from_dict(data=dico_traits_df)
    datatrait.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datatrait.to_csv('Traits_outputs_LongRunStep005' + str(d) + '_' + str(seed) + '.csv')
    # Creating distribution dataframe
    datadistrib = pd.DataFrame.from_dict(data=dico_distrib_df)
    datadistrib.insert(0, 'Time', VectimeDf, allow_duplicates=False)
    datadistrib.to_csv('Distribution_outputs_LongRunStep005' + str(d) + '_' + str(seed) + '.csv')

################### MULTIPROCESSING PART ###########


# Paramètres de multiprocessing
list_seeds = [1]
list_params =[1]
nbsims = len(list_seeds)

#Lancer un batch de nbsims simulations
#Attention a pas trop déconner avec, car le multisim annule la diffusion d'eventuels messages d'erreur
if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count()-2 #nombre de processeurs, moins deux par prudence. (et pour pouvoir faire d'autres choses en meme temps)
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb) # Je ne sais pas trop ce que ça fait
    for j in range(len(list_params)) :
        for i in range(nbsims):
            pool.apply_async(RunModel, args=(list_seeds[i],list_params[j])) #Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
            #RunModel(list_seeds[i],list_params[j]) #pour debug hors multisim (messages d'ereur visible)
    pool.close()
    pool.join()