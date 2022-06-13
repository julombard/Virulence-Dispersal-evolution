# Stochastic Simulation Algorithm for Metapopulation models
from copy import deepcopy
import classes
import fonctions
import numpy as np
import pandas as pd
import multiprocessing
import Params

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
#A story inspired by Modified Poisson Tau leap algorithm from cao et. al. (2006)
#Including features for metapopulations modelling, designed by Massol F., Lion S. and bibi

#Damn efficient
#Beware fun but useless commentaries (this one doesn't count)
def RunModel(seed, param) :

    #This part changes parameters value for long autonomous runs
    d = param
    classes.d = param
    fonctions.d = param

    # Stocker les sorties dans un dictionnaire
    dico_density_df = {}
    dico_distrib_df = {}
    dico_traits_df = {}

    #Simulation parameters
    # Some seeds doesnt work cause the 1st infected guy dies before infecting anyone
    # seed = 0
    print('Seed', seed)
    np.random.seed(seed) #Set seed for reproducibility

    sim_time = 0 # Simulation time (model time, not an iteration number)
    vectime = [0] # to keep t variable
    tmax = 100 # Ending time
    Nexactsteps = 20  # Number of steps to do if/when performing direct method
    # As Tauprime decreases less faster than 1/ax0 when adding event/site, DM is not performed anymore above a certain number of sites (6 ?)

    nbsite = 40 # Number de sites
    Taillepop = Params.k # Initial local population sizes

    #Define a trait which will be evolving
    #Evoltrait = classes.EvolvingTrait('alpha', True) # Name of the trait (str), Is mutation Enabled (bool)
    Evoltrait = classes.EvolvingTrait('None', False)

    #Define population as class instances
    ListSites = fonctions.SetMetapop(nbsite, Taillepop)

    #Event definition
    #Further expansion idea : build events from a unique model.txt file read by the program, in order to simulate whathever you want
    ReproductionS = classes.Event(name='Reproduction S',propensity='bi*(1- omega * (self.S+self.I)/k) * self.S', Schange='1', Ichange='0', order=1,EvolvingTrait=Evoltrait)
    DeathS = classes.Event(name='Death S',propensity='di*self.S', Schange='-1', Ichange='0', order=3,EvolvingTrait=Evoltrait)
    DispersalS = classes.Event(name='Dispersal S',propensity='d*self.S', Schange='-1', Ichange='0', order=1,EvolvingTrait=Evoltrait)
    DispersalI = classes.Event(name='Dispersal I',propensity='d*self.I', Schange='0', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    Extinction = classes.Event(name='Extinction',propensity='epsilon', Schange='-self.S', Ichange='-self.I', order=0,EvolvingTrait=Evoltrait)
    Infection = classes.Event(name='Infection',propensity='(beta0 * alpha / (1+ alpha )) *self.S*self.I', Schange='-1', Ichange='1', order=2,EvolvingTrait=Evoltrait)
    Recovery = classes.Event(name='Recovery',propensity='gamma*self.I', Schange='1', Ichange='-1', order=1,EvolvingTrait=Evoltrait)
    DeathI = classes.Event(name='Death I',propensity='alpha*self.I', Schange='0', Ichange='-1', order=1, EvolvingTrait=Evoltrait)

    #Event vector, cause tidying up things is always nice
    Events = [ReproductionS, Infection, DispersalI, DispersalS, DeathS, DeathI, Recovery, Extinction]

    #Initializing outputs storage
    Densities_out = [] # Collect densities outputs
    Propensities_out =[] # Collect propensities outputs
    Traits_Values_Out = []
    Distrib_out = []
    IsTrackPropensites = False #Set false/ true to not track/track propensities
    for i in range(len(Events)):
        Propensities_out.append([0]) #Store times series of propensities


    #Get initial lists and values in outputs
    #We want to get one list per Sx(t) and Ix(t) to store them easily in a dataframe at the end
    list_trait_site = []
    for i in ListSites:
        Densities_out.append([i.effectifS])
        Densities_out.append([i.effectifI])
        if i.effectifI > 0 :
            list_trait_site.append(sum(i.traitvalues)/i.effectifI)
            dico_traits_df[f"site{i}"] = sum(i.traitvalues) / i.effectifI
        else :
            list_trait_site.append('NA')
            dico_traits_df[f"site{i}"] = "NA"
    Traits_Values_Out.append(list_trait_site)

    #Same for individual values for distributions outputs
    list_distrib_init = []
    for i in fonctions.Possible_alpha_values:  # For each value defined
        count = 0
        for j in ListSites:  # We browse the different sites
            # We count the given value and sum it
            count += j.traitvalues.count(i)
        list_distrib_init.append(count)
        dico_distrib_df[f"site{i}"] = count
    Distrib_out.append(list_distrib_init)


    ############################# Model main Loop ######################################################################
    while sim_time < tmax :
        # Compute the propensities
        Propensities, Sum_propensities = fonctions.GetPropensites(ListSites,
                                                                  Events)  # Get a vector of propensities ordered by event and by sites
        SumS, SumI = fonctions.SumDensities(ListSites)  # Get total densities of species

        # Break the main loop if there are no infected remaining ( this happens essentially at start if the 1st infected dies)
        # print('SOMME I', SumI)
        if SumI == 0:
            print('WARNING : ABORTED SIMULATION, No infected remaining')
            break

        print( sim_time,'time') # Kind of a loading bar
        vectime.append(sim_time) #Update time vector

        #print('Les props', Propensities)
        #print('Les sommes',Sum_propensities)

        #Get Critical Reactions
        # maybe not useful since we use multinomial samples to identifiy sites where reactions occurs and number of reactions
        #In practice, we use brutal reset of population if becoming negative (barbare but working)
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
            #Evolution not implemented in direct method, because direct method is never performed with large number of sites
            print('Direct Method performed')
            Tau = fonctions.DoDirectMethodV2(Propensities,Sum_propensities,Nexactsteps, Events,ListSites)
        else:
            #print('Lets leap')
            #Here we do not compute TauPrimePrime to determine how much critical reactions occurs
            #We expect that random sample of the place of reactions will be equivalent and it is (almost)
            #As critical reactions in critical population will have low ocurrences
            Tau=TauPrime # You can divide this by whatever int you want, the more you divide ,the more it should be accurate (and time consuming)

            #print('Voici Tau', Tau)
            #So we directly sample the kjs (number of realisation of a given event during tau) from a poisson distribution
            Poisson_means = np.multiply(Tau,np.array(Sum_propensities)) #Get ajx * tau from which we will sample the kjs
            #print('Poisson', Poisson_means) # it's working so we are glad

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
                #print(event.name)
                if event.name == 'Infection' :
                    print(props)
                if event.name == 'DeathI' :
                    print(props)

                #print('site propensities',props)
                SumProp = sum(props) # We get the total propensity
                #print('SOMME', props)
                Probas = [i /SumProp for i in props] # We get probability of occurence in each site
                #print('les probas', Probas) #Good job boy

                if Noccur == 0 : #Case where the event can't happen
                    trigger_persite=[0 for i in range(nbsite)]
                else : # Normal cases
                    trigger_persite = np.random.multinomial(Noccur, Probas)
                #print('Occurrences per sites', trigger_persite)

                #The following part apply the effect of events in site populations
                for index, Site in enumerate(ListSites) :
                    if 'Dispersal' in event.name :
                        # Multiply the state change in population by the number of triggers
                        Site.effectifS += trigger_persite[index] * event.Schange
                        Site.effectifI += trigger_persite[index] * event.Ichange
                        nbmigrants = max(abs(trigger_persite[index] * event.Schange), abs(trigger_persite[index] * event.Ichange))
                        #print('Nombre de migrants', nbmigrants)

                        #Here we delete the trait value corresponding to dispersing individual
                        #For virulence evolution, hold only for I indviduals
                        dispersers_traitvalues = []
                        dispersers_beta = []
                        if Site.traitvalues and abs(event.Ichange) >0 : # If there are dispersers AND that those dispersers are infected

                            for i in range(trigger_persite[index]): #for each disperser
                                if Site.traitvalues : # In case we remove all trait values during the loop (induces error messages)

                                    disperser= np.random.choice(Site.traitvalues) #Get the value that is to be depleted and added to the receiving site
                                    Index_Disperser = Site.traitvalues.index(disperser) # Get his corresponding beta
                                    dispersers_traitvalues.append(disperser)
                                    dispersers_beta.append(Site.betaI[Index_Disperser])
                                    Site.traitvalues.remove(disperser) # Remove disperser from actual site
                                    Site.betaI.pop(Index_Disperser) # Remove corresponding beta
                                else : break
                        #print('Dispersers Trait values', dispersers_traitvalues)
                        #Here we apply dispersal cost to determine the number of successful migrants, where rho is a parameter defined at the top
                        SuccessfulMigrants = 0
                        for i in range(nbmigrants):
                            roll4urlife = np.random.uniform(0,1,1)
                            if roll4urlife > Params.rho : SuccessfulMigrants += 1
                        #Here we distribute successful migrants among neighboring sites
                        #print('Valeurs TRait dispersers', dispersers_traitvalues)
                        #print('Nombre de non morts', SuccessfulMigrants)
                        #This part can be improved as neighboring rules become more complex, using a specific class 'network' or a site attribute 'position' on a grid (matrix) to determine the neighbors
                        #Implemented this way, i should be easy to simulate on any network type, fingers crossed

                        for i in range(SuccessfulMigrants) :
                            if Site.traitvalues and abs(event.Ichange) > 0 :
                                #Get trait values of surviving individuals
                                SurvivorTrait = np.random.choice(dispersers_traitvalues) #All values are chosen with same probability
                                Index_survivor = dispersers_traitvalues.index(SurvivorTrait)  # Get index
                                SurvivorBeta = dispersers_beta[Index_survivor] # Get corresponding beta

                                dispersers_traitvalues.remove(SurvivorTrait)  # Chosen value is removed
                                dispersers_beta.remove(SurvivorBeta) # Remove it from the list

                            index_sites = deepcopy(sites_indexes) # working copy of site indexes vector
                            del index_sites[index]#Drop the current site from the list cause you can't emigrate to the place from which you departed

                            Index_destination = np.random.choice(index_sites) #Get index of destination site

                            #Add individual to destination
                            if abs(event.Schange) > 0 : #if S are dispersers
                                ListSites[Index_destination].effectifS += 1
                            elif abs(event.Ichange) > 0 : # if I are dispersers
                                #print('PYCHARM WAS HERE')
                                ListSites[Index_destination].effectifI += 1
                                # Add trait value of individual
                                #if Site.traitvalues : WHY ONLY IF TRAITSVALUES ????
                                ListSites[Index_destination].traitvalues.append(SurvivorTrait)
                                ListSites[Index_destination].betaI.append(SurvivorBeta)
                            else : print('ERROR : disperser is neither S nor I and that is very curious !') #The error never raises, good point
                            #Keep this error message cause it makes me think i'm a better programmer doing so
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

                            if trigger_persite[index] > 0 : #If the extinction has really occured
                                Site.traitvalues = []
                                Site.betaI = []
                        else :

                            if abs(event.Ichange) > 0: # If the event has an effect on 'evolving population'
                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange
                                #print('Avant',len(Site.traitvalues))

                                Site.traitvalues, Site.betaI = fonctions.ChooseTraitValue(Evoltrait,trigger_persite[index], event.Ichange,
                                                                                  Site.traitvalues, Site.betaI)
                            else:

                                Site.effectifS += trigger_persite[index] * event.Schange
                                Site.effectifI += trigger_persite[index] * event.Ichange


        #Update time
        sim_time += Tau
        #print('time increment', Tau)
        #print('Le temps passe si vite',sim_time)

        #Update the output tracking : THIS PART CAN BE OPTIMIZED BY USING A SINGLE LOOP FOR EVERYTHING
        #1. Densities
        indexlist = 0
        for i in ListSites :
            if i.effectifS < 0 : #Avoid negative population in the "big fat brute" way
                print('BOURINNAGE')
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
        #3. Trait Values
        indexlist2=0
        list_value = []
        for i in ListSites :
            if i.effectifI > 0 :

                SumTraitValues = sum(i.traitvalues)
                MeanTraitValue = SumTraitValues / i.effectifI
                list_value.append(MeanTraitValue)
                indexlist2 += 1
            else :
                MeanTraitValue = 'NA'
                list_value.append(MeanTraitValue)
        Traits_Values_Out.append(list_value)
        #print(Traits_Values_Out)
        #print(SumI, SumS)
        #4 Count the different phenotypes in the metapopulation, in order to follow their distribution over time
        indexlist3 = 0
        Possible_values = fonctions.Possible_alpha_values
        list_out_t = []
        for i in Possible_values: # For each value defined
            count = 0
            for j in ListSites : # We browse the different sites
                #We count the given value and sum it
                 count += j.traitvalues.count(i)
            list_out_t.append(count)
        Distrib_out.append(list_out_t)


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
        dataprop.to_csv('Propensities_outputs_Mutation_EDM'+str(d)+'_'+str(seed)+'.csv')


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
    data.to_csv('Metapop_OutputsOK_Mutation_EDM'+str(d)+'_'+str(seed)+'.csv')
    # NB : use relative path if output directory needs to be change


    #Creating Trait values DF
    #replace Site vector of traits values by a mean trait value
    datatrait = pd.DataFrame(data=dico_traits_df, index=vectime)
    # datatrait = pd.DataFrame(columns=['t'])
    # for i in range(nbsite):
    #     Colname = 'Site' + str(i)
    #     datatrait[Colname] = []
    #
    #
    # for i in range(len(vectime)) : # Fill df row by row
    #     timelist = [vectime[i]]
    #     newline_toadd = timelist+Traits_Values_Out[i]
    #     datatrait.loc[i] = newline_toadd
    # print(datatrait)
    datatrait.to_csv('Traits_outputsTDoff_MutationOK_EDM'+str(d)+'_'+str(seed)+'.csv')

    #Creating distribution dataframe
    datadistrib = pd.DataFrame(data=dico_distrib_df, index=vectime)
    # datadistrib = pd.DataFrame(columns=['t'])
    # for i in range(len((fonctions.Possible_alpha_values))):
    #     Colname = 'Alpha' + str(fonctions.Possible_alpha_values[i])
    #     datadistrib[Colname] = []
    # for i in range(len(vectime)):  # Fill df row by row
    #     timelist = [vectime[i]]
    #     newline_toadd = timelist + Distrib_out[i]
    #     datadistrib.loc[i] = newline_toadd
    print(datadistrib)
    datadistrib.to_csv('Distribution_outputsTDoffOK_Mutation_EDM' + str(d) + '_' + str(seed) + '.csv')


################### MULTIPROCESSING PART ###########


# Paramètres de multiprocessing
list_seeds = [1,2,3,4,5,6]
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