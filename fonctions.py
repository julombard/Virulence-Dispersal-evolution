# Stochastic Simulation Algorithm for Metapopulation models
import classes
import numpy as np
from copy import deepcopy
import Params
#Parameters extracted from param files
beta0 = Params.beta0 # Infectious contact rate
bi = Params.bi  # Per capita net growth rate
di = Params.di # Per capita natural death rate
omega = Params.omega # Strength of density dependance on births
k = Params.k  # Carrying capacity
d = 1 # Dispersal propensity
gamma = Params.gamma  # Parasite Clearance
alpha = Params.alpha  # Parasite Virulence
rho = Params.rho  # Dispersal Cost
epsilon = Params.epsilon  # Extinction rate

Possible_alpha_values = list(np.arange(0.01,3, 0.1)) # Here we get values from 1st to 2nd by 3rd


def SetMetapop(nbsite, taillepop): #Creates sites objects containing populations
    ListSites=[] # List that will contain all sites
    for i in range(nbsite): # Creates sites, the 1st will always contain one infected and the other 0
        if i == 0:
            newsite = classes.Site(effectifS=taillepop-5, effectifI=5)
            # Assign to each initialised infected individual a trait value for alpha
            for j in range(newsite.effectifI):
                # newsite.traitvalues.append(float(np.random.uniform(0,1,1))) # For random sample from a given law, here uniform
                newsite.traitvalues.append(1.7)  # For sampling from predefined trait vector
            newsite.betaI = GetBetaI(newsite.traitvalues)


            ListSites.append(newsite)
        else:
            newsite = classes.Site(effectifS=taillepop, effectifI=0)
            for j in range(newsite.effectifI):
                #newsite.traitvalues.append(float(np.random.uniform(0,1,1)))
                newsite.traitvalues.append(1.5)
            newsite.betaI = GetBetaI(newsite.traitvalues)
            ListSites.append(newsite)
    #print(ListSites)
    return ListSites

def GetPropensites (Sites, Events): # Compute the propensities
    Propensities = []
    for i in Events: # For each event
        PropEvent =[]
        for j in Sites : # And each site
            S, I = j.effectifS, j.effectifI # Get the xi
            traitvals = j.traitvalues
            betaI = j.betaI
            # print('VALEURS DE TRAITS', traitvals, betaI)
            Prop = i.UpdatePropensity(S,I, traitvals, betaI) #Compute propensity
            PropEvent.append(Prop)
        Propensities.append(PropEvent) # List by event AND by site

    sumpropensities = []
    for i in Propensities :
        sumpropensities.append(sum(i))
    return Propensities, sumpropensities


def GetPropensites_Per_sites (Sites, Events): # Compute the propensities
    Propensities_per_sites = []
    for i in Sites: # For each event
        PropSite =[]
        S, I = i.effectifS, i.effectifI # get the xi
        traitvals = i.traitvalues
        betaI = i.betaI
        for j in Events : # And each site
            Prop = j.UpdatePropensity(S,I,traitvals, betaI) #Compute propensity
            PropSite.append(Prop)
        Propensities_per_sites.append(PropSite)
    sumpropensities = []
    for i in Propensities_per_sites :
        sumpropensities.append(sum(i))
    return Propensities_per_sites, sumpropensities

def SumDensities(Sites) : # Count all S and I individual
    SumS = 0
    SumI = 0
    for i in Sites:
        SumS += i.effectifS
        SumI += i.effectifI
    return SumS, SumI

def GetCriticals(Propensities, Sites, Events): # Not used anymore but stored in case
    Crit_treshold = 11 #Critical number of individuals, should be between 2 and 20 (Cao et.al. 2006)
    Criticals = []
    for indexi ,i in enumerate(Events) :
        CriticalEvent = []
        for indexj,j in enumerate(Sites) :
            S, I = j.effectifS, j.effectifI  # Get the xi
            if 'epsilon' in i.formula : # Case of extinction which is always critic
                CriticalEvent.append(1)
                continue
            if i.Schange < 0 : # Case where an S individual is depleted
                if S < Crit_treshold : # If S subpop is low
                    if Propensities[indexi][indexj] > 0 : # But not zero
                        CriticalEvent.append(1) # Its critical
                    else: CriticalEvent.append(0) # Its not critical
                else : CriticalEvent.append(0)
            if i.Ichange < 0 : # Case where an I individual is depleted
                if I < Crit_treshold : # If I subpop is low
                    if Propensities[indexi][indexj] > 0 : # But not zero
                        CriticalEvent.append(1) # Its critical
                    else : CriticalEvent.append(0) # Its not critical
                else: CriticalEvent.append(0)
            elif i.Schange >0 and i.Ichange == 0 : CriticalEvent.append(0) # Case reproduction S
        Criticals.append(CriticalEvent)
    return Criticals

def ComputeMuNSigma(SumPropensities , Events):
    MuS_vector = []
    MuI_vector = []
    Output_vector = []

    for index, i in enumerate(Events) :
        if i.name == 'Extinction':
            continue
        elif i.Schange < 0 : # Case where an S individual is depleted
            vij = i.Schange
            PropS = SumPropensities[index]
            MuS_vector.append(vij*PropS)
        elif i.Ichange < 0 :
            vij = i.Schange
            PropI = SumPropensities[index]
            MuI_vector.append(vij * PropI)
        elif i.Schange > 0 and i.Ichange == 0 :
            vij = 1
            PropS = SumPropensities[index]
            MuS_vector.append(vij*PropS)
    Output_vector.append(sum(MuS_vector))
    Output_vector.append(sum(MuI_vector))
    out = np.array(Output_vector)

    return out

def GetEpsilonI(SS, SI):
    #Hor_S = 3 # Highest order reaction consuming S
    #Hor_I = 2 #Idem for I
    epsilon2 = 0.03
    Output_vector = []
    gi_s = 3/2*(2+1/(SS-1)) # Given in the paper
    gi_i = 2 # same
    Output_vector.append(epsilon2/gi_s * SS)
    Output_vector.append(epsilon2/gi_i * SI)

    out = np.array(Output_vector)
    return out

def GetTauPrime(Epsis, Mus):
    Tau_candidates = []
    Tau_finalists =[]
    for i in range(len(Mus)): # Get all Tau candidates (2 per specie)
        candidates = []
        upperterm = max(Epsis[i], 1) # Numerator
        upperterm_squared = upperterm**2 # Squared

        # Cases where all infected died in the early sim, we use a very small value for mu
        # (which leads to a huge value for tau which will never be selected but allows program to continue)
        #if Mus[i] == 0 :
            #Mus[i]= 1 / pow(10, 6)
            #print('WARNING : Avoiding division by zero by doing Trickster shit')
        #Abandonned because simulation is cancelled if all infected die


        candidates.append(upperterm/ abs(Mus[i]))
        candidates.append(upperterm_squared/abs(Mus[i])) # No need to square because mu = sigma here
        Tau_candidates.append(candidates)
    #print('Tau candidates', Tau_candidates)
    for i in range(len(Tau_candidates)): # get the tau 'finalists' (1 per specie)
        Tau_finalists.append(min(Tau_candidates[i]))
    TauPrime = min(Tau_finalists) # TauPrime is the minimum of all
    #print('Tau candidats', Tau_candidates)
    #print('Tau finalistes', Tau_finalists)
    #print('and the winner is', TauPrime)

    return TauPrime

def DoDirectMethod(Propensities, Sum_propensities, exactsteps, events, sites):
    rho = Params.rho  # Not very convenient to put it here but it has to....

    for i in range(exactsteps):
        r1 = np.random.uniform(0,1,1) #Draw random numbers
        a0 = sum(Sum_propensities) # Overall propensity

        Tau = (1/a0) * np.log(1/r1) #Time increment
        Probas = [i / a0 for i in Sum_propensities] #Individual propensity
        NextReaction = list(np.random.multinomial(1, Probas)) # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        #print('The next reaction Index', NextReactionIndex)

        #Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        sumprop = sum(props)
        proba_site = [i/sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1)
        #print('Next places', NextPlaceIndex)


        # This part apply the effect of events in site populations
        event = events[NextReactionIndex]
        site = sites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
            nbmigrants = 1
            # Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
            SuccessfulMigrants = 0

            roll4urlife = np.random.uniform(0, 1, 1)
            if roll4urlife > Params.rho: SuccessfulMigrants += 1
            # Here we distribute successful migrants among neighboring sites
            # This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
            if SuccessfulMigrants == 1:
                # Determine which site will receive the dispersing individual
                receiving_sites = deepcopy(sites)  # Create a working copy of sites
                # print('receivers', receiving_sites)
                del receiving_sites[NextPlaceIndex]  # removing departure site from the copy
                # print('receivers post suppression', receiving_sites)
                site_destination = np.random.choice(receiving_sites)  # destination is a site object
                # print('The destination is', site_destination)

                # add individual to destination
                if abs(event.Schange) > 0:  # if S are dispersers
                    site_destination.effectifS += 1
                elif abs(event.Ichange) > 0:
                    site_destination.effectifI += 1
                else:
                    pass
                    #print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            site.effectifS +=  event.Schange
            site.effectifI +=  event.Ichange
        return Tau


def DoDirectMethodV2(Propensities, Sum_propensities, exactsteps, events, sites, index_sites):
    rho = Params.rho  # Not very convenient to put it here but it has to....
    for i in range(exactsteps):
        r1 = np.random.uniform(0, 1, 1)  # Draw random numbers
        a0 = sum(Sum_propensities)  # Overall propensity

        Tau = (1 / a0) * np.log(1 / r1)  # Time increment
        Probas = [i / a0 for i in Sum_propensities]  # Individual propensity
        NextReaction = list(np.random.multinomial(1, Probas))  # List to get index easily
        NextReactionIndex = NextReaction.index(1)
        # print('The next reaction Index', NextReactionIndex)

        # Determine where the reaction is going to happen
        props = Propensities[NextReactionIndex]  # We get the propensity per sites
        sumprop = sum(props)
        proba_site = [i / sumprop for i in props]
        NextPlace = list(np.random.multinomial(1, proba_site))
        NextPlaceIndex = NextPlace.index(1)
        # print('Next places', NextPlaceIndex)

        # This part apply the effect of events in site populations
        event = events[NextReactionIndex]
        site = sites[NextPlaceIndex]
        if 'Dispersal' in event.name:
            # Multiply the state change in population by the number of triggers
            site.effectifS += event.Schange
            site.effectifI += event.Ichange
            nbmigrants = 1
            # Here we apply dispersal cost to determine the number of successful migrants, rho is defined at the top
            SuccessfulMigrants = 0

            roll4urlife = np.random.uniform(0, 1, 1)
            if roll4urlife > Params.rho: SuccessfulMigrants += 1
            # Here we distribute successful migrants among neighboring sites
            # This part can be improved as neighboring rules become more complex, using a specific class 'network' to determine the neighbors
            if SuccessfulMigrants == 1:

                sites_index = deepcopy(index_sites)  # working copy of site indexes vector
                del index_sites[
                    NextPlaceIndex]  # Drop the current site from the list cause you can't emigrate to the place from which you departed

                Index_destination = np.random.choice(index_sites)  # Get index of destination site

                # add individual to destination
                if abs(event.Schange) > 0:  # if S are dispersers
                    sites[Index_destination].effectifS += 1
                elif abs(event.Ichange) > 0:
                    sites[Index_destination].effectifI += 1
                else:
                    pass
                    # print('There was only one migrant, but he died. Nothing happend')
        else:
            # Multiply the state change in population by the number of triggers
            # Multiply the state change in population by the number of triggers
            if event.name == 'Extinction':
                # When extinction occur we need to retrieve densities values cause they're initialized to zero otherwise in class definition
                event.Schange = -site.effectifS  # I dont really understand why attribute is not protected but it's good news
                event.Ichange = -site.effectifI
                # print('Event Schange', trigger_persite[index])
                site.effectifS += event.Schange
                # print('Effectif après extinction', Site.effectifS)
                site.effectifI += event.Ichange
            else:
                site.effectifS += event.Schange
                site.effectifI += event.Ichange
        return Tau

def GetBetaI(traitvalues): # At initialization ONLY , When evolution with trade-off transmission virulence, get beta values from alpha
    IndividualValues = [] #Outlist
    if traitvalues : # If Not empty
        for i in range(len(traitvalues)): #For each value of alpha
            value = beta0 * traitvalues[i] / (1+ (traitvalues[i])) # Get an individual value of beta
            IndividualValues.append(value) # store them
    return IndividualValues


def ChooseTraitValue(EvolvingTrait,NbTrigger,Statechange, Traitvalues, BetaI) : #Function that choose wich individual is added / depleted during evolutive events (DOES NOT HOLD FOR DISPERSAL)
    newtraitsvalues = deepcopy(Traitvalues) # This line is probably useless, but the safer the wiser
    newBetaI = deepcopy(BetaI)
    for i in range(NbTrigger):
        if newtraitsvalues : # If the trait values vector is not empty (important : random sampling can't occur in empty objects)
            if Statechange > 0:  # If it's a birth
                if EvolvingTrait.TraitMutation == True : # If Mutation is allowed in simulation
                    probamut = 0.005 # We set a mutation probability (à passer en param global ou en attribut de classe)
                    roll2mutate = np.random.uniform(0,1,1)
                    if roll2mutate < probamut : #Here there is mutation
                        #Newalpha = float(np.random.uniform(0,1,1)) # sample the new alpha value
                        #Newalpha = float(np.random.choice(Possible_alpha_values,1, replace=True)) # To get Alpha in vector of discrete given values
                        #NewBeta = beta0 * Newalpha / (1+Newalpha) #Get the corresponding betai

                        #newtraitsvalues.append(Newalpha)
                        #newBetaI.append(NewBeta)
                        # sample the individual that gives birth
                        SumbetaI = sum(newBetaI)  # Get total reproduction propensity
                        IndividualsProbas = [i / SumbetaI for i in newBetaI]  # Get Individual probability to reproduce

                        Reproducer = list(
                            np.random.multinomial(1, IndividualsProbas))  # Choose one ( this return a beta_i)
                        Index_reproducer = Reproducer.index(1)  # Get the index of the reproducer

                        Parent_Value = newtraitsvalues[Index_reproducer] # We get the trait value of the parent

                        ChangeMut = [-0.1,0.1]
                        NewValue = Parent_Value + np.random.choice(ChangeMut, 1)
                        print('Nouvelle valeur après mutation')
                        newtraitsvalues.append(NewValue)
                        NewbetaValue = beta0 * NewValue / ( NewValue +1)
                        newBetaI.append(NewbetaValue)

                    else :
                        # sample the individual that gives birth
                        SumbetaI = sum(newBetaI) # Get total reproduction propensity
                        IndividualsProbas = [float(i / SumbetaI) for i in newBetaI] # Get Individual probability to reproduce

                        Reproducer = list(np.random.multinomial(1, IndividualsProbas)) # Choose one ( this return a beta_i)
                        Index_reproducer = Reproducer.index(1) # Get the index of the reproducer

                        newtraitsvalues.append(newtraitsvalues[Index_reproducer])
                        newBetaI.append(newBetaI[Index_reproducer])

                else : # If mutation not allowed
                    # sample the individual that gives birth
                    #print('coucou',Traitvalues)
                    SumbetaI = sum(newBetaI)  # Get total reproduction propensity

                    IndividualsProbas = [float(i / SumbetaI) for i in newBetaI]  # Get Individual probability to reproduce
                    Reproducer = list(np.random.multinomial(1, IndividualsProbas))  # Choose one ( this return a beta_i)

                    Index_reproducer = Reproducer.index(1)  # Get the index of the reproducer
                    newtraitsvalues.append(newtraitsvalues[Index_reproducer])
                    newBetaI.append(newBetaI[Index_reproducer])


            elif Statechange < 0:  # If it's a death

                SumTrait = sum(newtraitsvalues)  # We get the total of traits values
                IndividualsProbas = [float(i / SumTrait) for i in newtraitsvalues]  # We get individuals probabilities

                DepletedGuy = list(np.random.multinomial(1, IndividualsProbas))  # Sample the removed individual
                Index_depleted = DepletedGuy.index(1) #Which value was chosen
                newtraitsvalues.pop(Index_depleted)  # And remove it
                newBetaI.pop(Index_depleted) # and his corresponding beta value


    return newtraitsvalues, newBetaI