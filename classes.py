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

class Site(): #Site object containing (non explicit) individuals
    def __init__(self,effectifS,effectifI, *args): #First try with arg way to implement feature unsure
        self.effectifS = effectifS #S density
        self.effectifI = effectifI #I density
        self.traitvalues = [] #traitvalue Affect a trait to each individual (vector) without being individual-based
        self.betaI = []

        #self.pos = pos (tuple) : for future improvements, position of the site on the network grid, as matrix coordinates
        #self.neighbor = [] : for future, maybe including neighbors as an attribute (so it's computed only once)
class EvolvingTrait():
    def __init__(self, name, IsMutation): #str, bool
        self.Traitname = name
        self.TraitMutation = IsMutation

class Event():
    def __init__(self,name, propensity, Schange, Ichange, order,EvolvingTrait):
        self.name = name # Event name in letter and not in memory address, handful to identify what's happening
        self.S = 0 #Has to take density values to understand the maths
        self.I = 0
        self.formula = propensity # The unique formule (str) given by model construction
        self.propensity = eval(self.formula)#Convert string in maths instruction, very useful to externalise model building
        self.Ichange = eval(Ichange) # State Change due to event, Typically -1, 0, 1 except for extinctions
        self.Schange = eval(Schange)
        self.order = order # Reaction order, not really useful but we never know
        self.EvolvingTrait = EvolvingTrait

    def UpdatePropensity(self, S, I, TraitValues, BetaI): # Class method to compute propensities without creating new objects
        if self.EvolvingTrait.Traitname in self.formula :
            self.S = S
            self.I = I
            if self.EvolvingTrait.Traitname == 'alpha': # if the evolving trait is virulence
                if self.name == 'Death I' :
                    SumtraitValues = sum(TraitValues)  # Do some shit to take into account the sum of trait values instead of parameter value

                    # Comparaison avec "déterministe"
                    bleu = 1.5* self.I
                    #print('Somme VS Deter :', SumtraitValues, bleu)

                    self.propensity = SumtraitValues
                if self.name == 'Infection':

                    #Here we need first to compute the BetaI for each individual as we se a trade-off between alpha and beta
                    #BetaI =fonctions.GetBetaI(TraitValues)
                    #print('COUCOU', BetaI, type(BetaI))
                    SumBetaI = sum(BetaI)
                    rouge = (0.025* 1.5 /2.5) * self.I * self.S
                    #print('Somme VS Deter II :', SumBetaI, rouge)
                    self.propensity = SumBetaI * S #Proper calculation of propensity with individual beta
        else : #Otherwise take the fixed value given
            #print(self.name)
            self.S = S
            self.I = I
            self.propensity = eval(self.formula) # Changes propensity values while keeping formula untouched
        return self.propensity