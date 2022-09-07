# DEV SCRIPT FOR SPATIAL STUFF, NOT NEEDED IN THE MAIN PROGRAM
import math

import numpy as np
from numpy import random
import pandas as pd
from copy import deepcopy
# Based on adjacency list, geometry of 4-neighbor graph does not necessarily corresponds to a regular quare grid

# Set the List of sites (should be indexed by N )
from scipy.spatial import distance_matrix

Sites = [0,1,2,3,4,5,6,7,8,9]
NbSites = len(Sites)
Nb_Neighbors = 4

print(int(math.sqrt(25)))

### CONSTRUCTION D'UN RESEAU A PARTIR DE SITES ADJACENTS DANS UNE LISTE

def Build_Accordion_Neighboring(Listsites, Nb_neighbors ) : # ListSites is a list of sites objects, Nb_neighbors is an Int
    List_indexes = []
    for i in Listsites :
        List_indexes.append(i.index)
    if Nb_neighbors %2 != 0 :
        print("Please choose an even value for the number of neighbors")
        dico_adjacency = {}
    else :
        dico_adjacency = {}
        for i in List_indexes :
            dico_adjacency[f"{i}"]= []
            Index_current = List_indexes.index(i)
            for j in range(int(Nb_neighbors/2)+1) : # /5 necessary cause the loop takes previous and next neighbourgs according to indexes / +1 necessary cause range(2) = [0,1]
                if j != 0:
                    Previous_neighbor = Index_current - j
                    Next_neighbor = Index_current + j
                    if Previous_neighbor in List_indexes :
                        dico_adjacency[f"{i}"].append(Listsites[Index_current - j])
                    if Next_neighbor in List_indexes :
                        dico_adjacency[f"{i}"].append(Listsites[Index_current + j])
    return dico_adjacency


### CONSTRUCTION D'UN RESEAU ALEATOIRE AVEC DISTANCE EUCLIDIENNES OU MATRICE DE CONNECTIVITE
Sites = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
def Set_kneigh_random_network(List_sites) :

    Liste_x = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
    Liste_y = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
               29, 30]

    data_coord = []  # Liste qui contiendra la liste des coordonnées de chaque station

    # Dictionnaire qui contient en Key l'identifiant de la station et en valeur l'objet station
    dico_coordSites = {}

    for iteration in List_sites :
        #Tirage des coordonnées géographiques du site
        x_coord_l = random.choice(Liste_x,1)
        x_coord = x_coord_l[0]
        y_coord_l = random.choice(Liste_y,1)
        y_coord = y_coord_l[0]
        xy_list = [x_coord, y_coord]
        data_coord.append(xy_list)

        #On crée la station comme valeur du dictionnaire
        dico_coordSites[f"{iteration}"] = xy_list
        #print(dico_univers)


    liste_des_sites = list(dico_coordSites.keys())
    #print('data', data)
    #print("liste des noms de stations", liste_nom_stations)

    #On crée la matrice de distances entre les stations
    df = pd.DataFrame(data_coord, columns=['xcord', 'ycord'], index=liste_des_sites)
    Distances_entre_sites = pd.DataFrame(distance_matrix(df.values, df.values), index=df.index, columns=df.index)
    #print(Distances_entre_stations)

    ### Normaliser la matrice de distances (somme des lignes fait 1)
    Dico_normalised = {}
    for i in range(len(Distances_entre_sites)):
        row_values = Distances_entre_sites.iloc[i] # Get values per row
        Dico_normalised[f"{i}"] = [] # Enable storage
        sumrow = sum(row_values) # Get the sum of row values
        for j in row_values :
            Dico_normalised[f"{i}"].append(j/sumrow) # Normalize distances
    Normalised_distance_matrix =  pd.DataFrame.from_dict(data=Dico_normalised) # And Build a DF from distances dictionnary


    ### Transformer la matrice des distances en matrice de connectivité
    ### Chacun est connecté à ses n plus proches voisins (géographiques)
    Dico_Connectivity = {}
    kneigh = 4
    # Find the indexes of neighbors
    Matrix_size = len(List_sites)
    # List to fill with the matrix
    connectivity_matrix = [] # Create an empty list

    for site in liste_des_sites :
        ActualSite = site
        #Find the distance that is zero (distance of the site to itself, must be deleted)
        Distances_entre_sites_copy = deepcopy(Distances_entre_sites) # For local modification
        indexNames = Distances_entre_sites_copy[Distances_entre_sites_copy[ActualSite] == 0].index
        Distances_entre_sites_copy.drop(indexNames, inplace=True) # Deletion

        #Select distances to the actual place in a 'table'
        Distance_to_Actual = Distances_entre_sites_copy[ActualSite]

        #Find the mins : we simply sort the values and keep the kneigh first
        sorted_distances = Distance_to_Actual.sort_values()

        # Get in list the neighbourgs and associated distances
        otherSites_indexes = sorted_distances.index.tolist()
        Distance_to_others = sorted_distances.tolist()
        my_neighbors = []
        Distance_to_neigh = []
        for i in range(kneigh) :
            my_neighbors.append(otherSites_indexes[i])
            Distance_to_neigh.append(Distance_to_others[i])

        #Create a matrix row with only zeros
        row = np.zeros((Matrix_size))
        for i in my_neighbors : # For each neighbor index found
            row[int(i)] = 1 # Set connectivity to one
        connectivity_matrix.append(row)
        Array_matrix = np.array(connectivity_matrix)





    return dico_coordSites, liste_des_sites,  Distances_entre_sites, Normalised_distance_matrix, Array_matrix

dico_coordSites, liste_des_sites,  Distances_entre_sites, Normalised_distance_matrix, Conn_Mat = Set_kneigh_random_network(Sites)

#print("1",dico_coordSites)
#print("2",liste_des_sites)
#print("3",Distances_entre_sites)
#print("4", Normalised_distance_matrix)
#print(Normalised_distance_matrix.sum()) Attention ça ne fait 1 que sur la somme des lignes (Si A est PPV de B, la réciproque marche pas toujours)
#print(Conn_Mat)
# TOUT FONCTIONNE MAWAHAHAHAHAHAH


###### Algorithm for hexagonal lattice grid (each site IS an hexagon, six neighbors) #####
nb_sites = 25
SiteListe = []
for i in range(nb_sites):
    SiteListe.append(i)

# Explanation of the algorithm
# 1. Build a square lattice grid
# 2. For each even row, each site is linked to one on the next row, by his left-down diagonal
# 3. For each odd row, each site is link to one in the next row by his right-down diagonal
# 4. Deal with the boundary conditions if you want to.
# 5. That's all

# CONSTRUCTION D'UN RESEAU REGULIER A 4 VOISINS, LATTICE CARRE

def Build_Square_Lattice(n,p): # Takes as parameters a number of sites per row (n) et per column (p)
    nb_sites = n*p
    dico_coordSites = {}
    #First we explicit each site and assign them a vector of coordinates (x,y) that is their position on the grid
    Site_counter = 0
    for i in range(n) : # For each row
        for j in range(p) : # And each column
            dico_coordSites[Site_counter]= [i,j]
            Site_counter += 1
    #Now that the nodes are in place, we add edges as if we wanted to fill a square lattice (with borders, for now)
    row_edges = Pairwise_association(dico_coordSites, "Row", n)
    col_edges = Pairwise_association(dico_coordSites, "Column", p)

    #From this we build an adjacency dictionnary, which at each site makes correspond a list of neighbors
    dico_neighborhood = {}
    for i in range(nb_sites) :
        dico_neighborhood[f"{i}"]= []
        for pair in row_edges :
            if i in pair :
                paircopy = deepcopy(pair)
                paircopy.remove(i)
                dico_neighborhood[f"{i}"].append(paircopy[0])
        for pair in col_edges :
            if i in pair :
                paircopy = deepcopy(pair)
                paircopy.remove(i)
                dico_neighborhood[f"{i}"].append(paircopy[0])


    return dico_neighborhood

def Pairwise_association(Coordinates,Type, length) : # Coordinate is a Dict object with Site = (x,y), Type is an str() that takes "Rows" or "Columns" values, length is the (int) length associated
    pairs = []
    Coordinates_copy = deepcopy(Coordinates)
    Liste_keys = list(Coordinates.keys())
    if Type == "Row":
        for i in range(len(Coordinates)-1) :
            xpos = Coordinates[i][1]  # 1 because we look on rows
            xnext = Coordinates[i+1][1]
            if xpos == xnext - 1 : #If we look at an adjancent site (row)
                pairs.append([Liste_keys[i], Liste_keys[i+1]])

    if Type == "Column":
        for i in range(len(Coordinates)-length) :
            xpos = Coordinates[i][0]  # 0 because we look on columns
            xnext = Coordinates[i+length][0]
            if xpos == xnext - 1 : #If we look at an adjancent site (column)
                pairs.append([Liste_keys[i], Liste_keys[i+length]])

    if Type != "Row" and Type != "Column":
        print("Error, only Row or Column are suitable inputs")

    return pairs


Dico = Build_Square_Lattice(5,5)
print(Dico)
print("j'ai perdu mes cles", Dico.keys())


# CONSTRUCTION D'UN RESEAU REGULIER A 6 VOISINS, HEXAGONAL LATTICE

def Build_Hexagonal_Lattice(n,p): # Takes as parameters a number of sites per row (n) et per column (p)
    nb_sites = n*p
    dico_coordSites = {}
    #First we explicit each site and assign them a vector of coordinates (x,y) that is their position on the grid
    Site_counter = 0
    for i in range(n) : # For each row
        for j in range(p) : # And each column
            dico_coordSites[Site_counter]= [i,j]
            Site_counter += 1
    #Now that the nodes are in place, we add edges as if we wanted to fill a square lattice (with borders, for now)
    row_edges = Pairwise_association(dico_coordSites, "Row", n)
    col_edges = Pairwise_association(dico_coordSites, "Column", p)

    #Complete the diagonal edges in order to get the hexagonal lattice
    Diagonal_edges = []
    List_keys = list(dico_coordSites.keys())
    List_values = list(dico_coordSites.values())
    for i in range(len(dico_coordSites)): # For each node
        coord = dico_coordSites[i] # We get its coordinate
        print("HERE COME THE X", coord[0])
        if coord[0] % 2 == 0 : # If the x component is even
            diagonal_left_neighbor = [coord[0]+1, coord[1]-1] # we get the site in (x+1, y-1)
            if diagonal_left_neighbor in dico_coordSites.values() : # if it exists
                Index_neighbor = List_values.index(diagonal_left_neighbor)
                Diagonal_edges.append([i, List_keys[Index_neighbor]])

        if coord[0] % 2 != 0 : # If the x component is odd
            diagonal_right_neighbor = [coord[0] + 1, coord[1] + 1]  # we get the site in (x+1, y-1)
            if diagonal_right_neighbor in dico_coordSites.values():  # if it exists
                Index_neighbor = List_values.index(diagonal_right_neighbor)
                Diagonal_edges.append([i, List_keys[Index_neighbor]])



    #From this we build an adjacency dictionnary, which at each site makes correspond a list of neighbors
    dico_neighborhood = {}
    for i in range(nb_sites) :
        dico_neighborhood[f"{i}"]= []
        for pair in row_edges :
            if i in pair :
                paircopy = deepcopy(pair)
                paircopy.remove(i)
                dico_neighborhood[f"{i}"].append(paircopy[0])
        for pair in col_edges :
            if i in pair :
                paircopy = deepcopy(pair)
                paircopy.remove(i)
                dico_neighborhood[f"{i}"].append(paircopy[0])
        for pair in Diagonal_edges :
            if i in pair :
                paircopy = deepcopy(pair)
                paircopy.remove(i)
                dico_neighborhood[f"{i}"].append(paircopy[0])


    return dico_neighborhood

DicoH = Build_Hexagonal_Lattice(5,5)
print(DicoH)
