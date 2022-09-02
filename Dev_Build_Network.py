# Script that build regular network of k neighbours
from numpy import random
import pandas as pd
# Based on adjacency list, geometry of 4-neighbor graph does not necessarily corresponds to a regular quare grid

# Set the List of sites (should be indexed by N )
from scipy.spatial import distance_matrix

Sites = [0,1,2,3,4,5,6,7,8,9]
NbSites = len(Sites)
Nb_Neighbors = 4

def Build_Neighboring(Listsites, Nb_neighbors ) : # ListSites is a list of sites objects, Nb_neighbors is an Int
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

#Mon_raizo = Build_Neighboring(Sites, Nb_Neighbors)
#print(Mon_raizo)



### FONCTION QUI FAIT UN RESEAU ALEATOIRE AVEC POSSIBLES DISTANCES GEOGRAPHIQUES
Sites = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
def Set_random_network(List_sites) :

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
        row_values = Distances_entre_sites.iloc[i]
        Dico_normalised[f"{i}"] = []
        sumrow = sum(row_values)
        for j in row_values :
            Dico_normalised[f"{i}"].append(j/sumrow)
    Normalised_distance_matrix =  pd.DataFrame.from_dict(data=Dico_normalised)


    ### Transformer la matrice des distances en matrice de connectivité
    ### Chacun est connecté à ses n plus proches voisins (géographiques)
    Dico_Connectivity = {}
    n = 1

    # Copier la matrices des distances pour la bidouiller localement
    Distances_entre_sites_copy = deepcopy(Distances_entre_sites)



    return dico_coordSites, liste_des_sites,  Distances_entre_sites, Normalised_distance_matrix

dico_coordSites, liste_des_sites,  Distances_entre_sites, Normalised_distance_matrix = Set_random_network(Sites)

#print("1",dico_coordSites)
#print("2",liste_des_sites)
print("3",Distances_entre_sites)
print("4", Normalised_distance_matrix)
print(Normalised_distance_matrix.sum())



#### FIND k-order nearest neighbor
stationactu = position
    #print(stationactu)

    # Trouver la ligne de la station actuelle où la valeur de distance vaut zéro (car distance avec elle-même)
    indexNames = Distances_entre_stations_copy[Distances_entre_stations_copy[stationactu] == 0].index
    #print(indexNames)
    # Supprimer la ligne du DataFrame
    Distances_entre_stations_copy.drop(indexNames, inplace=True)
    #print(Distances_entre_stations_copy)
    # Sélection d'un sous tableau avec les valeurs d'intérets (distances au point actuel)
    Distance_a_stationactuelle = Distances_entre_stations_copy[stationactu]
    #print(Distance_a_stationactuelle)
    # Triage des valeurs dans l'ordre croissant
    sorted_distances = Distance_a_stationactuelle.sort_values()
    #print(sorted_distances)

    # Sort en liste les noms des voisins et leurs distances associées
    nom_Voisins = sorted_distances.index.tolist()
    distance_voisins = sorted_distances.tolist()
    #print(nom_Voisins)

    # pour un voyage possible au plus proche voisin d'ordre n
    n = 3
    destination = {}
    for i in range(n):
        #print(i, ". Vous pouvez voyager vers", nom_Voisins[i], "qui se trouve à une distance de",
              #distance_voisins[i], "UA.")
        text = "Voyage vers {} distance {} / station {} ".format(nom_Voisins[i],distance_voisins[i],dico_de_merde[nom_Voisins[i]].type)

        texte_affichage.append(text)
        destination[nom_Voisins[i]]= distance_voisins[i]