# Script that build regular network of k neighbours

# Based on adjacency list, geometry of 4-neighbor graph does not necessarily corresponds to a regular quare grid

# Set the List of sites (should be indexed by N )
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
                        dico_adjacency[f"{i}"].append(Sites[Index_current - j])
                    if Next_neighbor in List_indexes :
                        dico_adjacency[f"{i}"].append(Sites[Index_current + j]))
    return dico_adjacency

Mon_raizo = Build_Neighboring(Sites, Nb_Neighbors)
print(Mon_raizo)