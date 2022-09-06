#### SCRIPT USED (IN THE PAST) TO REDUCE THE SIZE OF OUPUTS WHEN NEEDED ####
#### USELESS : NOW PLEASE USE THE TIMESTEP SAVING IN THE MAIN PROGRAM #####

import os

import numpy as np
import pandas as pd

directory = os.getcwd() # C:\Users\Julien\PycharmProjects\Metapop_Model

#Get Files from direct method output
DirectFolder = os.path.join(directory, 'Sim Out')
os.chdir(DirectFolder)
print(DirectFolder)
ListDmDirectory = os.listdir(DirectFolder)
print(ListDmDirectory)
Dict_dfDM = {}
for index, file in enumerate(ListDmDirectory):
    if 'outputs' in file : #If the file is an output csv # If 'Output' is the good statement
        Dict_dfDM[str(file).replace('.csv','')] = pd.read_csv(str(file))
#Get Files from TauLeapMetapop outputs
#LeapFolder = os.path.join(directory, 'Leap out')
#os.chdir(LeapFolder)
#ListLeapDirectory = os.listdir(LeapFolder)
#Dict_dfLeap = {}
#for index, file in enumerate(ListLeapDirectory):
    #if 'Outputs' in file : #If the file is an output csv
        #Dict_dfLeap["LeapDataframe{0}".format(index)] = pd.read_csv(str(file))

#Create a mean dataframe of DM
#Round values of t à 3 digits (leaves 40.000 dots for RMSE computation)

#Prepare output destination
OutFolder = os.path.join(directory, 'Sim Out Reduced')
os.chdir(OutFolder)


#Fonction to find the nearest value of a given number in an array
def find_nearest(array, value):
    nearestvalue_index = (np.abs(array-value)).argmin() # return the index of the nearest vale
    return nearestvalue_index

keptvalues = [i for i in np.arange(0,1500,0.01)]
IndexesToKeepAll = []
for key, value in Dict_dfDM.items() : #Loop over all dfs, key is a name, Value is the df object

    t_serie = pd.Series(value['Time']) # Get the 't' series of a specific DF
    t_array = t_serie.to_numpy() # Convert the serie into array (needed)

    IndexesToKeep = []
    for i in keptvalues : #for each value we want
        print(i) # Keep track of 'loading' during execution
        index_to_keep = find_nearest(t_array, i) # Call the function that find nearest value from ex
        IndexesToKeep.append(index_to_keep) # List that contain the indexes to keep
    #Uptade the df by keeping only values of interest
    reducedDF = value.iloc[IndexesToKeep]
    print('Nouvelle taille supposée 150k', len(reducedDF))
    #Write a new csv file from reduced dataframe
    reducedDF.to_csv(str(key)+'_reduced.csv')



