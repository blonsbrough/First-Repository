#!/usr/bin/env python
# coding: utf-8

# In[1]:


import uproot3 as uproot
import hist
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import os
import subprocess
import pandas as pd


# In[2]:


def Event_Combination(Input_Directory, Background = False, Slice = False):
    Directories = os.listdir(Input_Directory)
    TotalEvents = 0
    HistogramArray = []
    PathArray = []
    Branches = []
    CrossSections = []
    MET = []
    j1PT = []
    mjj = []
    j1Eta = []
    j1Phi = []
    j2PT = []
    j2Eta = []
    j2Phi = []
    weight = []
    Scales = []
    MET_Test = []
    Output = {}

    #Test Directories to see if they actually contain a valid histogram file
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        if os.path.exists(Input_Directory+item+"/analysis/histograms.root") != True and os.path.exists(Input_Directory+item+"/analysis/SimpleAna.root") != True:
            Directories.remove(item)
            print("Error, Histogram not found for"+Input_Directory+item+"/analysis/histograms.root")
    #Test Directories to see if they contain a valid Cross Section Output
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        if type(float(CrossSectionOutput.stdout)) != float:
            Directories.remove(item)
            Statement = "File "+item+" Unable to be combined, could not find Cross Section"
            print(Statement)
    #Ignore any background directories not of the current slice
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        if Slice == 1:
            if "1000_4000" in item:
                a=1
                print("Graph for Slice 1, 1000-4000")
            else:
                Directories.remove(item)
        elif Slice == 2:
            if "4000_7000" in item:
                a=1
                print("Graph for Slice 2, 4000-7000")
            else:
                Directories.remove(item)
        elif Slice == 3:
            if "7000_10000" in item:
                a=1
                print("Graph for Slice 3, 7000-10000")
            else:
                Directories.remove(item)
            
        elif Slice == 4:
            if "10000_-1" in item:
                a=1
                print("Graph for Slice 4, 10000--1")
            else:
                Directories.remove(item)
        elif Slice == False:
            a=1
            print("Graph for no slices")
    #Add cross sections and valid histograms to files
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        CrossSections.append(float(CrossSectionOutput.stdout))
        if os.path.exists(Input_Directory+item+"/analysis/histograms.root") == True:
            HistogramArray.append(Input_Directory+item+"/analysis/histograms.root")
        elif os.path.exists(Input_Directory+item+"/analysis/SimpleAna.root") == True:
            HistogramArray.append(Input_Directory+item+"/analysis/SimpleAna.root")
    #Apply a mask cut and add remaining events to output arrays
    for item in HistogramArray:
        PathArray.append(uproot.open(item)['allev/hftree'])
    for item in PathArray:
        Branches.append(item.arrays())
    for item in Branches:
        mask = (item[b"mjj"] > 1000)&(item[b"MET"] > 200)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        for element in item[b"MET"][mask]:
            MET.append(element)
        for element in item[b"j1PT"][mask]:
            j1PT.append(element)
        for element in item[b"mjj"][mask]:
            mjj.append(element)
        for element in item[b"j1Eta"][mask]:
            j1Eta.append(element)
        for element in item[b"j1Phi"][mask]:
            j1Phi.append(element)
        for element in item[b"j2PT"][mask]:
            j2PT.append(element)
        for element in item[b"j2Eta"][mask]:
            j2Eta.append(element)
        for element in item[b"j2Phi"][mask]:
            j2Phi.append(element)
    #Take the weights and scale them by number of inputs
    i=0 
    for item in Branches:
        mask = (item[b"mjj"] > 1000)&(item[b"MET"] > 200)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        scale = (CrossSections[i]/sum(item[b"weight"]))
        i += 1
        Scales.append(scale)
        Scaled_weight = item[b"weight"][mask]*scale/len(Directories)
        for item in Scaled_weight:
            weight.append(item)
    TotalEvents = len(MET)
    #Use this Output function to implement and variables you want to calculate from this Combination Function
    Output = {"Directories":Directories, "Number of Events":TotalEvents,
              "Cross Sections": CrossSections, "MET": MET, "j1PT":j1PT, 
              "mjj":mjj, "j1Eta":j1Eta, "j1Phi":j1Phi, "j2PT":j2PT,
              "j2Eta":j2Eta, "j2Phi":j2Phi, "weight":weight}
    MET_ModifiedBranches = []
    for item in Branches:
        mask = (item[b"mjj"] > 1000)&(item[b"MET"] > 200)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        masked = item[b"MET"][mask]
        MET_ModifiedBranches.append(masked)
    Weights = []
    for item in Branches:
        mask = (item[b"mjj"] > 1000)&(item[b"MET"] > 200)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        maskedbackground = item[b"weight"][mask]
        Weights.append(maskedbackground)
    for i in range(len(MET_ModifiedBranches)):
        plt.hist(MET_ModifiedBranches[i],bins=100,range=(0,2000),weights = Weights[i]*Scales[i], alpha=0.05, color='black')
    plt.hist(MET,bins=100,range=(0,2000),weights = weight, alpha=0.3, color='purple')
    plt.yscale('log')
    plt.title("Met Comparison, Inputs (Grey) vs Averaged Output (Purple)")
    plt.show()
    return(Output)


# In[3]:


def Background(Input_Directory):
    Directories = []
    TotalEvents = 0
    CrossSections = []
    MET = []
    j1PT = []
    mjj = []
    j1Eta = []
    j1Phi = []
    j2PT = []
    j2Eta = []
    j2Phi = []
    weight = []
    Output = {}
    Slices = [1, 2, 3, 4]
    for item in Slices:
        partition = Event_Combination(Input_Directory,True, item)
        for item in partition["Directories"]:
            Directories.append(item)
        TotalEvents += partition["Number of Events"]
        for item in partition["Cross Sections"]:
            CrossSections.append(item)
        for item in partition["MET"]: 
            MET.append(item)
        for item in partition["j1PT"]: 
            j1PT.append(item)
        for item in partition["mjj"]: 
            mjj.append(item)
        for item in partition["j1Eta"]: 
            j1Eta.append(item)
        for item in partition["j1Phi"]: 
            j1Phi.append(item)
        for item in partition["j2PT"]: 
            j2PT.append(item)
        for item in partition["j2Eta"]: 
            j2Eta.append(item)
        for item in partition["j2Phi"]: 
            j2Phi.append(item)
        for item in partition["weight"]: 
            weight.append(item)
        
    Output = {"Directories":Directories, "Number of Events":TotalEvents,
              "Cross Sections": CrossSections, "MET": MET, "j1PT":j1PT, 
              "mjj":mjj, "j1Eta":j1Eta, "j1Phi":j1Phi, "j2PT":j2PT,
              "j2Eta":j2Eta, "j2Phi":j2Phi, "weight":weight}
    return(Output)


# In[4]:


#Find Background Sums.
SignalDirectory = "/data/users/jupyter-blonsbro/SUSY/Generations/13TeV/Signal/150mjj/"
EWKBackgroundDirectory = "/data/users/jupyter-blonsbro/SUSY/Generations/13TeV/Background/EWKBackground/"
QCDBackgroundDirectory = "/data/users/jupyter-blonsbro/SUSY/Generations/13TeV/Background/QCDBackground/"
A = Event_Combination(SignalDirectory)
B = Background(EWKBackgroundDirectory)
C = Background(QCDBackgroundDirectory)


# In[5]:


#Define Variables to Graph
S_j1PT = A["j1PT"]
B_j1PT = np.append(B["j1PT"],C["j1PT"])
S_weight = A["weight"]
B_weight = np.append(B["weight"],C["weight"])
S_MET = A["MET"]
B_MET = np.append(B["MET"],C["MET"])
S_mjj = A["mjj"]
B_mjj = np.append(B["mjj"],C["mjj"])
S_j1Eta = A["j1Eta"]
B_j1Eta = np.append(B["j1Eta"],C["j1Eta"])
S_j1Phi = A["j1Phi"]
B_j1Phi = np.append(B["j1Phi"],C["j1Phi"])
S_j2PT = A["j2PT"]
B_j2PT = np.append(B["j2PT"],C["j2PT"])
S_j2Eta = A["j2Eta"]
B_j2Eta = np.append(B["j2Eta"],C["j2Eta"])
S_j2Phi = A["j2Phi"]
B_j2Phi = np.append(B["j2Phi"],C["j2Phi"])


# In[6]:


##Plot Graphs
##Feel Free to define these variables above and use these graph presets as a guide.
#MET
plt.hist(S_MET, bins=100, range=(200,3000), weights=S_weight, alpha=0.5, color = 'purple')
plt.hist(B_MET, bins=100, range=(200,3000), weights=B_weight, alpha=0.5, color = 'black')
plt.yscale('log')
plt.title("MET")
plt.show()
#mjj
plt.hist(S_mjj,bins=100,range=(200,12000),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_mjj,bins=100,range=(200,12000),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("mjj")
plt.show()
#j1PT
plt.hist(S_j1PT,bins=100,range=(200,2000),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_j1PT,bins=100,range=(200,2000),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("j1PT")
plt.show()
#j1Eta
plt.hist(S_j1Eta,bins=100,range=(-10,10),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_j1Eta,bins=100,range=(-10,10),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("j1Eta")
plt.show()
#j1Phi
plt.hist(S_j1Phi,bins=100,range=(-10,10),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_j1Phi,bins=100,range=(-10,10),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("j1Phi")
plt.show()
#j2PT
plt.hist(S_j2PT,bins=100,range=(200,2000),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_j2PT,bins=100,range=(200,2000),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("j2PT")
plt.show()
#j2Eta
plt.hist(S_j2Eta,bins=100,range=(-10,10),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_j2Eta,bins=100,range=(-10,10),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("j2Eta")
plt.show()
#j2Phi
plt.hist(S_j2Phi,bins=100,range=(-10,10),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_j2Phi,bins=100,range=(-10,10),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("j2Phi")
plt.show()
#Etachange
plt.hist(S_Etachange,bins=100,range=(-10,10),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_Etachange,bins=100,range=(-10,10),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("Etachange")
plt.show()


# In[ ]:




