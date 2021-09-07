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


def Event_Combination(Input_Directory, Slice = False, Graph = False, optimize = False):
    """
    Takes an input directory with root files, then combines and weights them properly, arranging them in an output dictionary and applying a cut.
    Input_Directory: Defines the directory for which to take event generation files.
    Slice: Defines the slice for which to take data for. 
        Slice False does not slice the data, takes all files without divisions.
        Slice 1 Corresponds to mmjj 1000-4000
        Slice 2 Corresponds to mmjj 4000-7000
        Slice 3 Corresponds to mmjj 7000-10000
        Slice 4 Corresponds to mmjj 10000--1
    Graph: If True displays an example graph of MET average vs individual contribution
    """
    Directories = os.listdir(Input_Directory)
    TotalEvents = 0
    HistogramArray = []
    PathArray = []
    Branches = []
    CrossSections = []
    MET = []
    j1PT = []
    mjj = []
    mjj_13 = []
    mjj_23 = []
    j1Eta = []
    j1Phi = []
    j2PT = []
    j2Eta = []
    j2Phi = []
    Etachange = []
    j3PT = []
    j3Eta = []
    j3Phi = []
    weight = []
    Scales = []
    MET_Test = []
    Output = {}
    remove = []
    mjjoptimized = []
    message = "no message"
    #Ignore any background directories not of the current slice
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        #Feel Free to Change the slices.
        if Slice == 1:
            message = "Graph for Slice 1, 1000-4000"
            if "mmjj_1000_4000" in item:
                a=1
            else:
                remove.append(item)
        elif Slice == 2:
            message = "Graph for Slice 2, 4000-7000"
            if "mmjj_4000_7000" in item:
                a=1
            else:
                remove.append(item)
        elif Slice == 3:
            message = "Graph for Slice 3, 7000-10000"
            if "mmjj_7000_10000" in item:
                a=1
            else:
                remove.append(item)
        elif Slice == 4:
            message = "Graph for Slice 4, 10000--1"
            if "mmjj_10000_-1" in item:
                a=1
            else:
                remove.append(item)
        elif Slice == False:
            a=1
            message = "Graph for no slicing scheme"
    print(message)
    for item in set(remove):
        Directories.remove(item)
    remove = []
    #Test Directories to see if they actually contain a valid histogram file
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        if os.path.exists(Input_Directory+item+"/analysis/histograms.root") != True and os.path.exists(Input_Directory+item+"/analysis/SimpleAna.root") != True:
            remove.append(item)
            print("Error, Histogram not found for"+Input_Directory+item+"/analysis/histograms.root")
    #Test Directories to see if they contain a valid Cross Section Output
    for item in Directories:
        composite = ["""grep "Cross-section" """+Input_Directory+item+"/docker_mgpy.log"+"| tail -1 | awk '{print $8}'"]
        CrossSectionOutput = subprocess.run(composite, shell=True, capture_output=True)
        try: 
            z = float(CrossSectionOutput.stdout)
            if type(z) != float:
                remove.append(item)
                Statement = "File "+item+" Unable to be combined, could not find Cross Section"
                print(Statement)
                print(2)
        except:
            remove.append(item)
            Statement = "File "+item+" Unable to be combined, could not find Cross Section"
            print(Statement)
            print(3)
    

    
    for item in set(remove):
        Directories.remove(item)
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
    x=0
    
    for item in Branches:
        print(x)
        print(Directories[x])
        x+=1
        #mask = (item[b"mjj"] > 1000)&(item[b"MET"] > 200)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        mask = (item[b"mjj"] > 0)&(item[b"MET"] > 0)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
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
        Etachange = np.subtract(j1Eta,j2Eta)
        if optimize == True:
            for element in item[b"mjj_13"][mask]:
                mjj_13.append(element)
            for element in item[b"mjj_23"][mask]:
                mjj_23.append(element)
            for element in item[b"j3PT"][mask]:
                j3PT.append(element)
            for element in item[b"j3Eta"][mask]:
                j3Eta.append(element)
            for element in item[b"j3Phi"][mask]:
                j3Phi.append(element)
            for i in range(len(item[b"MET"][mask])):
                mjjarray = []
                mjjarray.append(mjj[i])
                mjjarray.append(mjj_13[i])
                mjjarray.append(mjj_23[i])
                mjjoptimized.append(np.max(mjjarray))
    #Take the weights and scale them by number of inputs
    i=0 
    for item in Branches:
        mask = (item[b"mjj"] > 0)&(item[b"MET"] > 0)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        scale = (CrossSections[i]/sum(item[b"weight"]))
        i += 1
        Scales.append(scale)
        Scaled_weight = item[b"weight"][mask]*scale/len(Directories)
        for item in Scaled_weight:
            weight.append(item)
    TotalEvents = len(MET)
    #Use this Output function to implement and variables you want to calculate from this Combination Function
    Output = {"Directories":Directories, "Number of Events":TotalEvents,
              "Cross Sections":CrossSections, "MET":MET, "j1PT":j1PT, 
              "mjj":mjj,"mjj_13":mjj_13, "mjj_23":mjj_23, "j1Eta":j1Eta,
              "j1Phi":j1Phi, "j2PT":j2PT, "j2Eta":j2Eta, "j2Phi":j2Phi,
              "j3Eta":j1Eta,"j3Phi":j1Phi, "j3PT":j3PT, "weight":weight,
              "mjjoptimized":mjjoptimized, "Etachange":Etachange}
    MET_ModifiedBranches = []
    for item in Branches:
        mask = (item[b"mjj"] > 0)&(item[b"MET"] > 0)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        masked = item[b"MET"][mask]
        MET_ModifiedBranches.append(masked)
    Weights = []
    for item in Branches:
        mask = (item[b"mjj"] > 0)&(item[b"MET"] > 0)&(item[b"njet"] >= 2)&(item[b"nElec"] == 0)&(item[b"nMuon"] == 0)
        maskedbackground = item[b"weight"][mask]
        Weights.append(maskedbackground)
    if Graph == True:
        for i in range(len(MET_ModifiedBranches)):
            plt.hist(MET_ModifiedBranches[i],bins=100,range=(0,2000),weights = Weights[i]*Scales[i], alpha=0.025, color='black')
        plt.hist(MET,bins=100,range=(0,2000),weights = weight, alpha=0.3, color='purple')
        plt.yscale('log')
        plt.title("Met Comparison, Inputs (Grey) vs Averaged Output (Purple); "+message)
        plt.show()
    return(Output)


# In[3]:


def Background(Input_Directory, Graph = False):
    """
    Takes a directory of background events, sorts them by slice, averages over each slice, then combines into one output. Applies a cut.
    Input_Directory: The directory you wish to input to average over, really only works with a background directory.
    Graph: If true displays a test graph of average vs individual contribution for each slice.
    """
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
        partition = Event_Combination(Input_Directory, Slice = item, Graph = Graph)
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


#Combine Events, Signal and Background.
SignalDirectory = "/data/users/jupyter-blonsbro/SUSY/Generations/13TeV/Signal/150mjj/"
EWKBackgroundDirectory = "/data/users/jupyter-blonsbro/SUSY/Generations/13TeV/Background/EWKBackground/"
QCDBackgroundDirectory = "/data/users/jupyter-blonsbro/SUSY/Generations/13TeV/Background/QCDBackground/"
A = Event_Combination(SignalDirectory, Graph = True)
B = Background(EWKBackgroundDirectory, Graph = True)
C = Background(QCDBackgroundDirectory, Graph = True)


# In[5]:


#Some Statistics
print("Statistics")
print("Number of Signal Events")
print(A["Number of Events"])

print("Number of Backround Events")
print(B["Number of Events"]+C["Number of Events"])


# In[6]:


print(B["Directories"])
print(C["Directories"])


# In[7]:


df1 = pd.DataFrame(A["Directories"])
df1


# In[8]:


df2 = pd.DataFrame(B["Directories"])
df2


# In[9]:


df3 = pd.DataFrame(C["Directories"])
df3


# In[10]:


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


# In[21]:


#Graph Slices as a stacked histogram
#First step below takes each slice of each background generation and averages them by slice
EWKSlice1 = Event_Combination(EWKBackgroundDirectory, 1)
EWKSlice2 = Event_Combination(EWKBackgroundDirectory, 2)
EWKSlice3 = Event_Combination(EWKBackgroundDirectory, 3)
EWKSlice4 = Event_Combination(EWKBackgroundDirectory, 4)
QCDSlice1 = Event_Combination(QCDBackgroundDirectory, 1)
QCDSlice2 = Event_Combination(QCDBackgroundDirectory, 2)
QCDSlice3 = Event_Combination(QCDBackgroundDirectory, 3)
QCDSlice4 = Event_Combination(QCDBackgroundDirectory, 4)


# In[22]:


#This part graphs the stacked histogram of the backround into one summed backround
#Darker colors means higher slice, blue is EWK, Yellow is QCD
parameter = "mjj"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(200,12000), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[23]:


parameter = "MET"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(200,1500), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[24]:


parameter = "j1PT"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(200,2000), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[25]:


parameter = "j1Eta"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(-5,5), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[26]:


parameter = "j1Phi"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(-3,3), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[27]:


parameter = "j2PT"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(200,2000), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[28]:


parameter = "j2Eta"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(-5,5), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[29]:


parameter = "j2Phi"
mjjData = [EWKSlice4[parameter], QCDSlice4[parameter], EWKSlice3[parameter], QCDSlice3[parameter], EWKSlice2[parameter], QCDSlice2[parameter], EWKSlice1[parameter], QCDSlice1[parameter]]
weightarray = [EWKSlice4["weight"], QCDSlice4["weight"], EWKSlice3["weight"], QCDSlice3["weight"], EWKSlice2["weight"], QCDSlice2["weight"], EWKSlice1["weight"], QCDSlice1["weight"]]
colors = ["navy", "darkorange", "royalblue", "goldenrod", "cornflowerblue", "gold", "lightsteelblue", "lemonchiffon"]
plt.figure()
plt.hist(mjjData, bins=100, range=(-3,3), stacked = True, weights = weightarray, color = colors )
plt.yscale('log')
plt.title("Stacked Background "+parameter)
plt.show()


# In[30]:


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
"""
#Etachange
plt.hist(S_Etachange,bins=100,range=(-10,10),weights=S_weight, alpha=0.5, color='purple')
plt.hist(B_Etachange,bins=100,range=(-10,10),weights=B_weight, alpha=0.5, color='black')
plt.yscale('log')
plt.title("Etachange")
plt.show()
"""


# In[ ]:




