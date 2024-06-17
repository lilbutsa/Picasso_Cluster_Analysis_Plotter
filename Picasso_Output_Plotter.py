# Open the cluster centers HDF5 file of interest
import h5py
import  tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
import os
import csv

pixelSizenM = 86;
milliSecondsPerFrame = 1

root = tk.Tk()
root.withdraw()

filename = filedialog.askopenfilename(filetypes=[("HDF5 Files",'*.hdf5')],title='Select the Binding Site Basic Output File')
f = h5py.File(filename, "r")
bindingEventData = list(f[list(f.keys())[0]])

filename = filedialog.askopenfilename(filetypes=[("HDF5 Files",'*.hdf5')],title='Select the DBScan Basic Output File')
f = h5py.File(filename, "r")
bindingSiteData = list(f[list(f.keys())[0]])


filename = filedialog.askopenfilename(filetypes=[("HDF5 Files",'*.hdf5')],title='Select the DBScan Cluster Centres File')
folder = os.path.dirname(filename)
f = h5py.File(filename, "r")
clusterCenterData = list(f[list(f.keys())[0]])

#%% Extract the relevent columns from the HDF5 files
allBindingEventTimes = [ milliSecondsPerFrame*x[11] for x in bindingEventData];
bindingEventSiteNumbers = [ x[14] for x in bindingEventData];

bindingEventPerSite = [ x[14] for x in bindingSiteData];
BindingSiteInputGroups = [ x[17] for x in bindingSiteData]
BindingSiteClusterGroups = [ x[18] for x in bindingSiteData]


meanBindingTimePerSite = [0]*len(bindingSiteData);
for j in range(len(BindingSiteInputGroups)):
    for i in range(len(bindingEventSiteNumbers)):
        if BindingSiteInputGroups[j]==bindingEventSiteNumbers[i]:
            meanBindingTimePerSite[j]+=allBindingEventTimes[i]
            
meanBindingTimePerSite = [ meanBindingTimePerSite[i]/bindingEventPerSite[i] for i in range(len(meanBindingTimePerSite))]


bindingSitesInClusters = [ x[14] for x in clusterCenterData];
areaOfClusters = [ pixelSizenM*pixelSizenM*x[15] for x in clusterCenterData]
densityOfClusters = [ x[14]/(pixelSizenM*pixelSizenM*x[15]) for x in clusterCenterData]
convexHullOfClusters = [ pixelSizenM*pixelSizenM*x[16] for x in clusterCenterData]
clusterGroupInputNumbers = [ x[17] for x in clusterCenterData];

meanBindingTimePerCluster = [0]*len(clusterGroupInputNumbers);
bindingEventPerCluster = [0]*len(clusterGroupInputNumbers);
countCheck = [0]*len(clusterGroupInputNumbers);
for j in range(len(clusterGroupInputNumbers)):
    for i in range(len(BindingSiteClusterGroups)):
        if clusterGroupInputNumbers[j]==BindingSiteClusterGroups[i]:
            meanBindingTimePerCluster[j]+=meanBindingTimePerSite[i]
            bindingEventPerCluster[j]+=bindingEventPerSite[i]
            countCheck[j] = countCheck[j]+1
            
meanBindingTimePerCluster = [ meanBindingTimePerCluster[i]/bindingSitesInClusters[i] for i in range(len(meanBindingTimePerCluster))]
bindingEventPerCluster = [ bindingEventPerCluster[i]/bindingSitesInClusters[i] for i in range(len(bindingEventPerCluster))]
#%% Output CSV of cluster parameters
outputFileName = os.path.join(folder,'ClusterAnalysis.csv');
file = open(outputFileName, 'w', newline='')
writer = csv.writer(file)
writer.writerow(['Binding Sites in Clusters', 'Area', 'Density', 'Mean Binding Events Per Site per Cluster', 'Mean Binding Time Per Site per Cluster'])
for i in range(len(bindingSitesInClusters)):
    writer.writerow([bindingSitesInClusters[i],areaOfClusters[i],densityOfClusters[i],bindingEventPerCluster[i],meanBindingTimePerCluster[i]])

file.close()
#%% Plot the number of binding events for each binding site
numOfBins = round((max(bindingEventPerSite)-min(bindingEventPerSite)/(2 * (np.percentile(bindingEventPerSite, 75) - np.percentile(bindingEventPerSite, 25))/ len(bindingEventPerSite)**(1/3))))
plt.figure(facecolor='white')
plt.hist(bindingEventPerSite, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,np.percentile(bindingEventPerSite, 99))
plt.xlabel("Binding Events Per Site",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'BindingEventsPerSiteHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'BindingEventsPerSiteHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the number of binding events for each binding site
numOfBins = round((max(meanBindingTimePerSite)-min(meanBindingTimePerSite)/(2 * (np.percentile(meanBindingTimePerSite, 75) - np.percentile(meanBindingTimePerSite, 25))/ len(meanBindingTimePerSite)**(1/3))))
plt.figure(facecolor='white')
plt.hist(meanBindingTimePerSite, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,np.percentile(meanBindingTimePerSite, 99))
plt.xlabel("Mean Binding Time Per Site (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'meanBindingTimePerSiteHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimePerSiteHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the mean number of binding events for each cluster
numOfBins = round((max(bindingEventPerCluster)-min(bindingEventPerCluster)/(2 * (np.percentile(bindingEventPerCluster, 75) - np.percentile(bindingEventPerCluster, 25))/ len(bindingEventPerCluster)**(1/3))))
plt.figure(facecolor='white')
plt.hist(bindingEventPerCluster, bins=numOfBins, color='#1f77b4', alpha=0.7)
#plt.xlim(0,np.percentile(bindingEventPerCluster, 99))
plt.xlabel("Mean Binding Events Per Cluster",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'BindingEventsPerClusterHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'BindingEventsPerClusterHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the binding time for each cluster
numOfBins = round((max(meanBindingTimePerCluster)-min(meanBindingTimePerCluster)/(2 * (np.percentile(meanBindingTimePerCluster, 75) - np.percentile(meanBindingTimePerCluster, 25))/ len(meanBindingTimePerCluster)**(1/3))))
plt.figure(facecolor='white')
plt.hist(meanBindingTimePerCluster, bins=numOfBins, color='#1f77b4', alpha=0.7)
#plt.xlim(0,np.percentile(meanBindingTimePerCluster, 99))
plt.xlabel("Mean Binding Time Per Cluster (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'meanBindingTimePerClusterHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimePerClusterHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the number of binding sites in each cluster
numOfBins = round((max(bindingSitesInClusters)-min(bindingSitesInClusters)/(2 * (np.percentile(bindingSitesInClusters, 75) - np.percentile(bindingSitesInClusters, 25))/ len(bindingSitesInClusters)**(1/3))))
plt.figure(facecolor='white')
plt.hist(bindingSitesInClusters, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,np.percentile(bindingSitesInClusters, 90))
plt.xlabel("Binding Sites In Clusters",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'BindingSitesInClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'BindingSitesInClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% and plot it on a log scale
logbins = np.logspace(np.log10(min(bindingSitesInClusters)),np.log10(max(bindingSitesInClusters)),30)
plt.figure(facecolor='white')
plt.hist(bindingSitesInClusters, bins=logbins, color='#1f77b4', alpha=0.7)
plt.xlim(min(bindingSitesInClusters),max(bindingSitesInClusters))
plt.xlabel("Binding Sites In Clusters",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
plt.xscale('log')
outputFileName = os.path.join(folder,'LogBindingSitesInClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'LogBindingSitesInClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot areas of clusters
numOfBins = round((max(areaOfClusters)-min(areaOfClusters))/(2 * (np.percentile(areaOfClusters, 75) - np.percentile(areaOfClusters, 25))/ len(areaOfClusters)**(1/3)))
plt.figure(facecolor='white')
plt.hist(areaOfClusters, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,np.percentile(areaOfClusters, 95))
plt.xlabel("Areas of Clusters",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'AreaOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'AreaOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% and on a log scale
logbins = np.logspace(np.log10(min(areaOfClusters)),np.log10(max(areaOfClusters)),100)
plt.figure(facecolor='white')
plt.hist(areaOfClusters, bins=logbins, color='#1f77b4', alpha=0.7)
plt.xlim(min(areaOfClusters),max(areaOfClusters))
plt.xlabel("Areas of Clusters ($nm^2$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
plt.xscale('log')
outputFileName = os.path.join(folder,'LogAreaOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'LogAreaOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% plot the density of clusters
numOfBins = round((max(densityOfClusters)-min(densityOfClusters))/(2 * (np.percentile(densityOfClusters, 75) - np.percentile(densityOfClusters, 25))/ len(densityOfClusters)**(1/3)))
plt.figure(facecolor='white')
plt.hist(densityOfClusters, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,np.percentile(densityOfClusters, 95))
plt.xlabel("Density of Clusters ($nm^{-2}$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'DensityOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'DensityOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% and on a log scale
plt.figure(facecolor='white')
plt.scatter(convexHullOfClusters,densityOfClusters, color='#1f77b4', alpha=0.7)
plt.xlabel("Area of Clusters",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Density of Clusters",fontdict = {'family':'serif','color':'black','size':18})
plt.xscale('log')
outputFileName = os.path.join(folder,'LogDensityOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'LogDensityOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
