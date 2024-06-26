# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 12:33:05 2024

@author: Shirin
"""

# Open the cluster centers HDF5 file of interest
import h5py
import  tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
import os
import csv

pixelSizenM = 86;
milliSecondsPerFrame = 150

root = tk.Tk()
root.withdraw()

filename1 = filedialog.askopenfilename(filetypes=[("HDF5 Files",'*.hdf5')],title='Select the Binding Site Basic Output File')
filename2 = filedialog.askopenfilename(filetypes=[("HDF5 Files",'*.hdf5')],title='Select the DBScan Basic Output File')
filename3 = filedialog.askopenfilename(filetypes=[("HDF5 Files",'*.hdf5')],title='Select the DBScan Cluster Centres File')

folder = os.path.dirname(filename1)

f = h5py.File(filename1, "r") # import basic output file
bindingEventData = f["locs"]

f = h5py.File(filename2, "r") # import DBscan output
bindingSiteData = f["locs"]

f = h5py.File(filename3, "r") # import DBscan clusters
clusterCenterData = f["locs"]

#%% Extract the relevent columns from the HDF5 files
allBindingEventTimes = list(bindingEventData["len"])
allBindingEventTimes = [milliSecondsPerFrame*x for x in allBindingEventTimes] # scale by frame rate
bindingEventSiteNumbers = list(bindingEventData["group"]);


bindingEventPerSite = bindingSiteData["n"];
BindingSiteInputGroups = bindingSiteData["group_input"];
BindingSiteClusterGroups = bindingSiteData["group"];

# calculate mean binding time per site
meanBindingTimePerSite = [0]*len(bindingSiteData);
groupNumberPositions = [-1]*(max(BindingSiteInputGroups)+1)
for i in range(len(BindingSiteInputGroups)):
    groupNumberPositions[BindingSiteInputGroups[i]]=i

for i in range(len(bindingEventSiteNumbers)):
    siteNum = bindingEventSiteNumbers[i];
    meanBindingTimePerSite[groupNumberPositions[siteNum]]+=allBindingEventTimes[i]    
    
meanBindingTimePerSite = [ meanBindingTimePerSite[i]/bindingEventPerSite[i] for i in range(len(meanBindingTimePerSite))]    


bindingSitesInClusters = clusterCenterData["n"]

areaOfClusters = clusterCenterData["area"]
areaOfClusters = [pixelSizenM*pixelSizenM*x for x in areaOfClusters] # scale by area of a single pixel

densityOfClusters = [bindingSitesInClusters[x]/areaOfClusters[x] for x in range(len(bindingSitesInClusters))]

convexHullOfClusters = clusterCenterData["convexhull"]
convexHullOfClusters = [pixelSizenM*pixelSizenM*x for x in convexHullOfClusters] # scale by area of a single pixel
clusterGroupInputNumbers = clusterCenterData["group"];

#calculate mean binding events per site in cluster

meanBindingTimePerCluster = [0]*len(clusterGroupInputNumbers);
bindingEventPerCluster = [0]*len(clusterGroupInputNumbers);
countCheck = [0]*len(clusterGroupInputNumbers);

groupNumberPositions = [-1]*(max(clusterGroupInputNumbers)+1)
for i in range(len(clusterGroupInputNumbers)):
    groupNumberPositions[clusterGroupInputNumbers[i]]=i
    
for i in range(len(BindingSiteClusterGroups)):
    siteNum = BindingSiteClusterGroups[i]
    meanBindingTimePerCluster[groupNumberPositions[siteNum]]+=meanBindingTimePerSite[i]
    bindingEventPerCluster[groupNumberPositions[siteNum]]+=bindingEventPerSite[i]
    countCheck[groupNumberPositions[siteNum]] = countCheck[groupNumberPositions[siteNum]]+1     
            
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
minPercentToPlot = 0
maxPercentToPlot = 99

toplot = bindingEventPerSite;
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))
if (max(toplot)-min(toplot))/numOfBins<1:
    numOfBins = (max(toplot)-min(toplot));
    
plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))
plt.xlabel("Binding Events Per Site",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'BindingEventsPerSiteHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'BindingEventsPerSiteHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the mean binding time per site
minPercentToPlot = 0
maxPercentToPlot = 99

toplot = meanBindingTimePerSite;
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = round((max(meanBindingTimePerSite)-min(meanBindingTimePerSite))/(2 * (np.percentile(meanBindingTimePerSite, 75) - np.percentile(meanBindingTimePerSite, 25))/ len(meanBindingTimePerSite)**(1/3)))

plt.figure(facecolor='white')
plt.hist(meanBindingTimePerSite, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,np.percentile(meanBindingTimePerSite, 99))
plt.xlabel("Mean Binding Time Per Site (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'meanBindingTimePerSiteHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimePerSiteHistogram.svg');
numOfBins = round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))
if numOfBins>100:
    numOfBins = 100;

plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))
plt.xlabel("Mean Binding Time Per Site (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'meanBindingTimePerSiteHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimePerSiteHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Mean binding time versus number of events as density
XMaxPercentToPlot = 95
YMaxPercentToPlot = 95
numOfBins = 30
colorMin = 0
colorMax = 40


toplot = bindingEventPerSite;
toplot2 = meanBindingTimePerSite;

plotRange = [[0 , np.percentile(toplot, XMaxPercentToPlot)],[0, np.percentile(toplot2, YMaxPercentToPlot)]]

plt.figure(facecolor='white')
plt.hist2d(toplot,toplot2, numOfBins,plotRange,False,None,colorMin,colorMax)
plt.xlabel("Binding Events Per Site",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Mean Binding Time\nPer Site (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.colorbar()
outputFileName = os.path.join(folder,'meanBindingTimeVsNumberOfEvents_DensityPlot.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimeVsNumberOfEvents_DensityPlot.svg');
plt.savefig(outputFileName)
plt.show()
#%% Mean binding time versus number of events as scatter plot
XMaxPercentToPlot = 100
YMaxPercentToPlot = 99.9

toplot = bindingEventPerSite;
toplot2 = meanBindingTimePerSite;

plotRange = [0 , np.percentile(toplot, XMaxPercentToPlot)]
plotRange2 = [0, np.percentile(toplot2, YMaxPercentToPlot)]

plt.figure(facecolor='white')
plt.scatter(toplot,toplot2, color='#1f77b4', alpha=0.7)
plt.xlabel("Binding Events Per Site",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Mean Binding Time\nPer Site (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.xlim(plotRange)
plt.ylim(plotRange2)
#plt.xscale('log')
#plt.yscale('log')
outputFileName = os.path.join(folder,'meanBindingTimeVsNumberOfEvents_Scatter.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimeVsNumberOfEvents_Scatter.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the mean number of binding events per site for each cluster
minPercentToPlot = 0
maxPercentToPlot = 99

toplot = bindingEventPerCluster;
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))


plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))
plt.xlabel("Mean Binding Events Per Cluster",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'BindingEventsPerClusterHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'BindingEventsPerClusterHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the binding time for each cluster
minPercentToPlot = 0
maxPercentToPlot = 99.9

toplot = meanBindingTimePerCluster;
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))

plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))
plt.xlabel("Mean Binding Time Per Cluster (ms)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'meanBindingTimePerClusterHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'meanBindingTimePerClusterHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% Plot the number of binding sites in each cluster
minPercentToPlot = 0
maxPercentToPlot = 99

toplot = bindingSitesInClusters;
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))

plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))
plt.xlabel("Binding Sites In Clusters",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'BindingSitesInClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'BindingSitesInClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% number of binding sites in each cluster on a log scale
numOfBins = 30;
toplot = bindingSitesInClusters;
logbins = np.logspace(np.log10(min(toplot)),np.log10(max(toplot)),numOfBins)
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
minPercentToPlot = 0
maxPercentToPlot = 90

toplot = [x/1000000 for x in areaOfClusters];
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = 2*round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))

plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))
plt.xlabel("Areas of Clusters ($μm^2$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'AreaOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'AreaOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% areas of clusters on a log scale
numOfBins = 30;
toplot = [x/1000000 for x in areaOfClusters];
logbins = np.logspace(np.log10(min(toplot)),np.log10(max(toplot)),numOfBins)
plt.figure(facecolor='white')
plt.hist(toplot, bins=logbins, color='#1f77b4', alpha=0.7)
plt.xlim(min(toplot),max(toplot))
plt.xlabel("Areas of Clusters ($μm^2$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
plt.xscale('log')
outputFileName = os.path.join(folder,'LogAreaOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'LogAreaOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% plot the density of clusters
minPercentToPlot = 0
maxPercentToPlot = 99.5

toplot = [x*1000000 for x in densityOfClusters];
plotRange = [np.percentile(toplot, minPercentToPlot) , np.percentile(toplot, maxPercentToPlot)]
toplot = [x for x in toplot if x>=plotRange[0] and x<=plotRange[1]]

numOfBins = round((max(toplot)-min(toplot))/(2 * (np.percentile(toplot, 75) - np.percentile(toplot, 25))/ len(toplot)**(1/3)))

plt.figure(facecolor='white')
plt.hist(toplot, bins=numOfBins, color='#1f77b4', alpha=0.7)
plt.xlim(0,max(toplot))

plt.xlabel("Density of Clusters ($μm^{-2}$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
outputFileName = os.path.join(folder,'DensityOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'DensityOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()
#%% density of clusters on a log scale
numOfBins = 30;

toplot = [x*1000000 for x in densityOfClusters];# delete this to plot only percentage interval

logbins = np.logspace(np.log10(min(toplot)),np.log10(max(toplot)),numOfBins)
plt.figure(facecolor='white')
plt.hist(toplot, bins=logbins, color='#1f77b4', alpha=0.7)
plt.xlim(min(toplot),max(toplot))
plt.xlabel("Density of Clusters ($μm^{-2}$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Count",fontdict = {'family':'serif','color':'black','size':18})
plt.xscale('log')
outputFileName = os.path.join(folder,'LogDensityOfClustersHistogram.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'LogDensityOfClustersHistogram.svg');
plt.savefig(outputFileName)
plt.show()

#%% Log Area versus density of clusters as a scatter plot
toplot = [x/1000000 for x in areaOfClusters];
toplot2 = [x*1000000 for x in densityOfClusters];
plt.figure(facecolor='white')
plt.scatter(toplot,toplot2, color='#1f77b4', alpha=0.7)
plt.xlabel("Areas of Clusters ($μm^2$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Density of Clusters ($μm^{-2}$)",fontdict = {'family':'serif','color':'black','size':18})
plt.xscale('log')
plt.yscale('log')
outputFileName = os.path.join(folder,'logAreaVsDensityOfClusters_Scatter.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'logAreaVsDensityOfClusters_Scatter.svg');
plt.savefig(outputFileName)
plt.show()
#%% Area versus density of clusters as a density plot
XMaxPercentToPlot = 95
YMaxPercentToPlot = 95
numOfBins = 20
colorMin = 0
colorMax = 10

toplot = [x/1000000 for x in areaOfClusters];
toplot2 = [x*1000000 for x in densityOfClusters];

plotRange = [[min(toplot) , np.percentile(toplot, XMaxPercentToPlot)],[max(toplot), np.percentile(toplot2, YMaxPercentToPlot)]]

logbins = [np.logspace(np.log10(min(toplot)),np.log10(max(toplot)),numOfBins),np.logspace(np.log10(min(toplot2)),np.log10(max(toplot2)),numOfBins)]

plt.figure(facecolor='white')
plt.hist2d(toplot,toplot2,logbins ,plotRange,False,None,colorMin,colorMax)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Areas of Clusters ($μm^2$)",fontdict = {'family':'serif','color':'black','size':18})
plt.ylabel("Density of Clusters ($μm^{-2}$)",fontdict = {'family':'serif','color':'black','size':18})
plt.colorbar()
outputFileName = os.path.join(folder,'LogAreaVsDensityOfClusters_Density.png');
plt.savefig(outputFileName, dpi=900)
outputFileName = os.path.join(folder,'LogAreaVsDensityOfClusters_Density.svg');
plt.savefig(outputFileName)
plt.show()

