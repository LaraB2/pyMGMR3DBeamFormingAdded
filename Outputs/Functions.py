#from pathlib import Path
import pandas as pd
import numpy as np
import h5py
import sys
import re
import matplotlib.pyplot as plt
import os
from contextlib import chdir
import scipy.fft as fft
import scipy.optimize as solve
import scipy.constants
#import NuRadioMC
#import NuRadioReco
import radiotools.atmosphere.models as refractiv
refractiv.default_model = 1
"""
Need to keep atm_model = 1 for compatibility with Kernelgenerator. 
"""
import math
import time
import tables


def LambaXoToRL(Lambda,X0,Xmax):
    X0Accent = abs(X0-Xmax)
    R = np.sqrt(Lambda/X0Accent)
    L = np.sqrt(X0Accent*Lambda)
    return R,L

    
def RLToLambaX0(R,L,Xmax):
    Lamba = R*L
    X0accent = np.power(L/(np.sqrt(R*L)),2)
    X0 = Xmax-X0accent
    return Lamba,X0

    
def GetcoordinatesfromVitalMGMR(hdf5name,accuracy = 5):
    xcoords = np.array([])
    ycoords = np.array([])
    InputFile = h5py.File(hdf5name)
    Observers = InputFile["observers"]
    for i in Observers.keys():
        xcoord = Observers[i].attrs['position'][0]
        ycoord = Observers[i].attrs['position'][1]
        xcoords = np.append(xcoords,xcoord)
        ycoords = np.append(ycoords,ycoord)
    distances = np.round(np.sqrt(np.add(np.power(xcoords,2),np.power(ycoords,2))),decimals = accuracy)
    angles = np.round(np.arcsin(ycoords/distances)*((180/np.pi)),decimals = accuracy)
    InputFile.close()
    return xcoords,ycoords,distances,angles

    
def GetUniqueValuesWithWeights(distances): 
    """
    distances is numpy array.
    """
    UniqueDistances = np.unique(distances)
    counts = np.array([])
    for i in UniqueDistances:
        count = np.count_nonzero(distances-i==0) 
        """
        Generates array where all the values = i are True hence 1
        All the values not equal to i are False hence 0
        """
        counts = np.append(counts,count)
    return UniqueDistances,counts

    
def GetDataFromKernelFile(PathToKernel):
    Datapd = pd.read_csv(PathToKernel,sep="\s+",skiprows = [0,1,2], header = None)
    Headerpd = pd.read_csv(PathToKernel,sep="\s+",nrows = 1, header = None)
    Datanp = np.array(Datapd)
    xmin = Headerpd[6][0]
    xmax = Headerpd[8][0]
    ymin = Headerpd[10][0]*1000
    ymax = Headerpd[12][0]*1000
    nx = Headerpd[2][0]
    ny = Headerpd[4][0]
    return xmin,xmax,ymin,ymax,nx,ny,Datanp

    
def PlotKernel(hdf5FileName,SpecificKernelName,xlabel,ylabel,title,valuemin = 0, valuemax = 0,multip = 1):
    InputFile = h5py.File(hdf5FileName)
    Dataset = InputFile[SpecificKernelName]
    plt.rcParams["figure.figsize"] = (8,8)
    Datanp = Dataset
    xmin = Dataset.attrs["tmin"]
    xmax = Dataset.attrs["tmax"]
    ymin = Dataset.attrs["DcMin"]
    ymax = Dataset.attrs["DcMax"]
    nx = int(np.round(Dataset.attrs["nt"]))
    ny = int(np.round(Dataset.attrs["nDc"]))
    xticks = np.round(np.linspace(xmin,xmax,nx),1)
    xlocs = np.arange(0,nx)
    ylocs = np.arange(0,ny)
    yticks = np.round(np.linspace(ymin,ymax,ny),1)
    if (valuemin != 0 and valuemax != 0):
        for i in range(len(Datanp)):
            for j in range(len(Datanp[i])):
                if Datanp[i][j] < valuemin:
                    Datanp[i][j] = valuemin
                if Datanp[i][j] > valuemax:
                    Datanp[i][j] = valuemax
        plt.contourf(Datanp,levels = 200)
    else:
        plt.contourf(Datanp,levels = 200)
    InputFile.close()
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(xlocs[::20],xticks[::20])
    plt.yticks(ylocs[::50],yticks[::50])
    plt.title(title,y= 1.05)
    plt.show()

    
def CalculatePKR(Kernel,LargestCalcDistance):
    hdf5Kernel = h5py.File(Kernel)
    SpecificKernel = hdf5Kernel[str(LargestCalcDistance)]
    tmin = SpecificKernel.attrs["tmin"]
    tmax = SpecificKernel.attrs["tmax"]
    ymin = SpecificKernel.attrs["DcMin"]
    ymax = SpecificKernel.attrs["DcMax"]
    ny = int(np.round(SpecificKernel.attrs["nDc"]))
    nt = int(np.round(SpecificKernel.attrs["nt"]))
    Dn = (ymax-ymin)/(ny-1)
    Nyspecific = int(np.round((LargestCalcDistance-ymin)/Dn))
    PKRNotRenom = SpecificKernel[:][Nyspecific]
    TotalSum = np.sqrt(np.sum(np.power(PKRNotRenom,2))*((np.abs(tmax-tmin)*1e-9)/nt))
    PKRRenom = PKRNotRenom/TotalSum
    return PKRRenom

    
def RunMGMR(NameMGMRoutput, extraFlags = ""):
    """
    NameMGMRoutput string of the name of the hdf5 file output from MGMR. Will replace a similarly named file
    """
    CurrentDirectory = os.getcwd()
    with chdir(CurrentDirectory):
        os.chdir("../run")
        currentdir = str(os.getcwd())
        os.system("./pyMGMR3D.sh "+currentdir+"/SIM000001.in " + str(extraFlags))
        os.replace("SIM000001.hdf5","../Outputs/"+NameMGMRoutput+".hdf5")

        
def RunBeamForming(Distances,Counts,zeta_K,Zenithangle_deg,BeamFormOutput):
    """
    Generates the beamforming Kernel to be used later for reconstruction.
    Distances is a list of all the distances of each antenna to the core. 
    Counts is the weight of each Antenna. 
    """
    CurrentDirectory = os.getcwd()
    with chdir(CurrentDirectory):
        os.chdir("../Beamforming/BeamformingRuns")
        N_antennas = len(Distances)
        with open(os.getcwd()+"/Kernel-SKA.in","w") as file: #generate the .in file to work wit fortran
            file.write(" &BeamPars\n")
            file.write("Base =   'SKx/Coarse'   ! 48 -- 300 MHz\n")
            file.write("DataBase = 'SKA/'\n")
            file.write("AntennaLayout = .true., DistMax = 500.\n")
            file.write("dN_tr=5\n")
            file.write("N_phir =  5 ! 9 ! 11  !no effect\n")
            file.write("d_h= 0.2  ! no real effect\n")
            file.write("N_h= 100 ! \n")
            file.write("ZenithAngle_deg = " + str(Zenithangle_deg)+ "\n")
            AllzetaKs = "zeta_K= "+ str(zeta_K[0])
            for ParticularHeight in zeta_K[1:]:
                AllzetaKs = AllzetaKs + ", " + str(ParticularHeight)
            file.write(AllzetaKs + "\n")
            file.write(" &end \n")
            file.write("   --- ---------- -------- ------------- -   nohup ./Beam-SKA.sh  >Beam-SKA.log 2>&1  & ")
        with open(os.getcwd()+"/AntennaPosInput.in","w") as file:
            file.write(str(N_antennas)+"\n")
            for i in range(len(Distances)):
                file.write(str(Distances[i])+","+str(Counts[i])+"\n")
        os.system("./RunKernelOnly.sh")
        os.chdir("SKx")
        OutputFile = h5py.File(str(BeamFormOutput)+".hdf5","w")
        OutputFile.attrs["focalheights"] = Distances
        for OutputHeight in zeta_K:
            if OutputHeight > 0: #the fortran code ignores negative values, so we need to ignore those here too. 
                Frontpart = int(Zenithangle_deg*np.power(10,9)+OutputHeight*np.power(10,6))
                if Zenithangle_deg <10:
                    FileName = "Coarsex-SrcKrnl-0"+str(Frontpart)+".z"
                else:
                    FileName = "Coarsex-SrcKrnl-"+str(Frontpart)+".z"
                TempData = GetDataFromKernelFile(FileName)
                dataset = OutputFile.create_dataset(str(int(OutputHeight*1000)),data= TempData[6])
                Attributestemp = ["tmin","tmax","DcMin","DcMax","nt","nDc"]
                for attriter in range(len(Attributestemp)):
                    dataset.attrs[Attributestemp[attriter]] = float(TempData[attriter])
                dataset.attrs["Axes"] = ["Dc","time"]
        OutputFile.close()
        os.replace(str(BeamFormOutput)+".hdf5","../../../Outputs/"+str(BeamFormOutput)+".hdf5")
        return None

        
def fractional_time_shift(signals, dt, shift,weights):
    """
    Vital code to do a time shift. Dt and shift need to be both in same time units. 
    """


    was_1d = signals.ndim == 1

    if was_1d:
        signals = signals[None, :]

    n_traces, N = signals.shape
    freqs = np.fft.fftfreq(N, d=dt)
    Weightsarray = np.zeros([len(signals),len(freqs)])
    for Weightsiter in range(len(weights)):
        Weightsarray[Weightsiter,:] = weights[Weightsiter]
    S = np.fft.fft(signals, axis=1)*Weightsarray
    phase = np.exp(1j * 2 * np.pi * freqs[None, :] * shift[..., None])
    shifted = np.fft.ifft(S * phase, axis=1).real

    return shifted[0] if was_1d else shifted

def BlockFreqFilter(signals,dt,Freqmin,Freqmax):
    """
    Frequency filtering of real signal.
    """
    was_1d = signals.ndim == 1

    if was_1d:
        signals = signals[None, :]
    n_traces, N = signals.shape
    freqs = np.fft.rfftfreq(N, d=dt)
    S = np.fft.rfft(signals,axis = 1)
    Smask = (freqs > Freqmin) & (freqs < Freqmax)
    filterfreqs = Smask*S
    EndSignals = np.fft.irfft(filterfreqs,axis = 1)
    return EndSignals[0] if was_1d else shifted
    
def Databeaming(MGMRhdf5,Beamformedhdf5,distances,freqmin = 50*10e6,freqmax = 300*10e6,ant_toffset = 0,Sample_Offset = 0,ConversionToSeconds = 1e-6):
    """
    the output will be in the order of the list of distances
    """
    MGMRData = h5py.File(MGMRhdf5)
    Beamformedhdf5 = h5py.File(Beamformedhdf5,"w")
    AtmosphereModel = refractiv.Atmosphere(model = 1, curved = False)
    Nantennas = len(MGMRData["observers"])
    Antennapos = np.zeros([int(Nantennas),3])
    amnttimesamples = 0
    timemin = 0
    timemax = 0
    timestep = 0
    Nyquistfreq = 0
    ZenithShowerAngle = np.radians(MGMRData["inputs"].attrs["Zen_sh"])
    for Antenna in MGMRData["observers"]:
        MGMRAntennafile = MGMRData["observers/"+Antenna]
        amnttimesamples = MGMRAntennafile.shape[0]
        timemin = MGMRAntennafile[0,0]*ConversionToSeconds
        timemax = MGMRAntennafile[amnttimesamples-1,0]*ConversionToSeconds
        timestep = (timemax-timemin)/(amnttimesamples-1) #in seconds.
        Nyquistfreq = 2/timestep
        break
        """
        Takes the first antenna and sets the lowest time, highest time and time step, all in seconds. Need seconds for the function fractional_time_shift which works in seconds. 
        """
    InitialAntennaData = np.zeros([int(Nantennas),int(amnttimesamples)],dtype="complex") #has form[Antenna, timebin]
    IterationStep = 0
    """
    To keep track of how often the Antenna loop is done
    """
    Keepingtrackcounter = 0
    for Antenna in MGMRData["observers"]:
        MGMRAntennafile = MGMRData["observers/"+Antenna]
        Antennapos[IterationStep] = [MGMRAntennafile.attrs["position"][0],MGMRAntennafile.attrs["position"][1],0]
        TempAntennaData = e_geo(np.transpose([MGMRAntennafile[:,2],MGMRAntennafile[:,1]]),MGMRAntennafile.attrs["position"][1],MGMRAntennafile.attrs["position"][0])*2*np.sqrt(np.power(Antennapos[IterationStep][0],2)+np.power(Antennapos[IterationStep][1],2))
        InitialAntennaData[IterationStep] =TempAntennaData
        IterationStep += 1
    IterationStep = 0
    ResultMatrix = np.zeros([len(distances),int(amnttimesamples)]) #has form [Focalpoint,timedata]
    for Focalpoint in range(len(distances)):
        Focalposition = np.zeros([int(Nantennas),3])
        Focalposition[:,2] = distances[Focalpoint]
        EffRefrac = AtmosphereModel.get_effective_refractivity(ZenithShowerAngle, distances[Focalpoint], 0)[0]
        DistanceFocalPointToAntenna = np.linalg.norm(Focalposition-Antennapos,axis = 1)*(EffRefrac+1)
        """
        Height above ground = 7 voor nederland, maar niet zeker of Beamforming code daarmee werkt of hoogte = 0. 
        """
        timedifference = (((DistanceFocalPointToAntenna - Focalposition[:,2]))/scipy.constants.c)
        ResultMatrix[Focalpoint] = fractional_time_shift(InitialAntennaData,timestep,timedifference,Weights).sum(axis = 0)
    """
    Above code generates the distances from focal point to antenna, then calculates the time difference in nanoseconds, and then take that timemax and time min 
    """
    for distance in range(len(distances)):
        dataset = Beamformedhdf5.create_dataset(str(distances[distance]),data = ResultMatrix[distance,:])
        dataset.attrs["distance"] = distances[distance]
    MGMRData.close()
    Beamformedhdf5.close()
    return None

def GaisserHillasCurrentProfile(X,Imax,R,L,Xmax):
    """
    Gaisser-Hillas current profile, X is numpy array, Imax,R,L,Xmax are floats. 
    1/ >= R/L*(Xmax-X) otherwise this formula will give an error. 
    """
    inbetweenstep = (1-(R/L)*(Xmax-X))
    Power = 1/np.power(R,2)
    count = 0
    for item in inbetweenstep:
        if item <0:
            inbetweenstep[count] = 0
        count += 1
    return Imax*np.power(inbetweenstep,Power)*np.exp((Xmax-X)/(L*R))

    
def CalculateWeightingFunction(kernelhdf5,LargestCalcDistance,WeightingFunctionhdf5):
    CurrentDirectory = os.getcwd()
    Kernel = h5py.File(kernelhdf5)
    hdf5WeightingFile = h5py.File(WeightingFunctionhdf5,"w")
    Attributestemp = ["DcMin","DcMax","nDc"]
    PKR = CalculatePKR(Kernel,LargestCalcDistance)
    DbCount = len(Kernel)
    for Db in Kernel:
        tempdata = np.dot(Kernel[Db],PKR)
        hdf5WeightingFile.create_dataset(Db,data = tempdata)
        for attribute in Attributestemp:
            hdf5WeightingFile[Db].attrs[attribute] = Kernel[Db].attrs[attribute]
    hdf5WeightingFile.close()

    
def MinimiseChisquare(MGMRhdf5,Kernelhdf5,LargestCalcDistance,WeightingFunctionhdf5,DataBeamedDatahdf5,ConversionToNanoS = 1000):
    """
    ConversionToNanoS is the conversion factor to get the units of time from the Datafile to Nanoseconds, which is used for the calculations in this function. 
    """
    CurrentDirectory = os.getcwd()
    with chdir(CurrentDirectory):
        os.chdir("../../Outputs")
        WeightingFunction = h5py.File(WeightingFunctionhdf5)
        Kernel = h5py.File(Kernelhdf5)
        DatabeamedData = h5py.File(DataBeamedDatahdf5)
        MGMRData = h5py.File(MGMRhdf5)
    AtmosphereModel = refractiv.Atmosphere(model = 1, curved = False)
    ZenithShowerAngle = np.radians(MGMRData["inputs"].attrs["Zen_sh"])
    MGMRtmin = 0
    MGMRtmax = 0
    MGMRnt = 0
    MGMRTimes = np.array([])
    for antenna in MGMRData["observers"]:
        MGMRAntennafile = MGMRData["observers/"+antenna]
        MGMRtmin = MGMRAntennafile[0,0]*ConversionToNanoS
        MGMRtmax = MGMRAntennafile[-1,0]*ConversionToNanoS
        MGMRnt = int(len(MGMRAntennafile[:,0]))
        MGMRTimes = MGMRAntennafile[:,0]*ConversionToNanoS
        break
    PKR = CalculatePKR(Kernel,LargestCalcDistance)
    distances = np.array([],dtype = "int")
    for item in Kernel:
        distances = np.append(distances,int(item))
    distances = np.sort(distances) #makes sure distances are sort from small to large
    Kerneltmin = Kernel[str(distances[0])].attrs["tmin"]
    Kerneltmax = Kernel[str(distances[0])].attrs["tmax"]
    Kernelnt = int(np.round(Kernel[str(distances[0])].attrs["nt"]))
    Kerneltimes = np.linspace(Kerneltmin,Kerneltmax,Kernelnt)
    """
    Assumes tmin,tmax,and nt is same for all datasets in Kernel. 
    """
    DatabeamDatainterpolFull = np.zeros([2,len(distances),len(Kerneltimes)]) # [x/y,Db,t]
    Px = np.zeros([2,len(distances)]) #[x/y,Db]
    FullDcsatmos = {}
    for distanceiter in range(len(distances)):
        DatabeamDatainterpolFull[0,distanceiter] = np.interp(Kerneltimes,MGMRTimes,DatabeamedData[str(distances[distanceiter])][0],left = 0, right = 0)
        DatabeamDatainterpolFull[1,distanceiter] = np.interp(Kerneltimes,MGMRTimes,DatabeamedData[str(distances[distanceiter])][1],left = 0, right = 0)
        Dcs = np.linspace(WeightingFunction[str(distances[distanceiter])].attrs["DcMin"],WeightingFunction[str(distances[distanceiter])].attrs["DcMax"],int(WeightingFunction[str(distances[distanceiter])].attrs["nDc"]))
        Dcsatmos = np.zeros(len(Dcs))
        for iteration in range(len(Dcs)):
            DcAtmos = AtmosphereModel.get_atmosphere(ZenithShowerAngle,h_low = Dcs[iteration])
            Dcsatmos[iteration] = DcAtmos
        FullDcsatmos[str(int(distances[distanceiter]))] = Dcsatmos
        Px[0,distanceiter] = np.dot(DatabeamDatainterpolFull[1,distanceiter],PKR)
        Px[1,distanceiter] = np.dot(DatabeamDatainterpolFull[1,distanceiter],PKR)
    def CalculateChiSquare(Fitparams):
        """
        Fitparams [Imax,R,L,Xmax]
        """
        Variation = 0
        for distanceiter in np.array([2,3,4,5,6]):
            DistanceStrForm = str(int(distances[distanceiter]))
            Guess = np.dot(WeightingFunction[DistanceStrForm],GaisserHillasCurrentProfile(FullDcsatmos[DistanceStrForm],Fitparams[0],Fitparams[1],Fitparams[2],Fitparams[3]))
            Variation += np.sum(np.power(Guess-Px[0,distanceiter],2))
        return Variation
    result = solve.minimize(CalculateChiSquare,np.array([500,0.25,250,680]),bounds = np.array([[None,None],[0.01,0.6],[100,500],[400,700]]))
    """
    General value for R is between 0.25 and 0.50, general value for L is between 150 and 400. Middle of these are good to start out with. 
    """
    return result

def MGMRBeamforminginput(MGMRINPUTLOC):
    RangeList = np.arange(1,76)+100000
    #MGMRHDF5 = h5py.File(NameMGMROUTPUT,"w")
    AntennaData = np.zeros([75,1094])
    for iteration in range(75):
        AntennaFileName = "th_"+str(RangeList[iteration])[1:]+".csv"
        AntennaFilelist = str(MGMRINPUTLOC)+ "/" + AntennaFileName
        Rawdata = pd.read_csv(AntennaFilelist, sep="\s+", engine='python', comment='!', header=None)
        RawdataArray = np.array(Rawdata)
        AntennaData[iteration] = RawdataArray[:,1]
    return AntennaData

def DatabeamingBeamforming(OutputName,distances,freqmin = 48*10e6,freqmax = 300*10e6,ant_toffset = 0,Sample_Offset = 0,ConversionToSeconds = 1e-9):
    """
    the output will be in the order of the list of distances
    """
    NanoSecondsToSeconds = 1e-9
    Beamformedhdf5 = h5py.File(OutputName,"w")
    AtmosphereModel = refractiv.Atmosphere(model = 1, curved = False)
    Nantennas = 50
    amnttimesamples = 1094
    Antennapos = np.zeros([int(Nantennas),3])
    Antennapos[:,0] = np.arange(1,51)
    RangeList = np.arange(1,len(distances)+1)+100000
    timestep = 1 # in [ns]
    ZenithShowerAngle = np.radians(43)
    InitialAntennaData = MGMRBeamforminginput("./SKA")[:50,:]
    #InitialAntennaData = np.zeros([int(Nantennas),int(amnttimesamples)],dtype="complex") #has form[Antenna, timebin]
    IterationStep = 0
    ResultMatrix = np.zeros([len(distances),int(amnttimesamples)]) #has form [Focalpoint,timedata]
    #plt.plot(InitialAntennaData[10,:])
    for Focalpoint in range(len(distances)):
        Focalposition = np.zeros([int(Nantennas),3])
        Focalposition[:,2] = distances[Focalpoint]
        EffRefrac = AtmosphereModel.get_effective_refractivity(ZenithShowerAngle, distances[Focalpoint], 0)[0]
        DistanceFocalPointToAntenna = np.linalg.norm(Focalposition-Antennapos,axis = 1)*(EffRefrac+1)
        """
        Height above ground = 7 voor nederland, maar niet zeker of Beamforming code daarmee werkt of hoogte = 0. 
        """
        WeightsArray = np.zeros(50)
        Norm = 0
        for AntennaIter in range(50):
            rc_a = 10*(AntennaIter+1)
            zeta_f = distances[Focalpoint]
            R_f = (EffRefrac+1)*np.sqrt(rc_a*rc_a+zeta_f*zeta_f)
            AntennaWeight = (10*(AntennaIter+1))
            Norm = Norm+1
            WeightsArray[AntennaIter] = 1
        WeightsArray = (WeightsArray)/Norm
        timedifference = (((DistanceFocalPointToAntenna - Focalposition[:,2]))/scipy.constants.c)*1e9 #in [ns]
        InitialAntennaData = InitialAntennaData
        ResultMatrix[Focalpoint] = (fractional_time_shift(InitialAntennaData,timestep,timedifference,WeightsArray)).sum(axis = 0)
    """
    Above code generates the distances from focal point to antenna, then calculates the time difference in nanoseconds, and then take that timemax and time min 
    """
    for distance in range(len(distances)):
        dataset = Beamformedhdf5.create_dataset(str(distances[distance]),data = ResultMatrix[distance,:])
        dataset.attrs["distance"] = distances[distance]
    MGMRData.close()
    Beamformedhdf5.close()
    return None