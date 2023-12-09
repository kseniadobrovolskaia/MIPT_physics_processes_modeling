import subprocess as sub
import numpy as np
import math
import matplotlib.pyplot as plt
import json


def getSteadyAmplitude(Gamma, OwnOmega, Step, Trajectory):
    # The time after which the amplitude of natural oscillations 
    # will die out exp(2) times
    Tau = 2.0 / (Gamma * OwnOmega)
    Start = int(Tau / Step)
    Max = 0
    for T in range(Start, len(Trajectory['X']) - Start):
        if (Trajectory['X'][T] > Max):
            Max = Trajectory['X'][T]
    return Max

def getAnalitycalACHH(Gamma, OwnOmega, F, Points):
    A = []
    for Point in Points:
        A0 = F / math.sqrt((OwnOmega**2 - Point**2)**2 + 4 * Gamma**2 * Point**2)
        A.append(A0)
    return A
        
def buildSimulator():
    sub.run("cd .. && cmake -B build", shell=True)
    sub.run("cd .. && cd build && make", shell=True)

def startSimulator(CfgName):
    FullCfgName = ".././Configs/" + CfgName
    sub.run([".././build/Simulator", FullCfgName])

def getTrajectory(FileName):
    TrajectoryTypes = np.dtype([('T', np.double), ('X', np.double), ('U', np.double)])
    Trajectory = np.fromfile(FileName, dtype=TrajectoryTypes);
    return Trajectory

def getEnergy(FileName):
    EnergyTypes = np.dtype([('T', np.double), ('E', np.double)])
    Energy = np.fromfile(FileName, dtype=EnergyTypes);
    return Energy

def createGraph(TitleName):
    plt.figure(figsize = (10, 10))
    plt.title(TitleName)
    plt.grid()

def showTrajectory(Trajectory, GraphicName):
    NameX = GraphicName + ' X'
    NameY = GraphicName + ' U'
    plt.plot(Trajectory['T'], Trajectory['X'], label = NameX)
    plt.plot(Trajectory['T'], Trajectory['U'], label = NameY)

def showX(Trajectory, GraphicName):
    NameX = GraphicName + ' X'
    plt.plot(Trajectory['T'], Trajectory['X'], label = NameX)

def showU(Trajectory, GraphicName):
    NameU = GraphicName + ' U'
    plt.plot(Trajectory['T'], Trajectory['U'], label = NameU)

def showPhaseDiagramm(Trajectory, GraphicName):
    NameX = GraphicName + ' X'
    NameY = GraphicName + ' U'
    plt.plot(Trajectory['X'], Trajectory['U'], label = NameX)

def showDeltaEnergy(Energy, Energy1, GraphicName):
    Name = GraphicName
    plt.plot(Energy['T'], Energy['E'] - Energy1['E'], label = Name)

def showDeltaTrajectory(Trajectory, Trajectory1, GraphicName):
    plt.plot(Trajectory['T'], Trajectory['X'] - Trajectory1['X'], label = GraphicName + 'Delta X')
    plt.plot(Trajectory['T'], Trajectory['U'] - Trajectory1['U'], label = GraphicName + 'Delta U')

def showDeltaX(Trajectory, Trajectory1, GraphicName):
    plt.plot(Trajectory['T'], Trajectory['X'] - Trajectory1['X'], label = GraphicName)

def showDeltaU(Trajectory, Trajectory1, GraphicName):
    plt.plot(Trajectory['T'], Trajectory['U'] - Trajectory1['U'], label = GraphicName)

def showEnergy(Energy, GraphicName):
    Name = GraphicName + ' E'
    plt.plot(Energy['T'], Energy['E'], label = Name)
    
def writeJsonInFile(Cfg, FileName):
    with open(".././Configs/" + FileName, 'w') as File:
        json.dump(Cfg, File)


def main():
    buildSimulator()
    startSimulator("Cfg.json")
    FileNameAnalyticMath   = 'AnaliticMath.bin'
    FileNameEilerMath      = 'EilerMath.bin'
    FileNameHeunMath       = 'HeunMath.bin'
    FileNameRungeKuttaMath = 'RungeKuttaMath.bin'
    AnalyticalMathTrajectory = getTrajectory(FileNameAnalyticMath)
    EilerMathTrajectory      = getTrajectory(FileNameEilerMath)
    HeunMathTrajectory       = getTrajectory(FileNameHeunMath)
    RungeKuttaMathTrajectory = getTrajectory(FileNameRungeKuttaMath)

    FileNameEilerPhys      = 'EilerPhys.bin'
    FileNameHeunPhys       = 'HeunPhys.bin'
    FileNameRungeKuttaPhys = 'RungeKuttaPhys.bin'
    EilerPhysTrajectory      = getTrajectory(FileNameEilerPhys)
    HeunPhysTrajectory       = getTrajectory(FileNameHeunPhys)
    RungeKuttaPhysTrajectory = getTrajectory(FileNameRungeKuttaPhys)

    FileNameAnalyticMathE   = 'AnaliticMathEnergy.bin'
    FileNameEilerMathE      = 'EilerMathEnergy.bin'
    FileNameHeunMathE       = 'HeunMathEnergy.bin'
    FileNameRungeKuttaMathE = 'RungeKuttaMathEnergy.bin'
    AnalyticalMathTrajectoryE = getEnergy(FileNameAnalyticMathE)
    EilerMathTrajectoryE      = getEnergy(FileNameEilerMathE)
    HeunMathTrajectoryE       = getEnergy(FileNameHeunMathE)
    RungeKuttaMathTrajectoryE = getEnergy(FileNameRungeKuttaMathE)

    FileNameEilerPhysE      = 'EilerPhysEnergy.bin'
    FileNameHeunPhysE       = 'HeunPhysEnergy.bin'
    FileNameRungeKuttaPhysE = 'RungeKuttaPhysEnergy.bin'
    EilerPhysTrajectoryE      = getEnergy(FileNameEilerPhysE)
    HeunPhysTrajectoryE       = getEnergy(FileNameHeunPhysE)
    RungeKuttaPhysTrajectoryE = getEnergy(FileNameRungeKuttaPhysE)

    createGraph("График")
    showTrajectory(EilerMathTrajectory, 'EilerMath')
    showTrajectory(EilerPhysTrajectory, 'EilerPhys')

    plt.legend()
    plt.savefig("График.png")
    CommandToShow = 'eog ' + "График.png"
    sub.run(CommandToShow, shell=True)

if __name__ == '__main__':
    main()
