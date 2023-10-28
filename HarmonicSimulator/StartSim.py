import subprocess as sub
import numpy as np
import matplotlib.pyplot as plt




def startSimulator():
    sub.run("cd .. && cmake -B build", shell=True)
    sub.run("cd .. && cd build && make", shell=True)
    sub.run([".././build/Simulator", ".././Configs/Cfg.json"])

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

def showPhaseDiagramm(Trajectory, GraphicName):
    NameX = GraphicName + ' X'
    NameY = GraphicName + ' U'
    plt.plot(Trajectory['X'], Trajectory['U'], label = NameX)

def showDeltaEnergy(Energy, Energy1, GraphicName):
    Name = GraphicName
    plt.plot(Energy['T'], Energy['E'] - Energy1['E'], label = Name)

def showDeltaTrajectory(Trajectory, Trajectory1):
    plt.plot(Trajectory['T'], Trajectory['X'] - Trajectory1['X'], label = 'Delta X')
    plt.plot(Trajectory['T'], Trajectory['U'] - Trajectory1['U'], label = 'Delta U')

def showEnergy(Energy, GraphicName):
    Name = GraphicName + ' E'
    plt.plot(Energy['T'], Energy['E'], label = Name)


def main():
    startSimulator()
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
