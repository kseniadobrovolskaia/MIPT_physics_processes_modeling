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

def createGraph():
    plt.figure(figsize = (10, 10))
    plt.title("График зависимости координат X и скорости V от времени Т")
    plt.xlabel('T, c')
    plt.grid()

def showTrajectory(Trajectory, GraphicName):
    NameX = GraphicName + ' X'
    NameY = GraphicName + ' U'
    plt.plot(Trajectory['T'], Trajectory['X'], label = NameX)
    plt.plot(Trajectory['T'], Trajectory['U'], label = NameY)

def showEnergy(Energy, GraphicName):
    Name = GraphicName + ' E'
    plt.plot(Energy['T'], Energy['E'], label = Name)


def main():
    startSimulator()
    FileGraphicsName = 'Graphic.png'
    FileNameAnalytic = 'Analitic.bin'
    FileNameEiler    = 'Eiler.bin'
    FileNameHeun     = 'Heun.bin'
    AnalyticalTrajectory = getTrajectory(FileNameAnalytic)
    EilerTrajectory      = getTrajectory(FileNameEiler)
    HeunTrajectory       = getTrajectory(FileNameHeun)

    createGraph()
    showTrajectory(AnalyticalTrajectory, 'Analitic')
    showTrajectory(EilerTrajectory, 'Eiler')
    showTrajectory(HeunTrajectory, 'Heun')

    plt.legend()
    plt.savefig(FileGraphicsName)
    CommandToShow = 'eog ' + FileGraphicsName
    sub.run(CommandToShow, shell=True)

if __name__ == '__main__':
    main()
