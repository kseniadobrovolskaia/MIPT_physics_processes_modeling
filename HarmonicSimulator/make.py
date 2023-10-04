import subprocess as sub
import numpy as np
import matplotlib.pyplot as plt




def startSimulator():
    sub.run("cd .. && cmake -B build", shell=True)
    sub.run("cd .. && cd build && make", shell=True)
    sub.run(".././build/Simulator", shell=True)

def getTrajectory(FileNameAnalytical):
    TrajectoryTypes = np.dtype([('T', np.double), ('X', np.double), ('U', np.double)])
    AnalyticalTrajectory = np.fromfile(FileNameAnalytical, dtype=TrajectoryTypes);
    return AnalyticalTrajectory

def createGraph():
    plt.figure(100)
    plt.title("График зависимости координат X и U маятника от времени T")
    plt.xlabel('T, c')
    plt.grid()

def showTrajectory(Trajectory, GraphicName):
    NameX = GraphicName + ' X'
    NameY = GraphicName + ' U'
    plt.plot(Trajectory['T'], Trajectory['X'], label = NameX)
    plt.plot(Trajectory['T'], Trajectory['U'], label = NameY)


def main():
    startSimulator()
    FileGraphicsName = 'Graphics.png'
    FileNameAnalytic = 'Analitic.bin'
    FileNameEiler = 'Eiler.bin'
    AnalyticalTrajectory = getTrajectory(FileNameAnalytic)
    EilerTrajectory = getTrajectory(FileNameEiler)

    createGraph()
    showTrajectory(AnalyticalTrajectory, 'Analitic')
    showTrajectory(EilerTrajectory, 'Eiler')

    plt.legend()
    plt.savefig(FileGraphicsName)
    CommandToShow = 'eog ' + FileGraphicsName
    sub.run(CommandToShow, shell=True)

if __name__ == '__main__':
    main()
