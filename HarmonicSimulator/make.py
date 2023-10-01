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

def showTrajectory(AnalyticalTrajectory):
    plt.figure(100)
    plt.title("График зависимости координат X и U маятника от времени T")
    plt.xlabel('T, c')
    plt.plot(AnalyticalTrajectory['T'], AnalyticalTrajectory['X'], label='Analytical X')
    plt.plot(AnalyticalTrajectory['T'], AnalyticalTrajectory['U'], label='Analytical U')
    plt.legend()
    plt.grid()
    plt.savefig("Analitycal.png")
    sub.run("eog Analitycal.png", shell=True)

def main():
    startSimulator()
    FileNameAnalytical = 'Solution.bin'
    AnalyticalTrajectory = getTrajectory(FileNameAnalytical)
    showTrajectory(AnalyticalTrajectory)

if __name__ == '__main__':
    main()
