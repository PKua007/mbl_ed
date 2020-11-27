#!/usr/bin/python3

import os
import sys

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: {} [nodes] [processes per node] [cores per process] [total] (mbl command)".format(sys.argv[0]))
        sys.exit(1)

    nodes = int(sys.argv[1])
    processesPerNode = int(sys.argv[2])
    coresPerProcess = int(sys.argv[3])
    totalParam = int(sys.argv[4])

    if nodes <= 0 or processesPerNode <= 0 or totalParam <= 0 or coresPerProcess <= 0:
        print("Incorrect arguments")
        sys.exit(1)

    totalProcesses = nodes * processesPerNode
    if totalProcesses > totalParam:
        print("To many tasks for the job")
        sys.exit(1)

    mblCommand = ""
    for arg in sys.argv[5:]:
        arg = arg.replace(" ", "\\ ")
        mblCommand += arg + " "

    simulationsPerProcess = totalParam/totalProcesses
    fromParam = 0
    fullCommand = ""
    for node in range(nodes):
        srunCommand = "srun -N 1 -n 1 -c {} bash -c '\n".format(processesPerNode * coresPerProcess)
        for process in range(processesPerNode):
            toParam = min(fromParam + simulationsPerProcess, totalParam)
            srunCommand += "    {} -P from={} -P to={} &\n".format(mblCommand, int(round(fromParam)), int(round(toParam)))
            fromParam += simulationsPerProcess
        srunCommand += "    wait\n' &\n\n"
        fullCommand += srunCommand
    fullCommand += "wait\n"
    os.system(fullCommand)

