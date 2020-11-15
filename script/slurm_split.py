#!python

import os
import sys

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: {} [nodes] [tasks per node] [total] (command)".format(sys.argv[0]))
        sys.exit(1)

    nodes = int(sys.argv[1])
    tasks_per_node = int(sys.argv[2])
    totalParam = int(sys.argv[3])

    if nodes <= 0 or tasks_per_node <= 0 or totalParam <= 0:
        print("Incorrect arguments")
        sys.exit(1)

    total_tasks = nodes * tasks_per_node
    if total_tasks > totalParam:
        print("To many tasks for the job")
        sys.exit(1)

    command = ""
    for arg in sys.argv[4:]:
        arg = arg.replace(" ", "\\ ")
        command += arg + " "

    simulations_per_task = totalParam/total_tasks
    fromParam = 0
    fullCommand = ""
    for node in range(nodes):
        srunCommand = "srun -N1 --ntasks={} bash -c '\n".format(tasks_per_node)
        for task in range(tasks_per_node):
            toParam = min(fromParam + simulations_per_task, totalParam)
            srunCommand += "    {} -P from={} -P to={} &\n".format(command, int(round(fromParam)), int(round(toParam)))
            fromParam += simulations_per_task
        srunCommand += "    wait\n' &\n\n"
        fullCommand += srunCommand
    fullCommand += "wait\n"
    os.system(fullCommand)

