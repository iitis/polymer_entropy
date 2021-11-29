"""
Reading module that prints some basic statistical information about analysed dataset.
Plots energy vs time in .png in /plots directory.
As first data points are consideered artifacts they are dropped. If more than one is dropped warning is displayed.
"""
import argparse
from experiment import Experiment

plotDirectory = 'plots'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reads tabular data for molecular dynamics")
    parser.add_argument('dirpath', type=str, help='path to a directory with data to use')
    parser.add_argument('--dropRange', type=float, default=-100.,
        help='drop first observation if their energy is lower than mean minus this many times stddev ')
    args = parser.parse_args()

    myExperiment = Experiment(args.dirpath)
    myExperiment.dropFirstObservations()
    #myExperiment.plotColumns(7)
    myExperiment.plotHistogram2D(8,9)
    print(myExperiment.get_entropy(8,9))
    print(myExperiment.get_entropy(9,10))
