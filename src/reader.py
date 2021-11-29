"""
Reading module that prints some basic statistical information about analysed
dataset.
Plots energy vs time in .png in /plots directory.
As first data points are considered artifacts they are dropped.
If more than one is dropped warning is displayed.
"""
import argparse

from experiment import Experiment, SetOfExperiments

plotDirectory = "plots"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reads tabular data for molecular dynamics"
    )
    parser.add_argument(
        "dirpath", type=str, help="path to a directory with data to use"
    )
    parser.add_argument(
        "--partial_dirpath",
        type=str,
        help="partial path to a directory to use many files",
    )
    parser.add_argument(
        "--dropRange",
        type=float,
        default=-100.0,
        help="drop first observation if their energy is lower than mean minus \
        this many times stddev ",
    )
    args = parser.parse_args()

    myExperiment = Experiment(args.dirpath)
    myExperiment.dropFirstObservations()
    myExperiment.plotColumns(8)
    myExperiment.plotColumns(9)
    myExperiment.plotColumns(99)  # maximal value 99
    myExperiment.plotHistogram2D(8, 9)

    mySetOfExperiments = SetOfExperiments(args.partial_dirpath)
    mySetOfExperiments.hist_of_entropy("_sidechain_Ca.tab", 8, 9, "hist")

    mySetOfExperiments.entropy_distribution_percentiles("analysis", "Ca", 8, 9)

    mySetOfExperiments.entropy_distribution_percentiles("analysis", "Ca", 10, 11)
