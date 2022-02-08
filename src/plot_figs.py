"""
Reading module that prints some basic statistical information about analysed
dataset.
Plots energy vs time in .png in /plots directory.
As first data points are considered artifacts they are dropped.
If more than one is dropped warning is displayed.
"""
import argparse
import os

from collections import namedtuple
from experiment import Experiment, SetOfExperiments

plotDirectory = "plots"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reads tabular data for molecular dynamics compute entropy"
    )

    parser.add_argument(
        "--datafolder",
        type=str,
        help="path to directory for folder with realisations",
        default="data",
    )

    parser.add_argument(
        "--plotdir", type=str, help="folder plots are saved", default="pics"
    )

    args = parser.parse_args()

    Point = namedtuple('Point', ['x', 'y'])

    for i in [1, 2]:
        file_path = os.path.join(args.datafolder, f"Albumin+HA_{i}_analysis_Ca.tab")

        myExperiment = Experiment(file_path)
        myExperiment.drop_first_observations()
        myExperiment.plot_columns(8, os.path.join(args.plotdir, f"{i}Ca_"))
        # myExperiment.plotColumns(99, args.plotdir)  # maximal value 99
        myExperiment.plot_histogram_2d(8, 9, os.path.join(args.plotdir, f"{i}Ca_"))
        myExperiment.plot_histogram_2d(12, 13, os.path.join(args.plotdir, f"{i}Ca_"))
    for i in [1, 2]:
        file_path = os.path.join(args.datafolder, f"Albumin+HA_{i}_sidechain_Ca.tab")

        myExperiment = Experiment(file_path)
        myExperiment.drop_first_observations()
        myExperiment.plot_columns(8, os.path.join(args.plotdir, f"{i}Ca_"))

        myExperiment.plot_histogram_2d(8, 32, os.path.join(args.plotdir, f"{i}Ca_"))
        myExperiment.plot_histogram_2d(9, 33, os.path.join(args.plotdir, f"{i}Ca_"))

    mySetOfExperiments = SetOfExperiments(args.datafolder, "Albumin+HA", "Ca", "analysis")

    x = 8
    y = 9
    mySetOfExperiments.hist_of_entropy(x, y, args.plotdir)
    y = 31

    mySetOfExperiments = SetOfExperiments(args.datafolder, "Albumin+HA", "Ca", "sidechain")
    mySetOfExperiments.hist_of_entropy(x, y, args.plotdir)

    for ion in ["Ca", "Mg", "Na"]:

        mySetOfExperiments = SetOfExperiments(args.datafolder, "Albumin+HA", ion, "analysis")
        startingPoints = [Point(8, 9), Point(10, 11), Point(8, 10)]

        for p in startingPoints:
            mySetOfExperiments.entropy_distribution_percentiles(p.x, p.y, args.plotdir)
            mySetOfExperiments.entropy_distribution_realisations(p.x, p.y, args.plotdir)

        mySetOfExperiments = SetOfExperiments(args.datafolder, "Albumin+HA", ion, "sidechain")
        startingPoints = [Point(8, 31), Point(8, 54), Point(31, 54)]

        for p in startingPoints:
            mySetOfExperiments.entropy_distribution_percentiles(p.x, p.y, args.plotdir)
            mySetOfExperiments.entropy_distribution_realisations(p.x, p.y, args.plotdir)
