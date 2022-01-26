"""
Reading module that prints some basic statistical information about analysed
dataset.
Plots energy vs time in .png in /plots directory.
As first data points are considered artifacts they are dropped.
If more than one is dropped warning is displayed.
"""
import argparse
import os

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

    for i in range(2):

        i += 1

        file_path = os.path.join(args.datafolder, f"Albumin+HA_{i}_analysis_Ca.tab")

        myExperiment = Experiment(file_path)
        myExperiment.dropFirstObservations()
        myExperiment.plotColumns(8, os.path.join(args.plotdir, f"{i}Ca_"))
        # myExperiment.plotColumns(99, args.plotdir)  # maximal value 99
        myExperiment.plotHistogram2D(8, 9, os.path.join(args.plotdir, f"{i}Ca_"))
        myExperiment.plotHistogram2D(12, 13, os.path.join(args.plotdir, f"{i}Ca_"))
    for i in range(2):

        i += 1

        file_path = os.path.join(args.datafolder, f"Albumin+HA_{i}_sidechain_Ca.tab")

        myExperiment = Experiment(file_path)
        myExperiment.dropFirstObservations()
        myExperiment.plotColumns(8, os.path.join(args.plotdir, f"{i}Ca_"))

        myExperiment.plotHistogram2D(8, 32, os.path.join(args.plotdir, f"{i}Ca_"))
        myExperiment.plotHistogram2D(9, 33, os.path.join(args.plotdir, f"{i}Ca_"))

    mySetOfExperiments = SetOfExperiments(os.path.join(args.datafolder, "Albumin+HA_"))

    x = 8
    y = 9
    mySetOfExperiments.hist_of_entropy("_analysis_Ca.tab", x, y, args.plotdir)
    y = 32
    mySetOfExperiments.hist_of_entropy("_sidechain_Ca.tab", x, y, args.plotdir)

    for ion in ["Ca", "Mg", "Na"]:

        type = "analysis"
        x_start = 8
        y_start = 9

        mySetOfExperiments.entropy_distribution_percentiles(
            type, ion, x_start, y_start, args.plotdir
        )

        mySetOfExperiments.entropy_distribution_realisations(
            type, ion, x_start, y_start, args.plotdir
        )

        x_start = 10
        y_start = 11

        mySetOfExperiments.entropy_distribution_percentiles(
            type, ion, x_start, y_start, args.plotdir
        )

        mySetOfExperiments.entropy_distribution_realisations(
            type, ion, x_start, y_start, args.plotdir
        )

        x_start = 8
        y_start = 10

        mySetOfExperiments.entropy_distribution_percentiles(
            type, ion, x_start, y_start, args.plotdir
        )

        mySetOfExperiments.entropy_distribution_realisations(
            type, ion, x_start, y_start, args.plotdir
        )

        type = "sidechain"
        x_start = 8
        y_start = 32

        mySetOfExperiments.entropy_distribution_percentiles(
            type, ion, x_start, y_start, args.plotdir
        )

        mySetOfExperiments.entropy_distribution_realisations(
            type, ion, x_start, y_start, args.plotdir
        )

        x_start = 8
        y_start = 56

        mySetOfExperiments.entropy_distribution_percentiles(
            type, ion, x_start, y_start, args.plotdir
        )

        mySetOfExperiments.entropy_distribution_realisations(
            type, ion, x_start, y_start, args.plotdir
        )

        x_start = 32
        y_start = 56

        mySetOfExperiments.entropy_distribution_percentiles(
            type, ion, x_start, y_start, args.plotdir
        )

        mySetOfExperiments.entropy_distribution_realisations(
            type, ion, x_start, y_start, args.plotdir
        )
