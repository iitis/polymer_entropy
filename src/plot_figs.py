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
        description="Reads tabular data for molecular dynamics compute entropy"
    )
    parser.add_argument(
        "--filepath",
        type=str,
        help="path to directory for single realisation",
        default="",
    )

    parser.add_argument(
        "--datafolder",
        type=str,
        help="path to directory for folder with realisations",
        default="../data/",
    )

    parser.add_argument(
        "--namecommon",
        type=str,
        help="common first part of files names",
        default="Albumin+HA_",
    )

    parser.add_argument(
        "--plotdir", type=str, help="folder where plots are saved", default="pics/"
    )

    args = parser.parse_args()

    if args.filepath != "":  # otherwise it performs only entropy analysis
        myExperiment = Experiment(args.filepath)
        myExperiment.dropFirstObservations()
        myExperiment.plotColumns(8, args.plotdir)
        myExperiment.plotColumns(9, args.plotdir)
        myExperiment.plotColumns(99, args.plotdir)  # maximal value 99
        myExperiment.plotHistogram2D(8, 9, args.plotdir)

    mySetOfExperiments = SetOfExperiments(args.datafolder + args.namecommon)

    x = 8
    y = 9
    mySetOfExperiments.hist_of_entropy("_analysis_Ca.tab", x, y, args.plotdir)

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
