"""
Reading module that prints some basic statistical information about analysed
dataset.
Plots energy vs time in .png in /plots directory.
As first data points are considered artifacts they are dropped.
If more than one is dropped warning is displayed.
"""
import argparse
import os
import itertools

from experiment import Experiment, SetOfExperiments

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
    parser.add_argument(
        '--ions', nargs='+', type=str, help="Ions to be considered", default=["Ca", "Mg", "Na"]
    )
    parser.add_argument(
        '--modes', nargs='+', type=str, help="'analysis' and 'sidechain' allowed'", default=["analysis", "sidechain"]
    )
    parser.add_argument(
        "--complex", type=str, help="Complex chosen", default="Albumin+HA"
    )
    parser.add_argument(
        "--bins", nargs='+', type=int, help="Number of bins to use for 2d histograms", default=[40,100,150]
    )
    parser.add_argument(
        "--realisations", type=int, help="Process this many first realisations", default=6
    )

    args = parser.parse_args()

    print(f"{args.complex=}")
    print(f"{args.ions=}")
    print(f"{args.modes=}")
    print(f"{args.bins=}")
    print(f"{args.realisations=}")

    for myMode in args.modes:
        for ion in args.ions:
            for i in range(1,args.realisations+1):
                file_path = os.path.join(args.datafolder, f"{args.complex}_{i}_{myMode}_{ion}.tab")

                myExperiment = Experiment(file_path)
                myExperiment.drop_first_observations()
                myExperiment.plot_columns(8, os.path.join(args.plotdir, f"realisation{i}_{ion}_"))
                angles = myExperiment.angles
                for bincount in args.bins:
                    myExperiment.plot_histogram_2d(f"{angles[0]} mers 1, 2", f"{angles[1]} mers 1, 2", os.path.join(args.plotdir, f"realisation{i}_{ion}_"), bincount)

    for myMode in args.modes:
        for ion in args.ions:
            mySetOfExperiments = SetOfExperiments(args.datafolder, args.complex, ion, myMode)

            angles = mySetOfExperiments.experiments[0].angles

            for angle1,angle2 in itertools.combinations(angles,2):
                mySetOfExperiments.entropy_distribution_percentiles(angle1, angle2, args.plotdir)
                mySetOfExperiments.entropy_distribution_realisations(angle1, angle2, args.plotdir)
                mySetOfExperiments.hist_of_entropy(f"{angle1} mers 1, 2", f"{angle2} mers 1, 2", args.plotdir)
