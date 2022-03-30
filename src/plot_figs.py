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
    parser.add_argument(
        '--ions', nargs='+', type=str, help="Ions to be considered", default=["Ca", "Mg", "Na"]
    )
    parser.add_argument(
        '--modes', nargs='+', type=str, help="'analysis' and 'sidechain' allowed'", default=["analysis", "sidechain"]
    )
    parser.add_argument(
        "--complex", type=str, help="Complex chosen", default="Albumin+HA"
    )

    args = parser.parse_args()

    Point = namedtuple('Point', ['x', 'y'])
 
    startingPoints = {
        'analysis' : [Point(8, 9), Point(10, 11), Point(8, 10)],
        'sidechain' : [Point(8, 31), Point(8, 54), Point(31, 54)]
    }

    numRealisations = 2 
    
    print(f"{args.complex=}")
    print(f"{args.ions=}")
    print(f"{args.modes=}")

    for myMode in args.modes:
        for ion in args.ions:
            for i in range(1,numRealisations+1):
                file_path = os.path.join(args.datafolder, f"{args.complex}_{i}_{myMode}_{ion}.tab")

                myExperiment = Experiment(file_path)
                myExperiment.drop_first_observations()
                myExperiment.plot_columns(8, os.path.join(args.plotdir, f"realisation{i}_{ion}_"))
                myExperiment.plot_histogram_2d(startingPoints[myMode][0].x, startingPoints[myMode][0].y, os.path.join(args.plotdir, f"realisation{i}_{ion}_"))

    for myMode in args.modes:
        for ion in args.ions:
            mySetOfExperiments = SetOfExperiments(args.datafolder, args.complex, ion, myMode)

            for p in startingPoints[myMode]:
                mySetOfExperiments.entropy_distribution_percentiles(p.x, p.y, args.plotdir)
                mySetOfExperiments.entropy_distribution_realisations(p.x, p.y, args.plotdir)
                mySetOfExperiments.hist_of_entropy(p.x, p.y, args.plotdir)
