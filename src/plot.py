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

from experiment import ExperimentalData

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
        '--chains', nargs='+', type=str, help="'analysis' and 'sidechain' allowed'", default=["analysis"]
    )
    parser.add_argument(
        "--complex", nargs='+', type=str, help="Complex chosen", default=["Albumin+HA"]
    )
    parser.add_argument(
        "--bins", nargs='+', type=int, help="Number of bins to use for 2d histograms", default=[100]
    )

    args = parser.parse_args()

    for k, v in vars(args).items():
        print(f"{k}={v}")

    myData = ExperimentalData(args.datafolder)

    criteria = {
                 'ion': args.ions,
                 'chain': args.chains,
                 'complex' : args.complex
               }

    myData.call_method_by_criteria('plot_column', criteria, 8, args.plotdir)

    for bincount in args.bins:
        for ion in args.ions:
            for chain in args.chains:
                for complex in args.complex:
                    myCriteria = { 'ion': ion, 'chain': chain, 'complex': complex }
                    myData.call_method_by_criteria('plot_angle_histogram', myCriteria, "ϕ₁₄","ψ₁₄", args.plotdir, bincount)
                    myData.call_method_by_criteria('plot_angle_histogram', myCriteria, "ϕ₁₃","ψ₁₃", args.plotdir, bincount)

    for complex in args.complex:
        for ion in args.ions:
            for chain in args.chains:
                myCriteria = { 'ion': ion, 'chain': chain, 'complex': complex }
                myData.entropy_distribution_percentiles(myCriteria, "ϕ₁₄","ψ₁₄", args.plotdir)
                myData.entropy_distribution_percentiles(myCriteria, "ϕ₁₃","ψ₁₃", args.plotdir)
                myData.entropy_distribution_realisations(myCriteria, "ϕ₁₄","ψ₁₄", args.plotdir)
                myData.entropy_distribution_realisations(myCriteria, "ϕ₁₃","ψ₁₃", args.plotdir)
                myData.hist_of_entropy(criteria, "ϕ₁₄ mers 1, 2", "ψ₁₄ mers 1, 2", args.plotdir)
                myData.hist_of_entropy(criteria, "ϕ₁₃ mers 1, 2", "ψ₁₃ mers 1, 2", args.plotdir)

    myData.plot21(criteria, args.plotdir)
