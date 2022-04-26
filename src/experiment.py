"""
A module holding Experiment class that represent entire dataset
 from a file along with basic operations on it.
"""
import os
import glob
import re
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.stats import entropy
from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable

class ColumnMeaning:
    """Describes contents of a column"""
    def __init__(self, description: str, unit: str):
        self.description = description
        self.unit = unit

    def __str__(self):
        return f"{self.description} [{self.unit}]"

class Experiment:
    """
    A class that represent an experiment of entropy calculation
    from angles from the single file
    """

    plot_dpi = 350

    angles_config = {
        ('Albumin+HA','analysis'): ["ϕ₁₄","ψ₁₄","ϕ₁₃","ψ₁₃"],
        #('Albumin+HA','side chain'): ["γ","ω","δ"],
        ('Albumin+CS6','analysis'): ["ϕ₁₄","ψ₁₄","ϕ₁₃","ψ₁₃"],
    }

    bind_energies = { #FIXME this should be handled in a more flexible way
        ('Albumin+CS6','analysis'): {
            1: 4.356,
            2: 4.354,
            3: 4.348,
            4: 4.347,
            5: 4.339,
            6: 4.322,
            7: 4.314,
            8: 4.308,
            9: 4.274,
            10: 4.249,
        },
        ('Albumin+HA','analysis'): {
            1: 5.416,
            2: 4.984,
            3: 4.945,
            4: 4.715,
            5: 4.69,
            6: 4.688,
            7: 4.687,
            8: 4.684,
            9: 4.635,
            10: 4.62,
        }
    }


    def __init__(self, filepath: str):
        split_items = re.split('_|\.',os.path.basename(filepath)) # pylint: disable=W1401
        self.complex, self.num_realisation, self.chain, self.ion, _ = split_items
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2, engine="python")
        self.angles = self.angles_config[(self.complex,self.chain)]
        self.no_mers = 24
        self.columns = [
            ColumnMeaning("Time", "ps"),                            #First column, index 0
            ColumnMeaning("Total energy of the system", "kJ/mol"),
            ColumnMeaning("Total bond energy", "kJ/mol"),
            ColumnMeaning("Total angle energy", "kJ/mol"),
            ColumnMeaning("Total dihedral energy", "kJ/mol"),
            ColumnMeaning("Total planar energy", "kJ/mol"),
            ColumnMeaning("Total Coulomb energy", "kJ/mol"),
            ColumnMeaning("Total Van der Waals energy", "kJ/mol"),  #Column index 7
        ]
        self.bind_energy = self.bind_energies[(self.complex, self.chain)][int(self.num_realisation)]

        if self.complex == 'Albumin+HA':
            if self.chain == 'analysis':
                for mer in range(self.no_mers):
                    for angle in self.angles:
                        if angle.endswith('₁₃'):
                            if mer+2 <= self.no_mers:
                                self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))
                        else:
                            self.columns.append(ColumnMeaning(f"{angle} mer {mer+1}", "deg"))
            elif self.chain == 'side chain':
                for angle in self.angles:
                    for mer in range(self.no_mers):
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))

        if self.complex == 'Albumin+CS6':
            if self.chain == 'analysis':
                for mer in range(self.no_mers):
                    for angle in self.angles:
                        if angle.endswith('₁₃'):
                            if mer+2 <= self.no_mers:
                                self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))
                        else:
                            self.columns.append(ColumnMeaning(f"{angle} mer {mer+1}", "deg"))
            elif self.chain == 'side chain':
                pass #not supported

        self.drop_first_observations() #Drop in constructor only so that we know data is stabilised

    def __str__(self):
        """
        Print experiment summary - for debug purposes
        """
        return f"{self.complex}_{self.num_realisation}_{self.chain}_{self.ion}"

    def get_colnum_by_meaning(self, meaning: str):
        """
        Returns column index when provided with column meaning
        """
        cols = [ x.description for x in self.columns]
        return cols.index(meaning)

    def drop_first_observations(self):
        """
        Drops first observations so that only points after stabilisation are
        further considered
        """
        self.dataframe.drop(index=self.dataframe.index[0], axis=0, inplace=True)

    def plot_column(self, ycol: int, plotdir: str = "plots"):
        """
        Plots column vs time. Saves to png file.
        """
        xcol = 0
        mytitle = f"{self.complex} {self.ion} R{self.num_realisation}"

        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]

        plt.plot(x, y)
        plt.xlabel(self.columns[xcol])
        plt.ylabel(self.columns[ycol])
        plt.title(mytitle)

        plot_filepath = os.path.join(plotdir,f"{self}_c{ycol}.png")
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.clf()
        plt.close()

    def plot_angle_histogram(self, angle_x, angle_y, plotdir: str = "plots", numbins = 100):
        """
        Plots histogram fo angles for all subsequent mers in experiment (realisation)
        """
        x_data = np.array([])
        y_data = np.array([])
        if angle_x.endswith('₁₃'):
            x_columns_of_interest = [ f"{angle_x} mers {mer+1}, {mer+2}" for mer in range(self.no_mers-1) ]
        else:
            x_columns_of_interest = [ f"{angle_x} mer {mer+1}" for mer in range(self.no_mers) ]

        if angle_y.endswith('₁₃'):
            y_columns_of_interest = [ f"{angle_y} mers {mer+1}, {mer+2}" for mer in range(self.no_mers-1) ]
        else:
            y_columns_of_interest = [ f"{angle_y} mer {mer+1}" for mer in range(self.no_mers) ]

        x_cols = [ self.get_colnum_by_meaning(x) for x in x_columns_of_interest ]
        y_cols = [ self.get_colnum_by_meaning(x) for x in y_columns_of_interest ]

        for c in x_cols:
            x = self.dataframe.iloc[:, c]
            x_data = np.concatenate((x_data, x))

        for c in y_cols:
            y = self.dataframe.iloc[:, c]
            y_data = np.concatenate((y_data, y))

        n_datapoints = min( len(x_data), len(y_data) )

        plt.subplots()
        plt.hist2d(x, y, bins=numbins, range=[[-180,180],[-180,180]], norm=mpl.colors.LogNorm(vmin=0.1, vmax=100), cmap=plt.cm.YlOrRd)
        plt.colorbar(format=mtick.ScalarFormatter())
        plt.xlabel(angle_x)
        plt.ylabel(angle_y)
        plt.title(f"{str(self).replace('analysis','')}, n={n_datapoints}")
        plot_filepath = os.path.join(plotdir,f"{self}_hist2D_{angle_x}_{angle_y}.png")
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.clf()
        plt.close()

    def get_entropy(self, xcol: int, ycol: int):
        """
        computes entropy from histogram of xcol vs. ycol
        """
        y = self.dataframe.iloc[:, ycol]
        x = self.dataframe.iloc[:, xcol]

        h = np.histogram2d(x, y, bins=10)
        h_vec = np.concatenate(h[0])
        h_norm = h_vec / sum(h_vec)
        # use molar gas constant R = 8.314
        return 8.314 * entropy(h_norm)

    def list_columns(self):
        """
        This is for debug purposes. Display column index and column meaning for all columns.
        """
        print(self)
        for index, column in enumerate(self.columns):
            print(f"{index+1}: {column}")

class ExperimentalData:
    """
    Gathers all tab files from a data directory, converts them to experiments and offers some manipulation capabilities on these.
    """
    plot_dpi = 350
    def __init__(self, data_path: str):
        experiment_file_names = glob.glob(f"{data_path}/*.tab")
        self.experiments = [Experiment(fp) for fp in experiment_file_names]

    def call_method_by_criteria(self, method, criteria, *args):
        """
        Calls method with arguments in underlying experiments if criteria are met
        """
        affected_experiments = self.choose_experiments(criteria)
        for e in affected_experiments:
            getattr(e,method)(*args)

    def choose_experiments(self, criteria):
        return [ e for e in self.experiments if all([ getattr(e,k) in criteria[k] for k in criteria.keys()]) ]

    def get_entropy_percentiles(self, criteria, angle1, angle2, percentiles):
        chosen_experiments = self.choose_experiments(criteria)
        no_mers = chosen_experiments[0].no_mers
        entropies = np.array([])
        if angle1.endswith('₁₃') and angle2.endswith('₁₃'):
            for mer in range(no_mers-1):
                entropies = np.concatenate((entropies, np.array(self.hist_of_entropy(criteria, f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}"))))
        else:
            for mer in range(no_mers):
                entropies = np.concatenate((entropies, np.array(self.hist_of_entropy(criteria, f"{angle1} mer {mer+1}", f"{angle2} mer {mer+1}"))))
        return [ np.percentile(entropies, p) for p in percentiles ]

    def entropy_distribution_percentiles(self, criteria, angle1: str, angle2: str, plotdir: str):
        """  compute percentiles of the histogram of entropies """
        chosen_experiments = self.choose_experiments(criteria)
        assert len(chosen_experiments)
        no_mers = chosen_experiments[0].no_mers

        median_entropy = []
        entropy_perc5 = []
        entropy_perc95 = []

        if angle1.endswith('₁₃') and angle2.endswith('₁₃'):
            for mer in range(no_mers-1):
                entropies = np.array(self.hist_of_entropy(criteria, f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}"))
                median_entropy.append(np.median(entropies))
                entropy_perc5.append(np.percentile(entropies, 5))
                entropy_perc95.append(np.percentile(entropies, 95))
                myxlabel = "mer n, n+1"
                x_axis_ticks = list(range(1, no_mers))
        else:
            for mer in range(no_mers):
                entropies = np.array(self.hist_of_entropy(criteria, f"{angle1} mer {mer+1}", f"{angle2} mer {mer+1}"))
                median_entropy.append(np.median(entropies))
                entropy_perc5.append(np.percentile(entropies, 5))
                entropy_perc95.append(np.percentile(entropies, 95))
                myxlabel = "mer n"
                x_axis_ticks = list(range(1, no_mers+1))

        mytitle = f"entropy perc {chosen_experiments[0].complex} {chosen_experiments[0].ion} {chosen_experiments[0].chain}"
        myylabel = f"entropy  {angle1} vs. {angle2}"

        plt.plot(x_axis_ticks, median_entropy, "o--", color="red", label="median")
        plt.plot(x_axis_ticks, entropy_perc5, ":", color="red", label="5 and 95 perc.")
        plt.plot(x_axis_ticks, entropy_perc95, ":", color="red", label=None)
        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plot_filepath = os.path.join(plotdir, f"{mytitle.replace(' ','_')}.png")
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.clf()

        return median_entropy

    def entropy_distribution_realisations(self, criteria, angle1: str, angle2: str, plotdir: str):
        """  compute percentiles of the histogram of entropies """
        chosen_experiments = self.choose_experiments(criteria)
        no_mers = chosen_experiments[0].no_mers
        no_struct = len(chosen_experiments)
        entropies = []

        if angle1.endswith('₁₃'):
            first_mers = list(range(1, no_mers))
            myxlabel = "mers x, x+1"
            for mer in range(no_mers-1):
                entropies.append(np.array(self.hist_of_entropy(criteria, f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}")))
        else:
            first_mers = list(range(1, no_mers+1))
            myxlabel = "mer"
            for mer in range(no_mers):
                entropies.append(np.array(self.hist_of_entropy(criteria, f"{angle1} mer {mer+1}", f"{angle2} mer {mer+1}")))

        mytitle = f"entropy reals {chosen_experiments[0].complex} {chosen_experiments[0].ion} {chosen_experiments[0].chain}"
        myylabel = f"entropy  {angle1} vs. {angle2}"

        for j in range(no_struct):
            color = 5 / 6 * (1 - j / no_struct) * np.array((1, 1, 1))
            e = [entropies[i][j] for i in range(len(entropies))]
            plt.plot(first_mers, e, label=f"{j+1}", color=color)

        plt.legend()
        plt.title(mytitle)
        plt.xlabel(myxlabel)
        plt.ylabel(myylabel)
        plot_filepath = os.path.join(plotdir, f"{mytitle.replace(' ','_')}.png")
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.clf()

    def hist_of_entropy(self, criteria, xcolumn: str, ycolumn: str, plotdir: str = None):
        """ compute histogram of entropy over realisations """
        chosen_experiments = self.choose_experiments(criteria)
        xcol = chosen_experiments[0].get_colnum_by_meaning(xcolumn)
        ycol = chosen_experiments[0].get_colnum_by_meaning(ycolumn)
        entropies = [experiment.get_entropy(xcol, ycol) for experiment in chosen_experiments]
        if plotdir:
            plotFile = os.path.join(plotdir,f"{chosen_experiments[0]}hist{xcol}_{ycol}.png")
            mytitle = f"{chosen_experiments[0]} {xcol} {ycol}"
            myxlabel = f"entropy {xcolumn} vs. {ycolumn}"
            myylabel = "frequency"

            plt.subplots()
            plt.hist(entropies, bins=5)
            plt.title(mytitle)
            plt.xlabel(myxlabel)
            plt.ylabel(myylabel)

            plt.savefig(plotFile, dpi=self.plot_dpi)
            plt.close()
        return entropies

    def aggregate_plot(self, criteria, plotdir):
        labels = []
        data = []
        assert len(criteria['complex']) == 1
        for ion in criteria['ion']:
            for chain in criteria['chain']:
                myCriteria = { 'ion': ion, 'chain': chain, 'complex': criteria['complex'][0] }
                for a1, a2 in [ ("ϕ₁₄","ψ₁₄") , ("ϕ₁₃","ψ₁₃") ]:
                    labels.append(f"{ion} {a1}{a2}")
                    data.append(self.get_entropy_percentiles( myCriteria, a1, a2, [5,50,95]))

        plot_filepath = os.path.join(plotdir,f"aggplot_{criteria['complex'][0]}.png")

        fig, ax = plt.subplots(1,1)

        plt.title(f"Aggregate Plot {criteria['complex'][0]}")

        x = range(len(labels))

        y = [ f[0] for f in data ]
        ax.plot(x, y, ".", color="red", label="p5")
        y = [ f[1] for f in data ]
        ax.plot(x, y, "o", color="red", label="median")
        y = [ f[2] for f in data ]
        ax.plot(x, y, ".", color="red", label="p95")

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation='vertical')
        ax.set_ylabel("Entropy")
        ax.set_xlabel("Ion")
        make_axes_area_auto_adjustable(ax)
        plt.legend()
        plt.savefig(plot_filepath, dpi=self.plot_dpi)
        plt.close()
