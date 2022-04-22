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
from scipy.stats import entropy

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

    initial_columns = [
        ColumnMeaning("Time", "ps"),                            #First column, index 0
        ColumnMeaning("Total energy of the system", "kJ/mol"),
        ColumnMeaning("Total bond energy", "kJ/mol"),
        ColumnMeaning("Total angle energy", "kJ/mol"),
        ColumnMeaning("Total dihedral energy", "kJ/mol"),
        ColumnMeaning("Total planar energy", "kJ/mol"),
        ColumnMeaning("Total Coulomb energy", "kJ/mol"),
        ColumnMeaning("Total Van der Waals energy", "kJ/mol"),  #Column index 7
    ]

    def __init__(self, filepath: str):
        split_items = re.split('_|\.',os.path.basename(filepath)) # pylint: disable=W1401
        self.complex, self.num_realisation, self.chain, self.ion, _ = split_items
        self.dataframe = pd.read_csv(filepath, sep=";", skipfooter=2, engine="python")
        self.angles = self.angles_config[(self.complex,self.chain)]
        self.no_mers = 23
        self.columns = self.initial_columns

        if self.complex == 'Albumin+HA':
            if self.chain == 'analysis':
                for mer in range(self.no_mers):
                    for angle in self.angles:
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))
            elif self.chain == 'side chain':
                for angle in self.angles:
                    for mer in range(self.no_mers):
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))

        if self.complex == 'Albumin+CS6':
            if self.chain == 'analysis':
                for mer in range(self.no_mers):
                    for angle in self.angles:
                        self.columns.append(ColumnMeaning(f"{angle} mers {mer+1}, {mer+2}", "deg"))
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

        x_columns_of_interest = [ f"{angle_x} mers {mer+1}, {mer+2}" for mer in range(self.no_mers) ]
        y_columns_of_interest = [ f"{angle_y} mers {mer+1}, {mer+2}" for mer in range(self.no_mers) ]

        x_cols = [ self.get_colnum_by_meaning(x) for x in x_columns_of_interest ]
        y_cols = [ self.get_colnum_by_meaning(x) for x in y_columns_of_interest ]

        for c in x_cols:
            x = self.dataframe.iloc[:, c]
            x_data = np.concatenate((x_data, x))

        for c in y_cols:
            y = self.dataframe.iloc[:, c]
            y_data = np.concatenate((y_data, y))

        plt.subplots()
        plt.hist2d(x, y, bins=numbins, range=[[-180,180],[-180,180]], cmap=plt.cm.Reds)
        plt.colorbar(format=mtick.FormatStrFormatter("%.1e")) #FIXME scientific scale
        plt.xlabel(angle_x)
        plt.ylabel(angle_y)
        plt.title(f"{self} - Subsequent mers angles")
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
        for mer in range(no_mers):
            entropies = np.concatenate((entropies, np.array(self.hist_of_entropy(criteria, f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}"))))
        return [ np.percentile(entropies, p) for p in percentiles ]

    def entropy_distribution_percentiles(self, criteria, angle1: str, angle2: str, plotdir: str):
        """  compute percentiles of the histogram of entropies """
        chosen_experiments = self.choose_experiments(criteria)
        assert len(chosen_experiments)
        no_mers = chosen_experiments[0].no_mers
        first_mers = list(range(1, no_mers+1))

        median_entropy = []
        entropy_perc5 = []
        entropy_perc95 = []

        for mer in range(no_mers):
            entropies = np.array(self.hist_of_entropy(criteria, f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}"))
            median_entropy.append(np.median(entropies))
            entropy_perc5.append(np.percentile(entropies, 5))
            entropy_perc95.append(np.percentile(entropies, 95))

        mytitle = f"entropy perc {chosen_experiments[0].complex} {chosen_experiments[0].ion} {chosen_experiments[0].chain}"
        myylabel = f"entropy  {angle1} vs. {angle2}"
        myxlabel = "mers x, x+1"

        plt.plot(first_mers, median_entropy, "o--", color="red", label="median")
        plt.plot(first_mers, entropy_perc5, ":", color="red", label="5 and 95 perc.")
        plt.plot(first_mers, entropy_perc95, ":", color="red", label=None)
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
        first_mers = list(range(1, no_mers+1))
        entropies = []

        for mer in range(no_mers):
            entropies.append(np.array(self.hist_of_entropy(criteria, f"{angle1} mers {mer+1}, {mer+2}", f"{angle2} mers {mer+1}, {mer+2}")))

        mytitle = f"entropy reals {chosen_experiments[0].complex} {chosen_experiments[0].ion} {chosen_experiments[0].chain}"
        myylabel = f"entropy  {angle1} vs. {angle2}"
        myxlabel = "mers x, x+1"

        for j in range(no_struct):
            color = 5 / 6 * (1 - j / no_struct) * np.array((1, 1, 1))
            e = [entropies[i][j] for i in range(no_mers)]
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

    def plot21(self, criteria): #TODO this method name is given after issue number, proper naming needed
        for myComplex in criteria['complex']:
            for ion in criteria['ion']:
                for chain in criteria['chain']:
                    myCriteria = { 'ion': ion, 'chain': chain, 'complex': complex }
                    for a1, a2 in [ ("ϕ₁₄","ψ₁₄") , ("ϕ₁₃","ψ₁₃") ]:
                        myLabel = f"{myComplex} {ion} {chain} {a1}{a2}"
                        print(f"{myLabel=}")
                        print(self.get_entropy_percentiles( myCriteria, a1, a2, [5,50,95]))
