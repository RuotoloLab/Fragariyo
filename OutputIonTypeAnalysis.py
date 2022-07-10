"""
Module for ion type analysis/comparison/etc of hits files.
#author: CRR (Based on DP's code from OutputAnalysis_v2)
#date: 1/25/2019
"""
import tkinter
from tkinter import filedialog
import pickle
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from PyQt5 import QtWidgets
from terminalFragmentor_Main import FragmentSite
import pyperclip

CONFIG_FILE = 'config.txt'  # config file for saving last directory for fancy filedialog

# Set the plot styling globally (Later on there will be more customization-CRR)
plt.style.use('ggplot')


def get_protein_seq(sitelist):
    """
    Search the list of sites to find the longest sequence string; this will be the protein sequence.
    Typically in the first position, but will search all just to make sure.
    :param sitelist: list of FragmentSite containers
    :type sitelist: list[FragmentSite]
    :return: (string) protein sequence
    """
    seq = ''
    for site in sitelist:
        if len(site.full_protein_seq) > len(seq):
            seq = site.full_protein_seq
    return seq


def get_data(config_file):
    """
    Load folders of data using custom FileDialog class
    :param config_file: path to the config file with the initial directory for the file chooser
    :return: list of strings of full system folder paths to the folders chosen, updated input_dir
    """
    input_dir = get_last_dir(config_file)

    app = QtWidgets.QApplication(sys.argv)
    ex = FileDialog(input_dir)
    ex.show()
    app.exec_()
    files = ex.selectedFiles()

    new_base_dir = os.path.dirname(files[0])
    save_config(config_file, new_base_dir)
    return files


def get_last_dir(config_file):
    """
    parse the config file for the last directory used, to use as the initial directory when
    opening the file chooser.
    :param config_file: text file with a single directory (full system path) and nothing else
    :return: (string) directory path
    """
    with open(config_file, 'r') as config:
        line = config.readline()
        return line


def save_config(config_file, new_base_dir):
    """
    Update the config file with a new directory name
    :param config_file: file path to update
    :param new_base_dir: information to save in the config file
    :return: void
    """
    with open(config_file, 'w') as config:
        config.write(new_base_dir)


class FileDialog(QtWidgets.QFileDialog):
    """
    File chooser for raw data, created after extensive searching on stack overflow
    """
    def __init__(self, input_dir, *args):
        QtWidgets.QFileDialog.__init__(self, *args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.DirectoryOnly)
        self.setDirectory(input_dir)

        self.tree = self.findChild(QtWidgets.QTreeView)
        self.tree.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)


def load_hits_file(filename):
    """
    Load a saved (pickled) .hits file and return the stored list of FragmtSites
    :param filename: full path to file to load
    :return: list of FragmtSite containers with hits information
    :rtype: list[FragmentSite]
    """
    with open(filename, 'rb') as loadfile:
        sitelist = pickle.load(loadfile)
    return sitelist


def unpack_hitsfile(files):
    """
    Obtains the site list form each hits file inputed by user (Max two)
    :param files: (list of strings) full system paths to .hits files to analyze
    :return: site list, protein_seq and output filename for each .hits file
    """
    # Intitialized an index so only two proteins will be compared
    hitsfileIndex = 1

    for hits_file in files:
        # For the first file
        if hitsfileIndex == 1:
            sitelist1 = load_hits_file(hits_file)
            protein_seq1 = get_protein_seq(sitelist1)
            output_filename1 = os.path.basename(hits_file).rstrip('.hits')
            hitsfileIndex += 1
        #  For the second file
        elif hitsfileIndex == 2:
            sitelist2 = load_hits_file(hits_file)
            protein_seq2 = get_protein_seq(sitelist2)
            output_filename2 = os.path.basename(hits_file).rstrip('.hits')
        else:
            print("No more than two files can be combined!!!")

    return sitelist1, protein_seq1, output_filename1, sitelist2, protein_seq2, output_filename2

def compute_combo_ions(proteinseq1, proteinseq2, site_list1, site_list2, ionTypes_list):
    """
    Combine and compute fragment ion types occurrences for two .hits files.
    :param protein_seq: str, from .hits files
    :param site_list: list of FragmtSite objects containing hit information
    :param ionTypes_list: What type of ions to analyzed for
    :return: dict, ion types (keys) and their occurrences (values)
    """

    assert(proteinseq1 == proteinseq2), "Different proteins cannot be combined!"

    # Initialize directory
    ionType_dict = {ion: [] for ion in ionTypes_list}

    # print(ionType_dict)

    for site1 in site_list1:
        for hit in site1.hits:
            if hit.exp_ion.pkar_cluster > 2000:
                # The fragment ions are being sorted by ion type, then by position
                if hit.thy_ion.ion_type[0] == 'a':
                    ionType_dict['a'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'b':
                    ionType_dict['b'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'y':
                    ionType_dict['y'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'x':
                    ionType_dict['x'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'c':
                    ionType_dict['c'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'z':
                    ionType_dict['z'].append(hit.thy_ion.ion_type_indx)

                # print(hit.pass_num)
                # print(hit.cal_error)
                # print(hit.thy_ion)
                # print(hit.thy_ion.iontype)
                # print(hit.thy_ion.mods)

    for site2 in site_list2:
        for hit in site2.hits:
            # The last pass of multipass search contains the "leftover ions"
            if hit.exp_ion.pkar_cluster > 2000:
                if hit.thy_ion.ion_type[0] == 'a':
                    ionType_dict['a'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'b':
                    ionType_dict['b'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'y':
                    ionType_dict['y'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'x':
                    ionType_dict['x'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'c':
                    ionType_dict['c'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'z':
                    ionType_dict['z'].append(hit.thy_ion.ion_type_indx)

    # A fragment ion can be added twice:
    # It was a hit with/without a modification or it had a different change state
    # Therefore the sequence position lists were made into sets
    for ion in ionTypes_list:
        ionType_dict[ion] = set(ionType_dict[ion])
    # The amount of fragment ion occurrences are the length of the unique sets of sequence positions per ion type
        ionType_dict[ion] = len(ionType_dict[ion])

    return ionType_dict


def compute_ions(site_list, ionTypes_list):
    """
    Compute fragment ion types occurrences for the experimental .hits files.
    :param site_list: list of FragmtSite objects containing hit information
    :param ionTypes_list: What type of ions to analyzed for
    :return: dict, ion types (keys) and their occurrences (values)
    """

    ionType_dict = {ion: [] for ion in ionTypes_list}

    # print(ionType_dict)

    for site in site_list:

        for hit in site.hits:
            # Small area clusters indicate false peak id
            #IONS = ['a','b','c','c-dot', 'c-1','x','y','z','z-dot', 'z+1', 'a+1']
            if hit.exp_ion.pkar_cluster > 2000:
                if hit.thy_ion.ion_type[0] == 'a':
                    ionType_dict['a'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'b':
                    ionType_dict['b'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'y':
                    ionType_dict['y'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'x':
                    ionType_dict['x'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'c':
                    ionType_dict['c'].append(hit.thy_ion.ion_type_indx)

                elif hit.thy_ion.ion_type[0] == 'z':
                    ionType_dict['z'].append(hit.thy_ion.ion_type_indx)


    for ion in ionTypes_list:
        ionType_dict[ion] = set(ionType_dict[ion])
        ionType_dict[ion] = len(ionType_dict[ion])

    return ionType_dict

def compute_mods(site_list, ionTypes_list):
    """
    Compute fragment modifications for the experimental .hits files.
    :param site_list: list of FragmtSite objects containing hit information
    :param ionTypes_list: What type of ions to analyzed for
    :return: dict, ion types (keys) and their occurrences (values)
    """

    mods_dict = {ion: [] for ion in ionTypes_list}

    # print(ionType_dict)

    for site in site_list:

        for hit in site.hits:
            # Small area clusters indicate false peak id
            if hit.exp_ion.pkar_cluster > 2000:

                #Add mods present in each hit per ion_type
                if hit.thy_ion.ion_type[0] == 'a' and len(hit.thy_ion.thy_mods) != 0:
                    mods_dict['a'].append(f"Ion: {hit.thy_ion.ion_type_indx} Mods: {hit.thy_ion.thy_mods}  Cysmods: {hit.thy_ion.cysmods}")

                elif hit.thy_ion.ion_type[0] == 'b' and len(hit.thy_ion.thy_mods) != 0:
                    mods_dict['b'].append(f"Ion: {hit.thy_ion.ion_type_indx} Mods: {hit.thy_ion.thy_mods}  Cysmods: {hit.thy_ion.cysmods}")

                elif hit.thy_ion.ion_type[0] == 'y' and len(hit.thy_ion.thy_mods) != 0:
                    mods_dict['y'].append(f"Ion: {hit.thy_ion.ion_type_indx} Mods: {hit.thy_ion.thy_mods}  Cysmods: {hit.thy_ion.cysmods}")

                elif hit.thy_ion.ion_type[0] == 'x' and len(hit.thy_ion.thy_mods) != 0:
                    mods_dict['x'].append(f"Ion: {hit.thy_ion.ion_type_indx} Mods: {hit.thy_ion.thy_mods}  Cysmods: {hit.thy_ion.cysmods}")

                elif hit.thy_ion.ion_type[0] == 'c' and len(hit.thy_ion.thy_mods) != 0:
                    mods_dict['c'].append(f"Ion: {hit.thy_ion.ion_type_indx} Mods: {hit.thy_ion.thy_mods}  Cysmods: {hit.thy_ion.cysmods}")

                elif hit.thy_ion.ion_type[0] == 'z' and len(hit.thy_ion.thy_mods) != 0:
                    mods_dict['z'].append(f"Ion: {hit.thy_ion.ion_type_indx} Mods: {hit.thy_ion.thy_mods}  Cysmods: {hit.thy_ion.cysmods}")


    # for ion in ionTypes_list:
    #     ionType_dict[ion] = set(ionType_dict[ion])
    #     ionType_dict[ion] = len(ionType_dict[ion])

    return mods_dict

def compute_neutloss(site_list, ionTypes_list):
    """
    Compute fragment modifications for the experimental .hits files.
    :param site_list: list of FragmtSite objects containing hit information
    :param ionTypes_list: What type of ions to analyzed for
    :return: dict, ion types (keys) and their occurrences (values)
    """

    nl_dict = {ion: [] for ion in ionTypes_list}

    # print(ionType_dict)

    for site in site_list:

        for hit in site.hits:
            # Small area clusters indicate false peak id
            if hit.exp_ion.pkar_cluster > 2000:

                #Add mods present in each hit per ion_type
                if hit.thy_ion.ion_type[0] == 'a' and len(hit.thy_ion.neutlosses) != 0:
                    nl_dict['a'].append(f"Ion:{hit.thy_ion.ion_type_indx} Neutral Loss = {hit.thy_ion.neutlosses}")

                elif hit.thy_ion.ion_type[0] == 'b' and len(hit.thy_ion.thy_mods) != 0:
                    nl_dict['b'].append(f"Ion:{hit.thy_ion.ion_type_indx} Neutral Loss = {hit.thy_ion.neutlosses}")

                elif hit.thy_ion.ion_type[0] == 'y' and len(hit.thy_ion.thy_mods) != 0:
                    nl_dict['y'].append(f"Ion:{hit.thy_ion.ion_type_indx} Neutral Loss = {hit.thy_ion.neutlosses}")

                elif hit.thy_ion.ion_type[0] == 'x' and len(hit.thy_ion.thy_mods) != 0:
                    nl_dict['x'].append(f"Ion:{hit.thy_ion.ion_type_indx} Neutral Loss = {hit.thy_ion.neutlosses}")

                elif hit.thy_ion.ion_type[0] == 'c' and len(hit.thy_ion.thy_mods) != 0:
                    nl_dict['c'].append(f"Ion:{hit.thy_ion.ion_type_indx} Neutral Loss = {hit.thy_ion.neutlosses}")

                elif hit.thy_ion.ion_type[0] == 'z' and len(hit.thy_ion.thy_mods) != 0:
                    nl_dict['z'].append(f"Ion:{hit.thy_ion.ion_type_indx} Neutral Loss = {hit.thy_ion.neutlosses}")


    # for ion in ionTypes_list:
    #     ionType_dict[ion] = set(ionType_dict[ion])
    #     ionType_dict[ion] = len(ionType_dict[ion])

    return nl_dict

# def compute_combo_mods(proteinseq1, proteinseq2, site_list1, site_list2, ionTypes_list):
#     """
#     Combine and compute fragment mods occurrences for two .hits files.
#     :param protein_seq: str, from .hits files
#     :param site_list: list of FragmtSite objects containing hit information
#     :param ionTypes_list: What type of ions to analyzed for
#     :return: dict, ion types (keys) and their occurrences (values)
#     """
#
#     assert(proteinseq1 == proteinseq2), "Different proteins cannot be combined!"
#
#     # Initialize directory
#     modscombo_dict = {ion: [] for ion in ionTypes_list}
#
#     # print(ionType_dict)
#
#     for site1 in site_list1:
#         for hit in site1.hits:
#             if hit.exp_ion.pkar_cluster > 2000:
#                 # The fragment ions are being sorted by ion type, then by position
#                 if hit.thy_ion.ion_type[0] == 'a':
#                     modscombo_dict['a'].append(f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'b':
#                     modscombo_dict['b'].append(f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'y':
#                     modscombo_dict['y'].append(f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'x':
#                     modscombo_dict['x'].append(f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'c':
#                     modscombo_dict['c'].append(f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'z':
#                     modscombo_dict['z'].append(f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 # print(hit.pass_num)
#                 # print(hit.cal_error)
#                 # print(hit.thy_ion)
#                 # print(hit.thy_ion.iontype)
#                 # print(hit.thy_ion.mods)
#
#     for site2 in site_list2:
#         for hit in site2.hits:
#             # The last pass of multipass search contains the "leftover ions"
#             if hit.exp_ion.pkar_cluster > 2000:
#                 if hit.thy_ion.ion_type[0] == 'a':
#                     modscombo_dict['a'].append(
#                         f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'b':
#                     modscombo_dict['b'].append(
#                         f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'y':
#                     modscombo_dict['y'].append(
#                         f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'x':
#                     modscombo_dict['x'].append(
#                         f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'c':
#                     modscombo_dict['c'].append(
#                         f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#                 elif hit.thy_ion.ion_type[0] == 'z':
#                     modscombo_dict['z'].append(
#                         f"Mods: {hit.thy_ion.thy_mods} Ion: {hit.thy_ion.ion_type_indx} Cysmods: {hit.thy_ion.cysmods}")
#
#     # A fragment ion can be added twice:
#     # It was a hit with/without a modification or it had a different change state
#     # Therefore the sequence position lists were made into sets
#     for ion in ionTypes_list:
#         modscombo_dict[ion] = set(modscombo_dict[ion])
#     # The amount of fragment ion occurrences are the length of the unique sets of sequence positions per ion type
#         modscombo_dict[ion] = len(modscombo_dict[ion])
#
#     return ionType_dict

def general_seq_plot(prop_dict, outputdir, plotitle, xtitle=None, ytitle=None, extension='.png'):
    """
    Plot any data as a function of protein sequence position (from N-terminus)
    :param prop_dict: A dictionary containing each ion type (a,b,y,c,z,x) and how many times do they ocurr in a hits file
    :param outputdir: directory in which to save output
    :param plotitle: title of plot, required, as file is saved with this name
    :param xtitle: x axis title (optional)
    :param ytitle: y axis title (optional)
    :param extension: With what extension will the plot be saved
    :return: none (void)
    """

    plt.clf()
    plt.figure(dpi=300)

    # The bar graph will have bars as the number of ion types analyzed
    indices = np.arange(len(prop_dict.keys()))
    # Make bar graph with the number of ions that were hits
    plt.bar(indices, prop_dict.values(), label= f"Total bonds broken = {sum(prop_dict.values())}")
    # Each bar's label
    plt.xticks(indices, list(prop_dict.keys()))

    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(plotitle)
    plt.legend(loc='best')

    # Named the output differently from the plots produced in the OutputAnalysis_v2 module
    plotname = plotitle + '_iontypes' + extension
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()

def general_seq_plot_combo(prop_dict, outputdir, plotitle, xtitle=None, ytitle=None, extension='.png'):
    """
    Plot any data as a function of protein sequence position (from N-terminus)
    :param prop_dict: A dictionary containing each ion type (a,b,y,c,z,x) and how many times do they ocurr in a hits file
    :param outputdir: directory in which to save output
    :param plotitle: title of plot, required, as file is saved with this name
    :param xtitle: x axis title (optional)
    :param ytitle: y axis title (optional)
    :param extension: With what extension will the plot be saved
    :return: none (void)
    """

    plt.clf()
    plt.figure(dpi=300)

    for dic in prop_dict:

        # The bar graph will have bars as the number of ion types analyzed
        indices = np.arange(len(dic.keys()))
        # Make bar graph with the number of ions that were hits
        plt.bar(indices, prop_dict.values(), label= f"Total bonds broken = {sum(prop_dict.values())}")
        # Each bar's label
        plt.xticks(indices, list(prop_dict.keys()))

    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(plotitle)
    plt.legend(loc='best')

    # Named the output differently from the plots produced in the OutputAnalysis_v2 module
    plotname = plotitle + '_iontypes-combo' + extension
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()

def ion_analysis(files, extension='.png', plotting = False, combination = False, outfile = False, outfile_extension=None):
    """
    Run ion typw analysis on the provided list of .hits files
    :param files: (list of strings) full system paths to .hits files to analyze
    :param extension: output plot extension (e.g. '.png)
    :param plotting: A plot will be created or only ion type calculations will be output in the stdout
    :param combination: weather or not to combine two hits files
    :param outfile: Outputs a file of type(see Outfile_extension). It has all the ion type counts in written form.
    :param outfile_extension: .csv
    :return: void
    """

    ionTypes = ['a', 'b', 'y', 'x', 'z', 'c']

    # Save output to where the .hits files are
    # outdir = os.path.dirname(files[0])

    # Save output where the user wants
    outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(outdir)

    #Create an output string with headers
    output_str = 'Sample name,a,b,y,x,c,z\n'
    outname = ""

    # Compute number of fragment ions and based on argument plot and/or combine
    if combination:

        sitelist1, protein_seq1, output_filename1, sitelist2, protein_seq2, output_filename2 = unpack_hitsfile(files)
        ion_type_dict_combo = compute_combo_ions(protein_seq1, protein_seq2, sitelist1, sitelist2,  ionTypes)
        output_title = output_filename1 + '_' + output_filename2
        print(output_title)
        print(ion_type_dict_combo)
        output_str += f"{output_title},{ion_type_dict_combo['a']},{ion_type_dict_combo['b']},{ion_type_dict_combo['y']}" \
            f",{ion_type_dict_combo['x']},{ion_type_dict_combo['c']},{ion_type_dict_combo['z']}\n"

        if plotting:
            general_seq_plot(ion_type_dict_combo, outputdir=outdir, plotitle=output_title, xtitle='ion type',
                             ytitle='Count', extension=extension)

        outname += f"{output_title}_iontypes"
    else:

        for index, hits_file in enumerate(files):
            print('analyzing file {} of {}'.format(index + 1, len(files)))
            sitelist = load_hits_file(hits_file)
            output_title = os.path.basename(hits_file).rstrip('.hits')
            print(output_title)
            ion_type_dict = compute_ions(sitelist, ionTypes)
            print(ion_type_dict)
            output_str += f"{output_title},{ion_type_dict['a']},{ion_type_dict['b']},{ion_type_dict['y']}," \
                f"{ion_type_dict['x']},{ion_type_dict['c']},{ion_type_dict['z']}\n"

            if plotting:
                general_seq_plot(ion_type_dict, outputdir=outdir, plotitle=output_title, xtitle='ion type',
                             ytitle='Count', extension=extension)

            outname = "Ion_Analysis"

    if outfile and outfile_extension is not None:
        output = open(outname + outfile_extension, 'w')
        output.write(output_str)
        output.close()
        print(output_str)


def mods_analysis(files, extension='.png', plotting = False, combination = False, outfile = False, outfile_extension=None):
    """
    Run ion typw analysis on the provided list of .hits files
    :param files: (list of strings) full system paths to .hits files to analyze
    :param extension: output plot extension (e.g. '.png)
    :param plotting: A plot will be created or only ion type calculations will be output in the stdout
    :param combination: weather or not to combine two hits files
    :param outfile:
    :param outfile_extension: .csv
    :return: void
    """

    ionTypes = ['a', 'b', 'y', 'x', 'c', 'z']

    # Save output to where the .hits files are
    # outdir = os.path.dirname(files[0])

    # Save output where the user wants
    outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(outdir)

    #Create an output string with headers
    output_str = 'Sample name\ta\tb\ty\tx\tc\tz\n'

    # Compute number of fragment ions and based on argument plot and/or combine
    if combination:

        sitelist1, protein_seq1, output_filename1, sitelist2, protein_seq2, output_filename2 = unpack_hitsfile(files)
        ion_type_dict_combo = compute_combo_ions(protein_seq1, protein_seq2, sitelist1, sitelist2,  ionTypes)
        output_title = output_filename1 + '_' + output_filename2 + 'combo'
        print(output_title)
        print(ion_type_dict_combo)
        output_str += f"{output_title}\t{ion_type_dict_combo['a']}\t{ion_type_dict_combo['b']}\t{ion_type_dict_combo['y']}\t" \
                f"{ion_type_dict_combo['x']}\t{ion_type_dict_combo['c']}\t{ion_type_dict_combo['z']}\n"


        if plotting:
            general_seq_plot(ion_type_dict_combo, outputdir=outdir, plotitle=output_title, xtitle='ion type',
                             ytitle='Count', extension=extension)
    else:

        for index, hits_file in enumerate(files):
            print('analyzing file {} of {}'.format(index + 1, len(files)))
            sitelist = load_hits_file(hits_file)
            output_title = os.path.basename(hits_file).rstrip('.hits')
            print(output_title)
            ion_type_dict = compute_mods(sitelist, ionTypes)
            print(ion_type_dict)
            output_str += f"{output_title}\t{ion_type_dict['a']}\t{ion_type_dict['b']}\t{ion_type_dict['y']}\t" \
                f"{ion_type_dict['x']}\t{ion_type_dict['c']}\t{ion_type_dict['z']}\n"

            if plotting:
                general_seq_plot(ion_type_dict, outputdir=outdir, plotitle=output_title, xtitle='ion type',
                             ytitle='Count', extension=extension)
    if outfile and outfile_extension is not None:
        output = open("Mods_analysis" + outfile_extension, 'w')
        output.write(output_str)
        output.close()
        print(output_str)

def neutloss_analysis(files, extension='.png', plotting = False, combination = False, outfile = False, outfile_extension=None):
    """
    Run ion typw analysis on the provided list of .hits files
    :param files: (list of strings) full system paths to .hits files to analyze
    :param extension: output plot extension (e.g. '.png)
    :param plotting: A plot will be created or only ion type calculations will be output in the stdout
    :param combination: weather or not to combine two hits files
    :param outfile:
    :param outfile_extension: .csv
    :return: void
    """

    ionTypes = ['a', 'b', 'y', 'x', 'c', 'z']

    # Save output to where the .hits files are
    # outdir = os.path.dirname(files[0])

    # Save output where the user wants
    outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(outdir)

    #Create an output string with headers
    output_str = 'Sample name\ta\tb\ty\tx\tc\tz\n'

    # Compute number of fragment ions and based on argument plot and/or combine
    if combination:

        sitelist1, protein_seq1, output_filename1, sitelist2, protein_seq2, output_filename2 = unpack_hitsfile(files)
        ion_type_dict_combo = compute_combo_ions(protein_seq1, protein_seq2, sitelist1, sitelist2,  ionTypes)
        output_title = output_filename1 + '_' + output_filename2 + 'combo'
        print(output_title)
        print(ion_type_dict_combo)
        output_str += f"{output_title}\t{ion_type_dict_combo['a']}\t{ion_type_dict_combo['b']}\t{ion_type_dict_combo['y']}\t" \
                f"{ion_type_dict_combo['x']}\t{ion_type_dict_combo['c']}\t{ion_type_dict_combo['z']}\n"


        if plotting:
            general_seq_plot(ion_type_dict_combo, outputdir=outdir, plotitle=output_title, xtitle='ion type',
                             ytitle='Count', extension=extension)
    else:

        for index, hits_file in enumerate(files):
            print('analyzing file {} of {}'.format(index + 1, len(files)))
            sitelist = load_hits_file(hits_file)
            output_title = os.path.basename(hits_file).rstrip('.hits')
            print(output_title)
            ion_type_dict = compute_neutloss(sitelist, ionTypes)
            print(ion_type_dict)
            output_str += f"{output_title}\t{ion_type_dict['a']}\t{ion_type_dict['b']}\t{ion_type_dict['y']}\t" \
                f"{ion_type_dict['x']}\t{ion_type_dict['c']}\t{ion_type_dict['z']}\n"

            if plotting:
                general_seq_plot(ion_type_dict, outputdir=outdir, plotitle=output_title, xtitle='ion type',
                             ytitle='Count', extension=extension)
    if outfile and outfile_extension is not None:
        output = open("NeutralLosses_analysis" + outfile_extension, 'w')
        output.write(output_str)
        output.close()
        print(output_str)

# Unused for now - CRR
def batch_main_seq_cov():
    """
    Batch processing method for main sequence coverage. Runs main_seq_cov on all .hits files
    in each folder selected in series
    :return: void
    """
    batch_folders = get_data(CONFIG_FILE)
    for index, batch_folder in enumerate(batch_folders):
        print('Starting batch {} of {}...'.format(index + 1, len(batch_folders)))
        files = [os.path.join(batch_folder, x) for x in os.listdir(batch_folder) if x.endswith('.hits')]

def input_PyMOL(files, replicates=False):
    """
    Takes in .hits files and outputs the hits that they have in common, so that these residues can be colored in pyMOL
    ONLY TO BE USED IN REPLICATES FILES!!!
    :param files: (list of strings) full system paths to .hits files to analyze
    :return: A set of residue numbers which are common hits in all the files
    """

    #Choose output folder
    # outdir = filedialog.askdirectory(title='Choose Output Folder')
    residueList_pyMOL = []
    for index, hits_file in enumerate(files):
        seqResidue_list = []
        print('analyzing file {} of {}'.format(index + 1, len(files)))
        sitelist = load_hits_file(hits_file)
        output_title = os.path.basename(hits_file).rstrip('.hits')
        print(f"Filename: {output_title}")
        protein_seq = get_protein_seq(sitelist)
        print(f"Protein sequence: {protein_seq}")
        print(f"Protein length: {len(protein_seq)}\n")
        for site in sitelist:
            #print(site)
            if site.term == 'C':
                for hit in site.hits:
                    #print("\t", hit)
                    # For the c-term extra math is needed to get the residue number (validated with Protein Prospector)
                    ctermResidue = len(protein_seq)-(hit.thy_ion.ion_type_indx)+1
                    #print("\t", f" Ion type: {hit.thy_ion.iontype}({hit.thy_ion.non_mod_length}), residue number: {ctermResidue}")
                    seqResidue_list.append(ctermResidue)

            elif site.term == 'N':
                for hit in site.hits:
                    #print("\t", hit)
                    #print("\t", f" Ion type: {hit.thy_ion.iontype}({hit.thy_ion.non_mod_length}), residue number: {hit.thy_ion.non_mod_length}")
                    seqResidue_list.append(hit.thy_ion.ion_type_indx)

            else:
                print('INVALID TERMINAL! Terminal was ' + site.term)

        residueList_pyMOL.append(seqResidue_list)
        resultstring = str(set(seqResidue_list))
        # print(resultstring)
        result_one = resultstring.replace(", ", "+")
        result_two = result_one.replace("{", " ")
        result = result_two.replace("}", " ")
        print(result + "\n")

    if replicates:
        #A list of each sample list of residue number covered
        print(residueList_pyMOL)
        #Select the first list and convert it to a set
        s = set(residueList_pyMOL[0])
        setCompare = residueList_pyMOL[1:]
        for x in setCompare:
            s.intersection_update(x)
        #the intial set now contains common elements between itself and the list
        #print(s)
        resultstring = str(s)
        #print(resultstring)
        result_one = resultstring.replace(", ", "+")
        result_two = result_one.replace("{", " ")
        result = result_two.replace("}", " ")
        print(result)




if __name__ == '__main__':
    root = tkinter.Tk()
    root.withdraw()

    hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])



    # ion_analysis(hitsfiles, extension='.png', plotting=True, combination=False, outfile=True, outfile_extension='.csv')
    # mods_analysis(hitsfiles, extension='.png', plotting=False, combination=False, outfile=True, outfile_extension='.tsv')
    # neutloss_analysis(hitsfiles, extension='.png', plotting=False, combination=False, outfile=True,
    #               outfile_extension='.tsv')
    #
    input_PyMOL(hitsfiles, replicates=True)



