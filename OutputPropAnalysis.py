"""
Module for fragmentation propensity analysis on output files. Moved out of generic output analysis
container after that ended up having too many conflicting purposes/etc.
#author: Dan Polasky
#date: 11/16/2018
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import itertools
import os
import OutputAnalysis_v2
import tkinter
from tkinter import filedialog


def main_frag_propensities(files, plot_tmp_bool, extension):
    """
    Method to compute fragmentation propensities by amino acid before/after (and pairs) and
    save outputs to graph and csv files. Uses .hits files as input.
    :param files: (list of strings) full system paths to .hits files to analyze
    :param plot_tmp_bool: (bool) if True, make summary ('prop effects') plot using the parse_tmp_filename method
    :param extension: output plot extension (e.g. '.png)
    :return: void
    """
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


    outputdir = filedialog.askdirectory(title='Choose Output Folder')
    # Save output where the user wants
    # outputdir = os.path.dirname(files[0])

    before_dicts, after_dicts, pair_dicts = [], [], []

    for index, hits_file in enumerate(files):
        print('Analyzing file {} of {}'.format(index + 1, len(files)))
        sitelist = OutputAnalysis_v2.load_hits_file(hits_file)
        protein_seq = OutputAnalysis_v2.get_protein_seq(sitelist)

        output_filename = os.path.basename(hits_file).rstrip('.hits')
        csv_filepath = os.path.join(outputdir, output_filename + '_props.csv')

        dict_before, dict_after, dict_pairs = compute_propensities(sitelist, protein_seq, amino_acids)
        before_dicts.append(dict_before)
        after_dicts.append(dict_after)
        pair_dicts.append(dict_pairs)
        plot_propensities_1d(dict_before, outputdir, output_filename, 'frag_Nterm_to_res', extension)
        plot_propensities_1d(dict_after, outputdir, output_filename, 'frag_Cterm_to_res', extension)
        save_prop_csv(dict_before, dict_after, dict_pairs, amino_acids, csv_filepath, output_filename)
        save_prop_file(dict_before, dict_after, dict_pairs, amino_acids, hits_file)

    if plot_tmp_bool:
        plot_prop_effects(before_dicts, after_dicts, files, outputdir, extension)

    save_prop_csv_list(before_dicts, after_dicts, files, amino_acids, outputdir)


def compute_propensities(sitelist, protein_seq, amino_acids):
    """
    Compute fragmentation propensity for each amino acid and amino acid pair and return
    as dictionary (key = AA or pair, value = number of hits)
    :param sitelist: list of Fragment Site containers with hits
    :type sitelist: list[FragmentSite]
    :param amino_acids: list of strings corresponding to all amino acids present in protein seqs
    :param protein_seq: protein sequence (string)
    :return: dictionary of (hit total intensity, num hits, norm'd intensity, norm'd by seq frequency) by AA before, after, and by pair
    """
    # initialize dictionaries to hold hit counts
    dict_before = {key: [0, 0, 0, 0] for key in amino_acids}
    dict_after = {key: [0, 0, 0, 0] for key in amino_acids}
    dict_pairs = {key: [0, 0, 0, 0] for key in itertools.product(amino_acids, amino_acids)}

    # compute total intensity for relative/percent calculations
    total_int = 0
    for site in sitelist:
        for hit in site.hits:
            total_int += float(hit.exp_ion.pkar_cluster)

    # get total intensity at this site and save
    for site in sitelist:
        site_total_int = 0
        for hit in site.hits:
            site_total_int += float(hit.exp_ion.pkar_cluster)

        if site_total_int > 0:
            # add hit intensities and counts to the appropriate AAs in the dict
            res_fragd_cterm, res_fragd_nterm = get_site_residues(site, protein_seq)
            if res_fragd_cterm is not None:
                dict_after[res_fragd_cterm][0] += site_total_int
                dict_after[res_fragd_cterm][1] += 1
                dict_after[res_fragd_cterm][2] += site_total_int / total_int
            if res_fragd_nterm is not None:
                dict_before[res_fragd_nterm][0] += site_total_int
                dict_before[res_fragd_nterm][1] += 1
                dict_before[res_fragd_nterm][2] += site_total_int / total_int
            if res_fragd_cterm is not None and res_fragd_nterm is not None:
                dict_pairs[(res_fragd_nterm, res_fragd_cterm)][0] += site_total_int
                dict_pairs[(res_fragd_nterm, res_fragd_cterm)][1] += 1
                dict_pairs[(res_fragd_nterm, res_fragd_cterm)][2] += site_total_int / total_int
            # if res_fragd_cterm is not None:
            #     dict_before[res_fragd_cterm][0] += site_total_int
            #     dict_before[res_fragd_cterm][1] += 1
            #     dict_before[res_fragd_cterm][2] += site_total_int / total_int
            # if res_fragd_nterm is not None:
            #     dict_after[res_fragd_nterm][0] += site_total_int
            #     dict_after[res_fragd_nterm][1] += 1
            #     dict_after[res_fragd_nterm][2] += site_total_int / total_int
            # if res_fragd_cterm is not None and res_fragd_nterm is not None:
            #     dict_pairs[(res_fragd_cterm, res_fragd_nterm)][0] += site_total_int
            #     dict_pairs[(res_fragd_cterm, res_fragd_nterm)][1] += 1
            #     dict_pairs[(res_fragd_cterm, res_fragd_nterm)][2] += site_total_int / total_int

    # generate data normalized by sequence frequency for each AA
    for residue, data_list in dict_before.items():
        res_freq_ratio = protein_seq.count(residue) / float(len(protein_seq))
        if res_freq_ratio == 0:
            data_list[3] = data_list[2]
        else:
            data_list[3] = data_list[2] / res_freq_ratio
    for residue, data_list in dict_after.items():
        res_freq_ratio = protein_seq.count(residue) / float(len(protein_seq))
        if res_freq_ratio == 0:
            data_list[3] = data_list[2]
        else:
            data_list[3] = data_list[2] / res_freq_ratio
    for res_pair, data_list in dict_pairs.items():
        residue = '{}{}'.format(res_pair[0], res_pair[1])
        res_freq_ratio = protein_seq.count(residue) / float(len(protein_seq))
        if res_freq_ratio == 0:
            data_list[3] = data_list[2]
        else:
            data_list[3] = data_list[2] / res_freq_ratio

    return dict_before, dict_after, dict_pairs


def get_site_residues(site, protein_seq):
    """
    Get the residues (single letter AA codes) before/after the cleavage site this Site represents.
    :param site: FragmentSite
    :type site: FragmentSite
    :param protein_seq: full protein sequence corresponding to this site
    :return: previous AA residue (N-term to cleavage site), next AA residue (C-term to cleavage site)
    """
    if site.term == 'N':
        prev_res = site.seq[-1]
        try:
            # next_res = protein_seq[site.seq_index + 1]
            next_res = protein_seq[site.seq_index]
        except IndexError:
            # Whole sequence - no next residue
            next_res = None
    else:
        try:
            prev_res = protein_seq[site.seq_index - 1]
        except IndexError:
            # whole sequence - no previous residue
            prev_res = None
        next_res = site.seq[0]
    return prev_res, next_res


def plot_prop_effects(before_dicts, after_dicts, files, outputdir, extension):
    """
    Make plots of acidic and proline effects for a list of hits files
    :param before_dicts: list of before dicts (frag to c-term)
    :param after_dicts: list of after dicts (frag to n-term)
    :param files: list of filepaths (to use for indexes)
    :param extension: output plot extension (e.g. '.png)
    :param outputdir: directory in which to save output
    :return: void
    """
    # Generate before and after plot information, then plot it
    lys_pct_tups = []
    acidic_pct_tups = []
    for index, before_dict in enumerate(before_dicts):
        # cv = index * 5 + 10
        num_tmp = parse_tmp_filename(files[index])

        acidic_int = before_dict['D'][0] + before_dict['E'][0]
        lys_int = before_dict['K'][0]
        total_int = np.sum([x[0] for x in list(before_dict.values())])
        acidic_pct = acidic_int / total_int * 100
        lys_pct = lys_int / total_int * 100

        lys_pct_tups.append((num_tmp, lys_pct))
        acidic_pct_tups.append((num_tmp, acidic_pct))

    # after plots
    pro_tups = []
    for index, after_dict in enumerate(after_dicts):
        # cv = index * 5 + 10
        num_tmp = parse_tmp_filename(files[index])

        proline_int = after_dict['P'][0]
        total_int = np.sum([x[0] for x in list(after_dict.values())])
        proline_pct = proline_int / total_int * 100

        pro_tups.append((num_tmp, proline_pct))

    prop_plot_helper(lys_pct_tups, '# TMPs', 'Percent Lysine', '_lysine' + extension, outputdir)
    prop_plot_helper(acidic_pct_tups, '# TMPs', 'Percent Acidic', '_acidic' + extension, outputdir)
    prop_plot_helper(pro_tups, '# TMPs', 'Percent Proline', '_proline' + extension, outputdir)

    prop_plot_helper_avg(lys_pct_tups, '# TMPs', 'Percent Lysine', '_lysine_avg' + extension, outputdir)
    prop_plot_helper_avg(acidic_pct_tups, '# TMPs', 'Percent Acidic', '_acidic_avg' + extension, outputdir)
    prop_plot_helper_avg(pro_tups, '# TMPs', 'Percent Proline', '_proline_avg' + extension, outputdir)


def prop_plot_helper(list_scatter_tups, xlabel, ylabel, plot_name, outputdir):
    """
    Actually generate output plot given a list of (x-value, y-value) tuples and associated
    plot information
    :param list_scatter_tups: list of (x-value, y-value) to plot (typically (# TMPs, % in channel))
    :param xlabel: x axis label
    :param ylabel: y axis label
    :param plot_name: plot filename (short)
    :param outputdir: path to directory in which to save output
    :return: void
    """
    plt.clf()
    plt.figure(dpi=300)
    x_vals = [plot_tup[0] for plot_tup in list_scatter_tups]
    y_vals = [plot_tup[1] for plot_tup in list_scatter_tups]
    plt.scatter(x_vals, y_vals, marker='s')

    # plt.ylim([30, 90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plot_path = os.path.join(outputdir, plot_name)
    plt.savefig(plot_path)
    plt.close()


def prop_plot_helper_avg(list_scatter_tups, xlabel, ylabel, plot_name, outputdir):
    """
    Alternative helper method to prop_plot helper. Instead of plotting each point individually,
    averages all points at same TMP (etc) value to generate single points with error bars.
    :param list_scatter_tups: list of (x-value, y-value) to plot (typically (# TMPs, % in channel))
    :param xlabel: x axis label
    :param ylabel: y axis label
    :param plot_name: plot filename (short)
    :param outputdir: path to directory in which to save output
    :return: void
    """
    # First, reorganize the input list into a dict by # TMPs for easy averaging/etc
    tmp_dict = {}
    for plot_tup in list_scatter_tups:
        if plot_tup[0] in tmp_dict.keys():
            tmp_dict[plot_tup[0]].append(plot_tup[1])
        else:
            tmp_dict[plot_tup[0]] = [plot_tup[1]]

    # Average all values at each # TMPs (/etc) and generate a plot with error bars
    plt.clf()
    plt.figure(dpi=300)
    x, y, xerr, yerr = [], [], [], []
    for tmp_num, val_list in sorted(tmp_dict.items(), key=lambda i: i[0]):
        avg = np.average(val_list)
        std_dev = np.std(val_list)
        x.append(tmp_num)
        y.append(avg)
        xerr.append(0)
        yerr.append(std_dev)
    plt.errorbar(x=x, y=y, xerr=xerr, yerr=yerr, marker='s', ls='none')

    plt.xlabel(xlabel, fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plot_path = os.path.join(outputdir, plot_name)
    plt.savefig(plot_path)
    plt.close()


def parse_tmp_filename(filename):
    """
    convenient handler for Mike/my TMP data. Not intended for general use
    :param filename: string
    :return: int number of TMPs
    """
    splits = os.path.basename(filename).split('_')
    # num_tmps = 0
    # for split in splits:
    #     if split.startswith('T') and len(split) == 3:
    #         num_tmps = int(split[1:])

    tmp_split = int(splits[4])
    # num_tmps = (tmp_split - 994) // 13            # Serf 8+
    # num_tmps = round((tmp_split - 1425) / 17.3)   # Ubq 6+
    # num_tmps = round((tmp_split - 2140) / 26.0)   # Ubq 4+
    # num_tmps = round((tmp_split - 1789) / 13.0)   # Lys 8+
    num_tmps = round((tmp_split - 1590) / 13.0)     # CytC 8+
    # num_tmps = round((tmp_split - 2022) / 11.6)     # B-lac 9+
    if num_tmps < 0:
        num_tmps = 0
    return num_tmps


def plot_propensities_1d(prop_dict, output_path, hitsfile_name, output_label, extension, norm=True):
    """
    Generate a bar graph of propensities for hits and total intensity for provided
    dictionary
    :param prop_dict: dict keys=AAs, values=[total int, hit count]
    :param output_path: directory in which to save output
    :param hitsfile_name: name of file to save
    :param output_label: before/after to append to filename
    :param extension: output plot extension (e.g. '.png)
    :param norm: if True, uses relative intensity rather than raw/total
    :return: void
    """
    plt.clf()
    plt.figure(dpi=300)

    indices = np.arange(len(prop_dict.keys()))

    # plot intensities
    if norm:
        intensities = [x[2] for x in list(prop_dict.values())]
        intensities_freq = [x[3] for x in list(prop_dict.values())]
    else:
        intensities = [x[0] for x in list(prop_dict.values())]
        intensities_freq = [x[0] for x in list(prop_dict.values())]
    # standard plot
    plt.bar(indices, intensities, label='Total Intensity')
    plt.xticks(indices, list(prop_dict.keys()))
    plotname = '{}_{}{}'.format(hitsfile_name, output_label, extension)
    plotfilename = os.path.join(output_path, plotname)
    plt.savefig(plotfilename)

    # normalized by frequence plot
    plt.bar(indices, intensities_freq, label='Intensity : AA Freq Ratio')
    plt.xticks(indices, list(prop_dict.keys()))
    plotname = '{}_{}{}'.format(hitsfile_name, output_label + '_freq', extension)
    plotfilename = os.path.join(output_path, plotname)
    plt.savefig(plotfilename)

    plt.close()


def save_prop_file(dict_before, dict_after, dict_pairs, amino_acids, hits_filename):
    """
    Pickle a .prop file for the given hits file by saving the dictionaries to file as a list
    :param dict_before: dict keys=AAs, values=[total int, hit count]
    :param dict_after: dict keys=AAs, values=[total int, hit count]
    :param dict_pairs: dict keys=(AA, AA), values=[total int, hit count]
    :param amino_acids: list of strings corresponding to all amino acids present in protein seqs
    :param hits_filename: hits file from which data was drawn
    :return: void
    """
    dict_list = [dict_before, dict_after, dict_pairs, amino_acids]
    prop_file = hits_filename.rstrip('.hits') + '.prop'
    with open(prop_file, 'wb') as savefile:
        pickle.dump(dict_list, savefile)


def load_prop_file(filename):
    """
    Load a .prop file and return the dictionaries
    :param filename: full path to file to load
    :return: dict before, after, pairs, amino acids
    """
    with open(filename, 'rb') as readfile:
        dict_list = pickle.load(readfile)

    return dict_list[0], dict_list[1], dict_list[2], dict_list[3]


def save_prop_csv(dict_before, dict_after, dict_pairs, amino_acids, output_filename, hits_filename):
    """
    Save all outputs to csv from dictionaries
    :param dict_before: dict keys=AAs, values=[total int, hit count]
    :param dict_after: dict keys=AAs, values=[total int, hit count]
    :param dict_pairs: dict keys=(AA, AA), values=[total int, hit count]
    :param amino_acids: list of strings corresponding to all amino acids present in protein seqs
    :param output_filename: full path to save output
    :param hits_filename: hits file from which data was drawn
    :return: void
    """
    output_strings = ['Hits File: {}'.format(hits_filename), 'Before AA', 'AA,Total Int,Hit Count']
    for key, val_tup in dict_before.items():
        line = '{},{},{}'.format(key, val_tup[0], val_tup[1])
        output_strings.append(line)
    output_strings.extend(['After AA', 'AA,Total Int,Hit Count'])
    for key, val_tup in dict_after.items():
        line = '{},{},{}'.format(key, val_tup[0], val_tup[1])
        output_strings.append(line)

    index = 0
    aa_index = 0
    header = 'Pairs,'
    header += ','.join(amino_acids)
    output_strings.append(header)
    while index < len(dict_pairs.items()):
        int_vals = [str(x[0]) for x in list(dict_pairs.values())[index: index + 20]]
        line = amino_acids[aa_index] + ','
        line += ','.join(int_vals)
        output_strings.append(line)
        index += 20
        aa_index += 1

    with open(output_filename, 'w') as outfile:
        for string in output_strings:
            outfile.write(string + '\n')


def save_prop_csv_list(list_dict_before, list_dict_after, list_hitsfiles, amino_acids, output_dir):
    """
    Similar to saving individual propensity csvs, but saves ALL into one file for easy excel manipulation
    :param list_dict_before: list of dicts, keys=AAs, values=[total int, hit count]
    :param list_dict_after: list of dicts, keys=AAs, values=[total int, hit count]
    :param amino_acids: list of strings corresponding to all amino acids present in protein seqs
    :param list_hitsfiles: list of hits file names from which data was drawn
    :param output_dir: directory in which to save output
    :return: void
    """
    output_file = os.path.join(output_dir, '_ComboProps.csv')
    header = 'AA,' + ','.join([os.path.basename(x) for x in list_hitsfiles]) + '\n'
    before_totals = {key: '' for key in amino_acids}
    before_counts = {key: '' for key in amino_acids}
    before_percents = {key: '' for key in amino_acids}
    after_totals = {key: '' for key in amino_acids}
    after_counts = {key: '' for key in amino_acids}
    after_percents = {key: '' for key in amino_acids}
    dict_list = [before_percents, before_totals, before_counts, after_percents, after_totals, after_counts]

    # assemble information sorted by amino acid and type
    # for file_index in range(len(list_dict_before)):
    for aa in amino_acids:
        for dict_before in list_dict_before:
            val_tup = dict_before[aa]
            before_totals[aa] += '{},'.format(val_tup[0])
            before_counts[aa] += '{},'.format(val_tup[1])
            before_percents[aa] += '{},'.format(val_tup[2])
        for dict_after in list_dict_after:
            val_tup = dict_after[aa]
            after_totals[aa] += '{},'.format(val_tup[0])
            after_counts[aa] += '{},'.format(val_tup[1])
            after_percents[aa] += '{},'.format(val_tup[2])

    # write to file
    with open(output_file, 'w') as outfile:
        outfile.write('Combined Propensities\n')
        for string_dict in dict_list:
            outfile.write(header)
            for key, val in sorted(string_dict.items()):
                outfile.write('{},{}\n'.format(key, val))
            outfile.write('\n')


def combine_props_to_csv(list_prop_files, outputdir):
    """
    Convenience method to combine a list of previously saved .props files into a single output
    csv
    :param list_prop_files: list of .props files to load
    :return: void
    """
    dicts_before, dicts_after = [], []
    aas = []
    # outputdir = os.path.dirname(list_prop_files[0])

    for prop_file in list_prop_files:
        before, after, pairs, aas = load_prop_file(prop_file)
        dicts_before.append(before)
        dicts_after.append(after)

    save_prop_csv_list(dicts_before, dicts_after, list_prop_files, aas, outputdir)


if __name__ == '__main__':
    root = tkinter.Tk()
    root.withdraw()

    hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])

    save_extension = '.png'
    # OutputAnalysis_v2.main_seq_cov(hitsfiles, save_extension)
    # main_frag_propensities(hitsfiles, plot_tmp_bool=True, extension=save_extension)
    main_frag_propensities(hitsfiles, plot_tmp_bool=False, extension=save_extension)

