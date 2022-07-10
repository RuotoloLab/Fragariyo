"""
Module for analysis/comparison/etc of hits files.
#author: DP
#date: 8/14/2018
Updated 06/01/20 to analyze output from Annotatorv3
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
from terminalFragmentor_Main import print_hits
#Not directyy used by the module, but used by FragmentSite
from terminalFragmentor_Main import ThyIon
from Parameter_Parser_terminal import load_hits_file
from terminalFragmentor_Main import save_fragments



CONFIG_FILE = 'config.txt'  # config file for saving last directory for fancy filedialog

def compare_main():
    """
    :return: 
    """

    file1 = filedialog.askopenfilename(title='Load CNTRL file', filetypes=[('Hits', '.hits')])
    file2 = filedialog.askopenfilename(title='Load TEMPO file', filetypes=[('Hits', '.hits')])

    # Transform paths into filenames
    file1_outname = file1.split('/')
    file1_outname = file1_outname[-1].strip('-avg.hits')
    file2_outname = file2.split('/')
    file2_outname = file2_outname[-1].strip('-avg.hits')
    # Create output filename
    output_title = f"{file2_outname}-{file1_outname}"

    # Choose output directory
    outdir = filedialog.askdirectory(title='Choose Output Folder')
    # Run comparison
    comparing_files(file1, file2, outdir, output_title=output_title)



def sitelist_to_sitedict(file):
    """
    :param:
    :return:
    """
    sitelist = load_hits_file(file)
    file = {}
    for site in sitelist:
        # print(site.hits)
        thyions = []
        for hit in site.hits:
            # print(hit.thy_ion)
            thyions.append(hit)

        file[site] = thyions
    protein_seq = get_protein_seq(sitelist)

    return file, protein_seq


def comparing_files(file1, file2, outdir, output_title=None):
    """
    Function to check what ions from a TEMPO file are in the CNTRL, whihc has been searched with TEMPO search space
    :param files: .hits files
    :param outdir: Directory to save the output files
    :param output_title: name of output file
    :return: void
    """
    # Initialize variables
    final_list = []

    print('Comparing Files...')

    dict1, protein_seq1 = sitelist_to_sitedict(file1)
    dict2, protein_seq2 = sitelist_to_sitedict(file2)

    dict3 = {}

    if protein_seq1 == protein_seq2:
        for site in dict1:
            print(site)
            # print(dict1[site])
            # print(dict2[site])

            thyion_ls = []
            for hit1 in dict1[site]:
                thyion_ls.append(hit1.thy_ion.mz_mono)


            print(f"To remove = {thyion_ls}")
            to_keep = []
            for hit2 in dict2[site]:
                print(hit2.thy_ion.mz_mono)
                if hit2.thy_ion.mz_mono not in thyion_ls:
                    print(hit2.thy_ion.mz_mono)
                    to_keep.append(hit2)
            print(f"To keep = {to_keep}")
            dict3[site] = to_keep



    else:
        print("Files must be from the same sequence!")

    print(dict3)

    for site in dict2:
        site.hits = dict3[site]
        final_list.append(site)


    output_name = os.path.join(outdir, output_title + '-avg.hits')
    csv_outname = output_name.strip('.hits')
    csv_outname += '.csv'
    print_hits(protein_seq1, final_list, csv_outname)

    print(output_name)
    with open(output_name, 'wb') as picklefile:
        pickle.dump(final_list, picklefile)
        picklefile.close()

def merging_files(files, outdir, output_title=None):
    """
    This function "merges" files that correspond to different passes of the same analysis (Same sample).
     Usually for complex searches (e.g albumin with disulfide bonds)
    :param files: .hits files
    :param outdir: Directory to save the output files
    :param output_title: name of output file
    :return: void
    """
    # Initilaize variables
    final_list = []
    counter = 0
    hit_counter =0
    mergefiles_dict = {}
    protein_seq = []

    for index, hits_file in enumerate(files):

        print('analyzing file {} of {}'.format(index + 1, len(files)))
        # print(hits_file)
        sitelist = load_hits_file(hits_file)
        print(f"Lenthg of sitelist: {len(sitelist)}")
        protein_seq.append(get_protein_seq(sitelist))

        for site in sitelist:
            # print(site.hits)
            thyions = []
            for hit in site.hits:
                # print(hit.thy_ion)
                thyions.append(hit)
                hit_counter += 1

            #Loop to initialize and empty dictionary
            if counter == 0:
                mergefiles_dict[site] = thyions
            else:
                mergefiles_dict[site].extend(thyions)
        counter += 1

    print(counter)
    print(f"before: {mergefiles_dict}")

    for site in mergefiles_dict:
        # print(site)
        # print(avgfiles_dict[site])
        # print(len(avgfiles_dict[site]))
        # print(site.hits)
        site.hits = list(mergefiles_dict[site])
        # print(site.hits)

        final_list.append(site)

    print(final_list)
    print(f"hit counter = {hit_counter}")

    output_name = os.path.join(outdir, output_title + '-merged.hits')
    csv_outname = output_name.strip('.hits')
    csv_outname += '.csv'
    print_hits(protein_seq[-1], final_list, csv_outname)

    print(output_name)
    with open(output_name, 'wb') as picklefile:
        pickle.dump(final_list, picklefile)
        picklefile.close()



def averaging_files(files, outdir, output_title=None):
    """
    This function "averages" replicate files into finding the ions that appeared across all replicate files
    :param files: .hits files
    :param outdir: Directory to save the output files
    :param output_title: name of output file
    :return: void
    """
    # Initilaize variables
    final_list = []
    counter = 0
    avgfiles_dict = {}
    protein_seq = []

    for index, hits_file in enumerate(files):

        print('analyzing file {} of {}'.format(index + 1, len(files)))
        # print(hits_file)
        sitelist = load_hits_file(hits_file)
        # print(sitelist)
        protein_seq.append(get_protein_seq(sitelist))

        for site in sitelist:
            # print(site.hits)
            thyions = []
            for hit in site.hits:
                # print(hit.thy_ion)
                thyions.append(hit)

            if counter == 0:
                avgfiles_dict[site] = thyions
            else:
                avgfiles_dict[site].extend(thyions)
        counter += 1

    print(counter)
    print(f"before: {avgfiles_dict}")

    for site in avgfiles_dict:
        """
        This code is used to count how many times an ion occurs in a site after aggregating all replicates files 
        """

        unique_ions = set(avgfiles_dict[site])
        # print(unique_ions)

        ions_remove = []
        for ion in unique_ions:
            # print(f"{site} {ion} = {avgfiles_dict[site].count(ion)}")

            # Ions need to be in all replicates to not be removed
            if avgfiles_dict[site].count(ion) < counter:
                ions_remove.append(ion)
        # print(ions_remove)

        avgfiles_dict[site] = unique_ions.difference(set(ions_remove))

    print(f"after: {avgfiles_dict}")

    for site in avgfiles_dict:
        # print(site)
        # print(avgfiles_dict[site])
        # print(len(avgfiles_dict[site]))
        # print(site.hits)
        site.hits = list(avgfiles_dict[site])
        # print(site.hits)

        final_list.append(site)

    print(final_list)

    output_name = os.path.join(outdir, output_title + '-avg.hits')
    csv_outname = output_name.strip('.hits')
    csv_outname += '.csv'
    print_hits(protein_seq[-1], final_list, csv_outname)

    print(output_name)
    with open(output_name, 'wb') as picklefile:
        pickle.dump(final_list, picklefile)
        picklefile.close()

def merging_batch():
    """
    Code is a bit strange because when I started I was taking stuff out of object and putting them in a dictionary
    In the process I ended up learning more about object and at the end I included the whole object in the dictionary
    :param batch: Boolean, if True a .csv file with the desired files to be average will be taken from desired input directories
    and will be saved in desired output folders with desired names
    :return:
    """

    batch_file = filedialog.askopenfilename(title='Merge_Batch Hits Files', filetypes=[('CSV File', '.csv')])

    files_list = []
    indir = ''
    with open(batch_file, 'r') as batch:
        lines = list(batch)

        for line in lines:

            if line.startswith("#"):
                files_sub = []
                files_list.append(files_sub)
                continue
            else:
                # print(line)
                line = line.strip("\n")
                splits = line.split(',')
                if splits[0]:
                    indir = splits[0]
                    outdir = splits[2]
                    output_title = splits[3]
                    files_sub.append(outdir)
                    files_sub.append(output_title)
                    files_sub.append(f"{indir}\{splits[1]}")
                else:
                    files_sub.append(f"{indir}\{splits[1]}")
            print(files_sub)

    print(indir)
    # p = os.access(indir, os.F_OK)
    # print(p)
    print(files_list)
    print(outdir)
    # print(output_title)

    for file_set in files_list:
        outdir = file_set[0]
        output_title = file_set[1]
        files = file_set[2:]
        merging_files(files, outdir, output_title)


def average_main(batch = False):
    """
    Code is a bit strange because when I started I was taking stuff out of object and putting them in a dictionary
    In the process I ended up learning more about object and at the end I included the whole object in the dictionary
    :param batch: Boolean, if True a .csv file with the desired files to be average will be taken from desired input directories
    and will be saved in desired output folders with desired names
    :return:
    """

    if batch:
        batch_file = filedialog.askopenfilename(title='Batch Hits Files', filetypes=[('CSV File', '.csv')])

        files_list = []
        indir = ''
        with open(batch_file, 'r') as batch:
            lines = list(batch)

            for line in lines:

                if line.startswith("#"):
                    files_sub = []
                    files_list.append(files_sub)
                    continue
                else:
                    # print(line)
                    line = line.strip("\n")
                    splits = line.split(',')
                    if splits [0]:
                        indir =  splits[0]
                        outdir = splits[2]
                        output_title = splits[3]
                        files_sub.append(outdir)
                        files_sub.append(output_title)
                        files_sub.append(f"{indir}\{splits[1]}")
                    else:
                        files_sub.append(f"{indir}\{splits[1]}")
                print(files_sub)

        print(indir)
        # p = os.access(indir, os.F_OK)
        # print(p)
        print(files_list)
        print(outdir)
        # print(output_title)

        for file_set in files_list:
            outdir = file_set[0]
            output_title = file_set[1]
            files = file_set[2:]
            averaging_files(files, outdir, output_title)

    else:
        hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])

        outdir = filedialog.askdirectory(title='Choose Output Folder')


        output_title = os.path.basename(hitsfiles[0]).rstrip('.hits')

        averaging_files(hitsfiles, outdir, output_title)


def main_seq_cov(files, outputdir, extension='.png', iontype_exclusion=False, errors=False, outfile = False, outfile_extension=None, cz_iontypes_only=False):
    """
    Run sequence coverage analysis on the provided list of .hits files
    :param files: (list of strings) full system paths to .hits files to analyze
    :param extension: output plot extension (e.g. '.png)
    :return: void
    """
    # outputdir = os.path.dirname(files[0])
    # It used to set the input directory as the output directory...Me no like it
    # It gets too confusing. Modified to code so the user can choose the output directory
    # outputdir = filedialog.askdirectory(title='Choose Output Folder')

    output_header = 'Sample name, seq_cov, max_int\n'
    output_str = f"{output_header}"


    for index, hits_file in enumerate(files):
        print('analyzing file {} of {}'.format(index + 1, len(files)))
        sitelist = load_hits_file(hits_file)
        protein_seq = get_protein_seq(sitelist)

        output_filename = os.path.basename(hits_file).rstrip('.hits')
        csv_filepath = os.path.join(outputdir, output_filename + '_analysis.csv')

        # Compute sequence coverage percentage, need to set seqcovnum argument in general_seq_plot
        seqcovnum = str(seq_cov_number(protein_seq, sitelist))

        # Compute sequence coverage and plot/write to file
        n_ints, c_ints, max_int = compute_seq_cov(protein_seq, sitelist,  norm = False, iontype_exclusion=iontype_exclusion, cz_iontypes_only=cz_iontypes_only)

        output_str += f"{output_filename},{seqcovnum}, {round(max_int)}\n"

        # Plots sequence position vs normalized intensity of sites covered
        general_seq_plot([n_ints, c_ints], ['n', 'c'], outputdir, plotitle=output_filename, seqcovnumber=seqcovnum, xtitle='Sequence position',
                         ytitle='Normed Int', extension=extension)

        # Outputs a .csv file with the info necessary to output the graph above
        general_seq_write([n_ints, c_ints], ['n', 'c'], seq=protein_seq, title='Sequence coverage',
                          new_filename=csv_filepath)

        seq_plot_cterm([c_ints], ['c'], outputdir, plotitle=output_filename)
        seq_plot_nterm([n_ints], ['n'], outputdir, plotitle=output_filename)

        if errors:
            error_list, error_lists_all = seq_error_analysis(protein_seq, sitelist)
            # - There is a ZeroDivisionError
            error_scatter_plot(error_list, outputdir, plotitle=output_filename, xtitle="Sequence Position", ytitle="Avg. Error(ppm)")
            # error_scatter_plot(error_lists_all, outputdir, plotitle=output_filename, xtitle="Sequence Position",
        #                    ytitle="Error(ppm)")

        if outfile and outfile_extension is not None:
            output = open("Sequence_cov" + outfile_extension, 'w')
            output.write(output_str)
            output.close()
        print(output_str)


def seq_cov_overlay(list_hits_files, norm_bool):
    """
    Generate an overlay plot given a set of .hits files. Recommended to use less than 10 or so
    hits files to keep things from getting too crowded
    :param list_hits_files: list of strings (full system paths to hits files to overlay)
    :param norm_bool: whether to normalize the data or not
    :return: void
    """
    outputdir = os.path.dirname(list_hits_files[0])
    plt.clf()
    plt.figure(dpi=600)
    ax1 = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(212, sharex=ax1)

    for hits_file in list_hits_files:
        sitelist = load_hits_file(hits_file)
        protein_seq = get_protein_seq(sitelist)

        output_filename = os.path.basename(hits_file).rstrip('.hits')

        # Compute sequence coverage and plot/write to file
        n_ints, c_ints = compute_seq_cov(protein_seq, sitelist, norm=norm_bool)

        # Plot sequence position vs normalized intensity of sites covered
        seq_indices = np.arange(1, len(n_ints) + 1)
        ax1.plot(seq_indices, n_ints, linewidth=0.8)
        ax1.axes.xaxis.set_visible(False)

        cv_label = parse_funcs_filename(output_filename)
        ax2.plot(seq_indices, c_ints, linewidth=0.8, label=cv_label)

    # rearrange subplots to fit legend at the side
    box1 = ax1.get_position()
    ax1.set_position([box1.x0, box1.y0, box1.width * 0.95, box1.height])
    box2 = ax2.get_position()
    ax2.set_position([box2.x0, box2.y0, box2.width * 0.95, box2.height])
    plt.legend(loc='center left', bbox_to_anchor=(1, 1), fancybox=True, fontsize=8)

    plt.xlabel('Seq Position')
    plt.ylabel('Intensity')

    plotname = 'overlay.png'
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()


def parse_funcs_filename(filename):
    """
    Convert function number from IMTBX to CV (specifically for DP fragment data)
    :param filename: string
    :return: string of format '<voltage> V'
    """
    splits = filename.split('.')
    func_num = int(splits[0][-3:])
    cv = 10 + (func_num - 1) * 5
    return '{}V'.format(cv)


def main_seqcov_iontype(files, outputdir, neutloss = False):
    """
    Run sequence coverage analysis by ion_type on the provided list of .hits files
    :param files: (list of strings) full system paths to .hits files to analyze
    :param outputdir: Directory to save output
    :param neutloss: If neutral losses want to be considered (Pyteomics calculates neut losses both both CID and ExD ions)
    :return: void
    """

    #loop over files
    for index, hits_file in enumerate(files):
        print('analyzing file {} of {}'.format(index + 1, len(files)))
        # Read the pickel file and get protein sequence
        sitelist = load_hits_file(hits_file)
        protein_seq = get_protein_seq(sitelist)

        output_filename = os.path.basename(hits_file).rstrip('.hits')
        csv_filepath = os.path.join(outputdir, output_filename + '_analysis.csv')

        # Compute sequence coverage percentage
        seqcovnum = str(seq_cov_number(protein_seq, sitelist))

        # Compute sequence coverage and plot/write to file
        if neutloss:
            a_ratios, b_ratios, c_ratios, x_ratios, y_ratios, z_ratios = compute_seq_cov_iontype(
                protein_seq, sitelist)
        else:
            a_ratios, b_ratios, c_ratios, x_ratios, y_ratios, z_ratios = compute_seq_cov_iontype_noneutloss(
                protein_seq, sitelist)

        # Plots sequence position vs normalized intensity of sites covered
        general_seq_plot([a_ratios, b_ratios, c_ratios, x_ratios, y_ratios, z_ratios], ['a', 'b','c', 'x', 'y', 'z'], outputdir, plotitle=output_filename, seqcovnumber = seqcovnum,
                         xtitle='Sequence position', ytitle='Normed Int')
        # Outputs a .csv file with the info necessary to output the graph above
        general_seq_write([a_ratios, b_ratios, c_ratios, x_ratios, y_ratios, z_ratios], ['a', 'b','c', 'x', 'y', 'z'], seq=protein_seq, title='Sequence coverage',
                          new_filename=csv_filepath)

        # error_list, error_lists_all = seq_error_analysis(protein_seq, sitelist) - There is a ZeroDivisionError
        # error_scatter_plot(error_list, outputdir, plotitle=output_filename, xtitle=None, ytitle=None)


def main_seq_cov_overlayn(files, outputdir):
    """
    Run sequence coverage analysis on the provided list of .hits files (only n-term half of the protein is considered)
    :param files: (list of strings) full system paths to .hits files to analyze
    :return: void
    """
    # Initialize frame outside of files loop in order to overlay
    plt.figure('overalyn', dpi=200)
    # Clear frame settings, just in case
    plt.clf()
    # each file gets a counter which will shift the x-values so that the overlaying works
    counter = 0
    # loop over files
    for hits_file in files:
        outputfilename = os.path.basename(hits_file).rstrip('.hits')
        # Set counter
        counter += .45
        # Read the pickel file and get protein sequence
        sitelist = load_hits_file(hits_file)
        protein_seq = get_protein_seq(sitelist)
        # Compute sequence coverage percentage
        seqcovnum = str(seq_cov_number(protein_seq, sitelist[:len(sitelist)//2]))
        # Create legend
        legend_str = os.path.basename(hits_file).rstrip('.hits')
        # Creating legend for graph - still work in progress
        #Currently skipping the date....Needs to be more general
        legend = legend_str[13:]

        # Compute sequence coverage and plot/write to file
        n_ints, c_ints = compute_seq_cov(protein_seq, sitelist)

        # It works if output is enclosed in [], I guess making a list
        intensities = [n_ints]
        # Initialize x-values
        indices = np.arange(0, len(intensities[0]))
        #Creat graph
        for data in intensities:
            plt.bar(indices[:len(indices) // 2] + counter, data[:len(data) // 2], width=0.45,  align = 'edge', label=[legend + ' {}%'.format(seqcovnum)], alpha=0.75)

    #Graph details
    plt.xlabel('Sequence Position')
    plt.ylabel('Relative Intensity')
    plt.title('N-term Overlay')
    plt.legend(loc = 'best')

    # Plotname will be asked to the user, but to simplify things...not for now
    plotname = outputfilename + 'Ntermoveraly' + '.png'
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    # plt.show() - It works awesomely (able to zoom in and such), however once I close the plot the program is stuck...- CRR
    plt.close()


def main_seq_cov_overlayc(files, outputdir, outputfilename):
    """
    Run sequence coverage analysis on the provided list of .hits files (only c-term half of the protein is considered)
    :param files: (list of strings) full system paths to .hits files to analyze
    :return: void
    """
    # Initialize frame out side of files loop in order to overlay
    plt.figure('overalyn', dpi=200)
    # Clear frame settings, just in case
    plt.clf()
    # each file gets a counter which will shift the x-values so that the overlaying works
    counter = 0
    # loop over files
    for hits_file in files:
        # Set counter
        counter += 0.45
        # Read the pickel file and get protein sequence
        sitelist = load_hits_file(hits_file)
        protein_seq = get_protein_seq(sitelist)
        # Compute sequence coverage percentage
        seqcovnum = str(seq_cov_number(protein_seq, sitelist[len(sitelist)//2:]))
        # Create legend
        legend_str = os.path.basename(hits_file).rstrip('.hits')
        # Creating legend for graph - still work in progress
        legend = legend_str[13:]

        # Compute sequence coverage and plot/write to file
        n_ints, c_ints = compute_seq_cov(protein_seq, sitelist)

        # It works is out put is enclosed in [], I guess making a list
        intensities = [c_ints]
        # Initialize x-values
        indices = np.arange(0, len(intensities[0]))
        for data in intensities:
            plt.bar(indices[len(indices) // 2:] + counter, data[len(data) // 2:], width=0.45,  align = 'edge', label=[legend + ' {}%'.format(seqcovnum)], alpha=0.75)


    plt.xlabel('Sequence Position')
    plt.ylabel('Relative Intensity')
    plt.title('C-term Overlay')
    plt.legend(loc = 'best')

    # Plotname will be asked to the user, but to simplify things...not for now
    plotname = outputfilename + 'Ctermoveraly' + '.png'
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    # plt.show() - It works awesomely, however once I close the plot the program is stuck... - CRR
    plt.close()

def main_seq_cov_combo(files, outputdir):
    """
    Run sequence coverage analysis on combined lists of two .hits files
    :param files: (list of strings) full system paths to .hits files to analyze
    :return: void
    """
    #From the two .hits files selected obtain specific information to perform the combination
    sitelist1, protein_seq1, output_filename1, sitelist2, protein_seq2, output_filename2 = unpack_hitsfile(files)

    #Combine the two datasets and obtaining the new sitelist and sequence coverag percentage
    sitelist3, seqcovnum3 = combo_sitelists(sitelist1, sitelist2, protein_seq1, protein_seq2)



    output_filename = f'{output_filename1}-{output_filename2}_combo'
    csv_filepath = os.path.join(outputdir, output_filename + '-analysis.csv')

    # Compute sequence coverage and plot/write to file
    n_ints, c_ints, max_int = compute_seq_cov(protein_seq1, sitelist3)
    # Plots sequence position vs normalized intensity of sites covered
    general_seq_plot([n_ints, c_ints], ['n', 'c'], outputdir, plotitle=output_filename, seqcovnumber=seqcovnum3,
                     xtitle='Sequence position', ytitle='Normed Int')
    # Outputs a .csv file with the info necessary to output the graph above
    general_seq_write([n_ints, c_ints], ['n', 'c'], seq=protein_seq1, title='Sequence coverage',
                      new_filename=csv_filepath)

    #Save combined hits files
    save_fragments(sitelist3, output_filename + '.hits')


def combo_sitelists(site_list1, site_list2, protein_seq1, protein_seq2):
    """
    Merge site_list from two different hits_file in order to combine sequencing datasets
    :param sitelist1: Sitelist from hits-file1
    :param sitelist2: Sitelist from hits-file2
    :param protein_seq1: from hist_file1
    :param protein_seq2:  from hist_file2
    :return:
    """
    #Only allow combination of two datasets of the same protein
    if protein_seq1 == protein_seq2:
        # Make a copy of the first site_list
        site_list3 = site_list1[:]

        for site3 in site_list3:
            for site2 in site_list2:
                if site3.term == site2.term:
                    if site3.seq_index == site2.seq_index:
                        site3.hits += site2.hits
        seqcovnum3 = str(seq_cov_number(protein_seq1, site_list3))
    else:
        print("Protein Sequences should be the same!")

    return site_list3, seqcovnum3

def unpack_hitsfile(files):
    """
    Obtains the site list form each hits file inputed by user (Max two)
    :param files: (list of strings) full system paths to .hits files to analyze
    :return: site list, protein_seq and output filename for each .hits file
    """
    hitsfileIndex = 1
    for hits_file in files:
        if hitsfileIndex == 1:
            sitelist1 = load_hits_file(hits_file)
            protein_seq1 = get_protein_seq(sitelist1)
            output_filename1 = os.path.basename(hits_file).rstrip('.hits')
            hitsfileIndex += 1
        elif hitsfileIndex == 2:
            sitelist2 = load_hits_file(hits_file)
            protein_seq2 = get_protein_seq(sitelist2)
            output_filename2 = os.path.basename(hits_file).rstrip('.hits')
        else:
            print("No more than two files can be combined!!!")

    return sitelist1, protein_seq1, output_filename1, sitelist2, protein_seq2, output_filename2

def compute_seq_cov(protein_seq, site_list, norm=False, iontype_exclusion=False, cz_iontypes_only=False):
    """
    Compute and plot sequence coverage for a protein sequence given experimental hit list.
    :param protein_seq: protein amino acid sequence string
    :param site_list: list of FragmtSite objects containing hit information
    :param norm: set to True to return absolute intensities rather than normalized
    :return: n-terminal normalized intensities, c-terminal normalized intensities. Both in lists from n-terminal (ready
    for plotting).
    """
    bad_iontypes = ["z+1", "z+2", "z+3", "z-H2O", "x-H2O", "c-H2O", "z-NH3", "c-NH3", "x-NH3"]
    cz_iontypes = ['z', 'c']
    # find max intensity for normalization
    max_int = get_max_int(site_list)
    if max_int == 0:
        max_int = 1e-10

    # compute normalized intensity at each site and in list by site index
    frag_sites_len = len(protein_seq) + 1  # Extra entry to include the complete sequence
    c_norm_ints = np.zeros(frag_sites_len)
    n_norm_ints = np.zeros(frag_sites_len)
    for site in site_list:

        # get total intensity at this site and normalize it against max found in any site
        total_int = 0

        for hit in site.hits:
            if iontype_exclusion:
                if hit.thy_ion.iontype not in bad_iontypes:
                    total_int += float(hit.exp_ion.pkar_cluster)
            elif cz_iontypes_only:
                if hit.thy_ion.iontype[0] in cz_iontypes:
                    total_int += float(hit.exp_ion.pkar_cluster)
            else:
                total_int += float(hit.exp_ion.pkar_cluster)

        if not norm:
            normed_int = total_int
        else:
            normed_int = float(total_int) / max_int

        # get terminal and correct so C- and N-terminal ions are using same numbering scheme
        site_num = site.get_seq_index(len(protein_seq))
        if site.term == 'N':
            # numbering starts at 1 and goes upwards. Map to starting at 0 and going up (subtract 1)
            n_norm_ints[site_num] = normed_int
        elif site.term == 'C':
            # length of the sequence is referenced to the opposite terminal and needs to be flipped
            c_norm_ints[site_num] = normed_int
        else:
            print('INVALID TERMINAL! Terminal was ' + site.term)

    return n_norm_ints, c_norm_ints, max_int

def compute_seq_cov_iontype_noneutloss(protein_seq, site_list, norm=False):
    """
    Compute and plot sequence coverage by iontype for a protein sequence given experimental hit list.
    :param protein_seq: protein amino acid sequence string
    :param site_list: list of FragmtSite objects containing hit information
    :return: n-terminal normalized intensities, c-terminal normalized intensities. Both in lists from n-terminal (ready
    for plotting).
    """
    # find max intensity for normalization
    max_int = get_max_int(site_list)

    #Create an array for each ion type
    a_ratios = np.zeros(len(protein_seq))
    b_ratios = np.zeros(len(protein_seq))
    c_ratios = np.zeros(len(protein_seq))
    x_ratios = np.zeros(len(protein_seq))
    y_ratios = np.zeros(len(protein_seq))
    z_ratios = np.zeros(len(protein_seq))
    for site in site_list:

        site_num = site.get_seq_index(len(protein_seq))

        if site.term == 'C':
            x_int = 0
            y_int = 0
            z_int = 0

            for hit in site.hits:
                # to no considering neutral losses don't index at 0
                if hit.thy_ion.iontype == 'x':
                    x_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype == 'y':
                    y_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype == 'z' or hit.thy_ion.iontype == 'z-dot':
                    z_int += float(hit.exp_ion.pkar_cluster)
                else:
                    print('INVALID ION-TYPE! Ion-type was {}'.format(hit.thy_ion.iontype))

            if not norm:
                x_ratio = x_int
                x_ratios[site_num - 1] = x_ratio
                y_ratio = y_int
                y_ratios[site_num - 1] = y_ratio
                z_ratio = z_int
                z_ratios[site_num - 1] = z_ratio
            else:
                # get total intensity at site per ion type and normalize it against max found in any site
                x_ratio = float(x_int) / max_int
                x_ratios[site_num - 1] = x_ratio
                y_ratio = float(y_int) / max_int
                y_ratios[site_num - 1] = y_ratio
                z_ratio = float(z_int) / max_int
                z_ratios[site_num - 1] = z_ratio




        elif site.term == 'N':
            a_int = 0
            b_int = 0
            c_int = 0
            # check for a and b fragments and compare cluster total areas
            for hit in site.hits:
                if hit.thy_ion.iontype == 'a':
                    a_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype == 'b':
                    b_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype == 'c' or hit.thy_ion.iontype == 'c-dot':
                    c_int += float(hit.exp_ion.pkar_cluster)
                else:
                    print('INVALID ION-TYPE! Ion-type was {}'.format(hit.thy_ion.iontype))

            if not norm:
                a_ratio = a_int
                a_ratios[site_num - 1] = a_ratio
                b_ratio = b_int
                b_ratios[site_num - 1] = b_ratio
                c_ratio = c_int
                c_ratios[site_num - 1] = c_ratio
            else:
                a_ratio = float(a_int) / max_int
                a_ratios[site_num-1] = a_ratio
                b_ratio = float(b_int) / max_int
                b_ratios[site_num-1] = b_ratio
                c_ratio = float(c_int) / max_int
                c_ratios[site_num-1] = c_ratio
        else:
            print('INVALID TERMINAL! Terminal was ' + site.term)


    return a_ratios, b_ratios, c_ratios, x_ratios, y_ratios, z_ratios

def compute_seq_cov_iontype(protein_seq, site_list, norm =False):
    """
    Compute and plot sequence coverage by iontype for a protein sequence given experimental hit list.
    :param protein_seq: protein amino acid sequence string
    :param site_list: list of FragmtSite objects containing hit information
    :return: n-terminal normalized intensities, c-terminal normalized intensities. Both in lists from n-terminal (ready
    for plotting).
    """
    # find max intensity for normalization
    max_int = get_max_int(site_list)

    #Create an array for each ion type
    a_ratios = np.zeros(len(protein_seq))
    b_ratios = np.zeros(len(protein_seq))
    c_ratios = np.zeros(len(protein_seq))
    x_ratios = np.zeros(len(protein_seq))
    y_ratios = np.zeros(len(protein_seq))
    z_ratios = np.zeros(len(protein_seq))
    for site in site_list:

        site_num = site.get_seq_index(len(protein_seq))

        if site.term == 'C':
            x_int = 0
            y_int = 0
            z_int = 0

            for hit in site.hits:
                # to no considering neutral losses don't index at 0
                if hit.thy_ion.iontype[0] == 'x':
                    x_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype[0] == 'y':
                    y_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype[0] == 'z':
                    z_int += float(hit.exp_ion.pkar_cluster)
                else:
                    print('INVALID ION-TYPE! Ion-type was {}'.format(hit.thy_ion.iontype))
            if not norm:
                x_ratio = x_int
                x_ratios[site_num - 1] = x_ratio
                y_ratio = y_int
                y_ratios[site_num - 1] = y_ratio
                z_ratio = z_int
                z_ratios[site_num - 1] = z_ratio
            else:
                # get total intensity at site per ion type and normalize it against max found in any site
                x_ratio = float(x_int) / max_int
                x_ratios[site_num - 1] = x_ratio
                y_ratio = float(y_int) / max_int
                y_ratios[site_num - 1] = y_ratio
                z_ratio = float(z_int) / max_int
                z_ratios[site_num - 1] = z_ratio

        elif site.term == 'N':
            a_int = 0
            b_int = 0
            c_int = 0
            # check for a and b fragments and compare cluster total areas
            for hit in site.hits:
                if hit.thy_ion.iontype[0] == 'a':
                    a_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype[0] == 'b':
                    b_int += float(hit.exp_ion.pkar_cluster)
                elif hit.thy_ion.iontype[0] == 'c':
                    c_int += float(hit.exp_ion.pkar_cluster)
                else:
                    print('INVALID ION-TYPE! Ion-type was {}'.format(hit.thy_ion.iontype))
            if not norm:
                a_ratio = a_int
                a_ratios[site_num - 1] = a_ratio
                b_ratio = b_int
                b_ratios[site_num - 1] = b_ratio
                c_ratio = c_int
                c_ratios[site_num - 1] = c_ratio
            else:
                a_ratio = float(a_int) / max_int
                a_ratios[site_num - 1] = a_ratio
                b_ratio = float(b_int) / max_int
                b_ratios[site_num - 1] = b_ratio
                c_ratio = float(c_int) / max_int
                c_ratios[site_num - 1] = c_ratio
        else:
            print('INVALID TERMINAL! Terminal was ' + site.term)


    return a_ratios, b_ratios, c_ratios, x_ratios, y_ratios, z_ratios


def seq_error_analysis(protein_seq, site_list):
    """
    Generate lists of errors (ppm) at each sequence position to plot
    :param protein_seq: protein amino acid sequence string
    :param site_list: list of FragmtSite objects containing hit information
    :return: error list by sequence position (ready for general seq plot/write), error lists not collapsed
    (all errors listed at each site, rather than averaged)
    """
    # compute normalized intensity at each site and in list by site index
    error_list = np.zeros(len(protein_seq))
    #Per-term - the error list not per term populates nterm first and then when the cter is populated it does not add to an index if it already filled with nterm hits
    error_listn = np.zeros(len(protein_seq))
    error_listc = np.zeros(len(protein_seq))
    error_lists_all = [[] for x in np.arange(len(protein_seq))]
    for site in site_list:
        # print(site)
        # get total intensity at this site and normalize it against max found in any site
        total_int = 0
        errors = []

        if site.hits:
            for hit in site.hits:
                total_int += float(hit.cal_error)
                errors.append(hit.cal_error)
            avg_error=total_int / len(site.hits)
        else:
            avg_error=0


        # get terminal and correct so C- and N-terminal ions are using same numbering scheme
        site_num = (site.get_seq_index(len(protein_seq)))

        # print(site_num)
        # print(site.hits)
        # print(errors)
        # print(avg_error)
        # print(len(errors))
        # print(error_list)
        # print(error_lists_all)

        #If only error_list[site_num-1] = avg_error is used the c-term sites overwrites the n-term sites in the array
        if site.term == "N":
            error_list[site_num-1] = avg_error
            error_listn[site_num - 1] = avg_error
        else:
            if error_list[site_num-1] == 0:
                error_list[site_num-1] = avg_error
            error_listc[site_num - 1] = avg_error

        error_lists_all[site_num-1] = errors
        # print(len(error_lists_all))
        # print(len(error_list))
    # print(error_list)
    # print(len(error_list))

    # print(error_lists_all)
    return error_list, error_lists_all


def general_seq_plot(datalists, labels, outputdir, plotitle, xtitle=None, ytitle=None, seqcovnumber=None, errorlists=None, extension='.png'):
    """
    Plot any data as a function of protein sequence position (from N-terminus)
    :param datalists: list of lists of data to be stacked *in order* as a function of seq. MUST all be same length
    :param labels: list of labels for legend
    :param errorlists: list of lists of error bars (optional). Must be in same order as datalists, must all be same len
    :param outputdir: directory in which to save output
    :param xtitle: x axis title (optional)
    :param ytitle: y axis title (optional)
    :param seqcovnumber: int, The percentage of sequence coverage
    :param plotitle: title of plot, required, as file is saved with this name
    :return: none (void)
    """
    plt.clf()
    plt.figure(dpi=300)

    indices = np.arange(0, len(datalists[0]))
    index = 0
    for data in datalists:
        if index == 0:
            if errorlists is not None:
                plt.bar(indices, data, label=labels[index], yerr=errorlists[index])
            else:
                plt.bar(indices, data, label=labels[index])
        else:
            if errorlists is not None:
                plt.bar(indices, data, label=labels[index], yerr=errorlists[index], bottom=datalists[index - 1])
            else:
                # need to sum all previously plotted lists to get correct bottom coordinates
                summedlist = [sum(x) for x in zip(*datalists[0:index])]
                plt.bar(indices, data, label=labels[index], bottom=summedlist)
        index += 1
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    if seqcovnumber is not None:
        plt.title(plotitle + '\n {}%'.format(seqcovnumber))
    else:
        plt.title(plotitle)
    plt.legend(loc='best')

    plotname = plotitle + extension
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()


def seq_plot_nterm(datalists, labels, outputdir, plotitle, xtitle=None, ytitle=None, errorlists=None):
    """
    Plot any data as a function of protein sequence position (from N-terminus)
    :param datalists: list of lists of data to be stacked *in order* as a function of seq. MUST all be same length
    :param labels: list of labels for legend
    :param errorlists: list of lists of error bars (optional). Must be in same order as datalists, must all be same len
    :param outputdir: directory in which to save output
    :param xtitle: x axis title (optional)
    :param ytitle: y axis title (optional)
    :param plotitle: title of plot, required, as file is saved with this name
    :return: none (void)
    """
    plt.clf()
    plt.figure(dpi=200)

    indices = np.arange(0, len(datalists[0]))
    index = 0
    for data in datalists:
        if index == 0:
            if errorlists is not None:
                plt.bar(indices[0:len(indices)//2], data[0:len(data)//2], label=labels[index], yerr=errorlists[index])
            else:
                plt.bar(indices[0:len(indices)//2], data[0:len(data)//2], label=labels[index])
        else:
            if errorlists is not None:
                plt.bar(indices[0:len(indices)//2], data[0:len(data)//2], label=labels[index], yerr=errorlists[index], bottom=datalists[index - 1])
            else:
                # need to sum all previously plotted lists to get correct bottom coordinates
                summedlist = [sum(x) for x in zip(*datalists[0:index])]
                plt.bar(indices[0:len(indices)//2], data[0:len(data)//2], label=labels[index], bottom=summedlist)
        index += 1

    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(plotitle)
    plt.legend(loc='best')

    plotname = plotitle + 'Nterm' + '.png'
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()

def seq_plot_cterm(datalists, labels, outputdir, plotitle, xtitle=None, ytitle=None, errorlists=None):
    """
    Plot any data as a function of protein sequence position (from N-terminus)
    :param datalists: list of lists of data to be stacked *in order* as a function of seq. MUST all be same length
    :param labels: list of labels for legend
    :param errorlists: list of lists of error bars (optional). Must be in same order as datalists, must all be same len
    :param outputdir: directory in which to save output
    :param xtitle: x axis title (optional)
    :param ytitle: y axis title (optional)
    :param plotitle: title of plot, required, as file is saved with this name
    :return: none (void)
    """
    plt.clf()
    plt.figure(dpi=200)

    indices = np.arange(0, len(datalists[0]))
    index = 0
    for data in datalists:
        if index == 0:
            if errorlists is not None:
                plt.bar(indices[len(indices)//2:], data[len(data)//2:], color='tab:orange', label=labels[index], yerr=errorlists[index])
            else:
                plt.bar(indices[len(indices)//2:], data[len(data)//2:],  color='tab:orange',label=labels[index])
        else:
            if errorlists is not None:
                plt.bar(indices[len(indices)//2:], data[len(data)//2:],  color='tab:orange',label=labels[index], yerr=errorlists[index], bottom=datalists[index - 1])
            else:
                # need to sum all previously plotted lists to get correct bottom coordinates
                summedlist = [sum(x) for x in zip(*datalists[0:index])]
                plt.bar(indices[len(indices)//2:], data[len(data)//2:],  color='tab:orange',label=labels[index], bottom=summedlist)
        index += 1

    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(plotitle)
    plt.legend(loc='best')

    plotname = plotitle + 'Cterm' + '.png'
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()

def general_seq_write(datalists, labels, seq=None, title='New analysis:', prev_file=None, new_filename=None):
    """
    Method to append sequence information to an existing file, or write a new file with sequence information.
    **One of prev_file OR output_dir must be provided** (if both, defaults to prev_file)
    :param datalists: list of lists of data to be written to file as columns in output table
    :param labels: list of column labels
    :param seq: (optional) protein sequence (string). Write one protein seq letter for each row if provided
    :param title: title to place at top of table
    :param prev_file: file (full path) to append to
    :param new_filename: ONLY needed if prev_file is None. Full path to save if writing new file
    :return: none (void)
    """

    if prev_file is not None:
        myfile = prev_file
        mode = 'a'
    elif new_filename is not None:
        myfile = new_filename
        mode = 'w'
    else:
        print('No file path provided for seq_write. Data will not be saved')
        return

    with open(myfile, mode) as outfile:
        # write title and headers
        outfile.write('\n' + title + '\n')
        label_line = ''
        if seq is not None:
            label_line = 'Seq Position,AA,'

        # for label in labels:
        label_args = ['{}'.format(x) for x in labels]
        label_line = label_line + ','.join(label_args)
        outfile.write(label_line + '\n')

        # write data
        test = list(zip(*datalists))
        # index = 1
        index = 0
        for data in test:
            data_args = ['{}'.format(x) for x in data]
            line = ','.join(data_args)
            if seq is not None:
                if index == 0:
                    outfile.write(str(index) + ',-, ' + line + '\n')
                else:
                    outfile.write(str(index) + ',' + seq[index - 1] + ', ' + line + '\n')
            else:
                outfile.write(line + '\n')
            index += 1


def error_scatter_plot(list_of_error_lists, outputdir, plotitle, xtitle=None, ytitle=None):
    """
    Plot error (ppm) as a function of sequence position in a scatter plot
    :param list_of_error_lists: list of lists, length is the length of the protein sequence as each sub-list contains
    all error values for ions from that position
    :param outputdir: save directory
    :param plotitle: plot title
    :param xtitle: (optional) x axis title
    :param ytitle: (optional) y axis title
    :return: void
    """
    plt.clf()

    indices = np.arange(0, len(list_of_error_lists))
    for x_val, y_vals in zip(indices, list_of_error_lists):
        plt.scatter(x_val, y_vals)

    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(plotitle)
    plt.legend(loc='upper right')

    plotname = plotitle + 'errors.png'
    plotfilename = os.path.join(outputdir, plotname)
    plt.savefig(plotfilename)
    plt.close()

def seq_cov_number(protein_seq, site_list, iontype_exclusion=False):
    """
    :param protein_seq: protein amino acid sequence string
    :param site_list: list of FragmtSite objects containing hit information
    :return: percent of sequence coverage, number of fragsites/number of peptidebonds*100. It automatically filters extra ion types
    """
    bad_iontypes = ["z+1", "z+2", "z+3", "z-H2O", "x-H2O", "c-H2O", "z-NH3", "c-NH3", "x-NH3"]
    # initialize counter of FragmentSites with hits
    hitcount = 0
    #Creating list to store the seq_index of each hit
    seqindexListC = []
    seqindexListN= []
    for site in site_list:
        # For n-terminal fragments: If the Fragsite has at least one hit, it is added to the hit counter, also its
        #seq_index is added so that the hits from the c-term that were already covered by a n-term-fragment are not added
        if site.term == 'N' and len(site.hits) != 0:

            if site.hits == 1:
                print(site.hits)
                for x in site.hits:
                    print(f"x cal_error = {x.cal_error}")
                    if x.thy_ion.iontype in bad_iontypes:
                        continue
                    else:
                        seqindexListN.append(site.seq_index)
                        hitcount += 1
            else:
                print(f"Not == 1 {site.hits}")

                seqindexListN.append(site.seq_index)
                hitcount += 1


        if site.term == "C" and len(site.hits) != 0:
            if site.hits == 1:
                print(site.hits)
                for x in site.hits:
                    print(f"x cal_error = {x.cal_error}")
                    if x.thy_ion.iontype in bad_iontypes:
                        continue
                    else:
                        seqindexListC.append(site.seq_index)
                        hitcount += 1
            else:
                print(f"Not == 1 {site.hits}")
                seqindexListC.append(site.seq_index)
                hitcount += 1
            #Removing if the hit correspond to an index already covered by n-term fragments
            if site.seq_index in seqindexListN:
                hitcount -=1

    # sequence coverage percentage is calculated
    # Number of Fragsites over number of peptide bonds in a protein
    seqcov = (hitcount / (len(protein_seq) - 1))*100
    return round(seqcov)


def get_max_int(sitelist):
    """
    computes the maximum intensity found for all hits at a particular site
    :param sitelist: list of fragment sites with hit information initialized
    :return: float - maximum intensity found in the *sum* of all hits for a given site
    """
    max_int = 0
    for site in sitelist:
        site_total = 0
        for hit in site.hits:
            try:
                site_total += float(hit.exp_ion.pkar_cluster)
            except AttributeError:
                # mmass data - use peak height
                site_total += float(hit.exp_ion.pkht_cluster)
        if site_total > max_int:
            max_int = site_total
    return max_int


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
        main_seq_cov(files)


def frips_main_seq_cov(hits_files, seqcov=False, seqcov_iontype=False, overlay=False, combination=False):
    """
    Function created specifically for my FRIPS datasets sequence coverage analysis - CRR
    :param hits_files: (list of strings) full system paths to .hits files to analyze
    :param seqcov: Boolean, only sequence coverage analysis by terminus
    :param seqcov_iontype: Boolean, only sequence coverage analysis by ion type
    :param overlay: Boolean, sequence coverage comparison of two datasets by terminus
    :param combination: Boolean, sequence coverage combination of two datasets
    :return:
    """

    # Save output where the user wants
    outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(outdir)

    if seqcov:
        main_seq_cov(hits_files, outputdir = outdir, extension=".pdf",  iontype_exclusion=False, cz_iontypes_only=False, errors=True, outfile=True, outfile_extension=".csv")
        print("Done!")
    else:

        if seqcov_iontype:
            main_seqcov_iontype(hits_files, outdir, neutloss=True)
        elif overlay:
            main_seq_cov_overlayn(hits_files, outdir)
            main_seq_cov_overlayc(hits_files, outdir)
        elif combination:

            main_seq_cov_combo(hits_files, outdir)


if __name__ == '__main__':
    root = tkinter.Tk()
    root.withdraw()

    # Seq_cov functions
    hitsfiles = filedialog.askopenfilenames(title='Load Hits Files', filetypes=[('Hits', '.hits')])

    # main_seq_cov(hitsfiles)

    # seq_cov_overlay(hits_files, norm_bool=False)

    # batch_main_seq_cov()

    frips_main_seq_cov(hitsfiles, seqcov=False, seqcov_iontype=False, combination=True)

    #TODO: start including internal fragments for sequence coverage calculations



    #Reducing FDR files
    # average_main(batch=True)
    # compare_main()
    #
    #For complicated searches where passes must be separated
    # merging_batch()
    


