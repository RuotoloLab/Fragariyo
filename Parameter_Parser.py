"""
Author: Most code from DP
Date: Feb 19, 2020
"""
from pyteomics import mass
from pyteomics import parser
import pandas as pd
from tkinter import filedialog
from tkinter import simpledialog
import tkinter as tk
import os
# import pyperclip
import pickle
import combination
import re
from tkinter import messagebox
from PyQt5 import QtWidgets
import tkinter
from tkinter import filedialog
from tkinter import messagebox
import os
import sys
import RenameIMTBXoutputs

# parameter position in template file dictionary for quick reference. Key = parameter attribute name, value = position in splits
ppos = {'Analysis Num':0,
        'Analysis Name': 1,
        'seq': 2,
        'iontypes': 3,
        'mincharge':4,
        'maxcharge': 5,
        'min_ifrag_len': 6,
        'max_ifrag_len': 7,
        'noncys_mods': 8,
        'mods_array': 9,
        'r': 10,
        'disulfides': 11,
        'uniprot_offset': 12,
        'ss_allowbroken': 13,
        'disulfides_ls': 14,
        'naturally_redcys':15
        }


INTERIONS = ['c-z', 'c-zdot', 'c-y', 'cdot-y', 'a-z', 'a-zdot', 'a-y', 'b-y', 'x-c', 'x-cdot']

def parse_disulf_ls(disulf_str):
    """
    :param disulf_str: a list of strings of cysteine location involve in disulfide bonds (e.g ['15-18', ''])
    :return: A list of cystine locations that are involve in disulfide bonds (set), e.g. [{18, 15}, set()]
    """

    spl = disulf_str.split(";")

    # print(spl)
    ssls = []
    for ssbond in spl:
        # print(ssbond)
        enlace = ssbond.split("-")
        bondset = set()
        for num in enlace:
            #If statement helps to get rid off and empty set
            if num:
                # print(f"num = {num}")
                # print(f"num type = {type(num)}")
                intnum = int(num)
                bondset.add(intnum)
                # print(f"bondset = {bondset}")
            # else:
            #     print("Not a number!")

        ssls.append(bondset)

    # print(f"ssls = {ssls}")
    return ssls

def parse_param_template(param_file, terminal = None):
    """
    Read template csv file for all parameter information.
    :param param_file: (string) full system path to csv template file
    """
    with open(param_file, 'r') as pfile:
        for line in list(pfile):
            if line.startswith('#'):
                continue

            splits = line.rstrip('\n').split(',')
            # print(splits)

            #Initilize params object
            params = Parameters()
            params.analysisName = splits[ppos['Analysis Name']]
            params.seq = splits[ppos['seq']].strip()
            params.mincharge = int(splits[ppos['mincharge']])
            params.maxcharge = int(splits[ppos['maxcharge']])

            #Ion types
            iontypes_str = splits[ppos['iontypes']]
            iontypes_ls = []
            # print(iontypes_str)
            iontypes_strsplit = iontypes_str.split(';')
            for type in iontypes_strsplit:
                # print(type)
                # for x in IONS:
                #     if x == type:
                #         iontypes_ls.append(x)
                #     else:s
                #         print("Incorrect Ion type!")
                if type in INTERIONS:
                        iontypes_ls.append(type)

            # print(iontypes_ls)
            params.iontypes = iontypes_ls

            #Internal fragment length
            params.min_len = int(splits[ppos['min_ifrag_len']])
            params.max_len = int(splits[ppos['max_ifrag_len']])

            # Parameteres for modification permutations for disulfide breakage
            modstr = splits[ppos['mods_array']]
            params.arr = mods_fromstr_tols(modstr)

            params.r = splits[ppos['r']]

            #Disulfide_analysis
            params.disulfide_bool = parse_bool(splits[ppos['disulfides']])
            params.uniprot_offset = int(splits[ppos['uniprot_offset']])
            params.ss_allowbroken = int(splits[ppos['ss_allowbroken']])


            #Parsing disulfides
            disulfides_str = splits[ppos['disulfides_ls']]
            params.disulfide_ls = parse_disulf_ls(disulfides_str)


    return params


def mods_fromstr_tols(modstr):
    """
    :param modstr: A string of modifications
    :return:A lsit of modifications
    """
    # Parameteres for modification permutations for disulfide breakage

    modssplit = modstr.split(';')
    mods_ls = []
    for mod in modssplit:
        # print(f'the mod is {mod} and mod.isalum():{mod.isalnum()}')
        if mod.isalnum():
            mods_ls.append((mod))

    # print(f"In the Parser this is modbool! {mods_ls} and its length is {len(mods_ls)}")
    return mods_ls

def parse_param_template_batch(param_file):
    """
    Read template csv file for all parameter information.
    :param param_file: (string) full system path to csv template file
    """

    params_dict = {}

    with open(param_file, 'r') as pfile:
        for line in list(pfile):
            if line.startswith('#'):
                continue

            splits = line.rstrip('\n').split(',')
            # print(splits)

            #Initilize params object
            params = Parameters()
            params.analysisName = splits[ppos['Analysis Name']]
            params.seq = splits[ppos['seq']].strip()
            params.mincharge = int(splits[ppos['mincharge']])
            params.maxcharge = int(splits[ppos['maxcharge']])

            #Ion types
            iontypes_str = splits[ppos['iontypes']]
            iontypes_ls = []
            # print(iontypes_str)
            iontypes_strsplit = iontypes_str.split(';')
            for type in iontypes_strsplit:
                print(type)
                if type in INTERIONS:
                    iontypes_ls.append(type)
                else:
                    print("Incorrect Ion type!")
            # print(iontypes_ls)
            params.iontypes = iontypes_ls

            #Internal fragment length
            params.min_len = int(splits[ppos['min_ifrag_len']])
            params.max_len = int(splits[ppos['max_ifrag_len']])

            #Parameteres for modification permutations
            modstr = splits[ppos['mods_array']]
            modssplit = modstr.split(';')
            mods_ls = []
            for mod in modssplit:
                mods_ls.append((mod))
            params.arr = mods_ls
            params.r = splits[ppos['r']]

            #Disulfide_analysis
            params.disulfide_bool = parse_bool(splits[ppos['disulfides']])
            params.uniprot_offset = int(splits[ppos['uniprot_offset']])
            params.ss_allowbroken = int(splits[ppos['ss_allowbroken']])


            #Parsing disulfides
            disulfides_str = splits[ppos['disulfides_ls']]
            params.disulfide_ls = parse_disulf_ls(disulfides_str)

            params_dict[params.analysisName] = params


    return params_dict

def parse_param_template_batch_multipass(param_file, terminal=None):
    """
    Read template csv file for all parameter information.
    :param param_file: (string) full system path to csv template file
    """

    params_dict = {}

    with open(param_file, 'r') as pfile:
        processed_analysis = 0
        for line in list(pfile):
            # print(line)

            # print(f"Processed: {processed_analysis}")
            if line.startswith('#'):
                continue

            splits = line.rstrip('\n').split(',')
            # print(splits)
            # print(len(splits))


            current_analysis = splits[ppos['Analysis Num']]
            # print(f"Current: {current_analysis}")
            if current_analysis != processed_analysis:
                params_dict[current_analysis] = []
            #Initilize params object
            params = Parameters()
            params.analysisName = splits[ppos['Analysis Name']]
            params.analysisNum = splits[ppos['Analysis Num']]
            params.seq = splits[ppos['seq']].strip()
            params.mincharge = int(splits[ppos['mincharge']])
            params.maxcharge = int(splits[ppos['maxcharge']])

            #Ion types
            iontypes_str = splits[ppos['iontypes']]
            iontypes_ls = []
            # print(iontypes_str)
            iontypes_strsplit = iontypes_str.split(';')
            for type in iontypes_strsplit:
                # print(type)
                if type in INTERIONS:
                    iontypes_ls.append(type)
                else:
                    print("Incorrect Ion type!")
            # print(iontypes_ls)
            params.iontypes = iontypes_ls

            #Internal fragment length
            params.min_len = int(splits[ppos['min_ifrag_len']])
            params.max_len = int(splits[ppos['max_ifrag_len']])

            #Parameteres for modification permutations
            # Parameteres for modification permutations for disulfide breakage
            modstr = splits[ppos['mods_array']]
            params.arr = mods_fromstr_tols(modstr)
            # print(f"params.arr= {params.arr}")
            modstr = splits[ppos['noncys_mods']]
            params.noncysmods = mods_fromstr_tols(modstr)
            params.r = splits[ppos['r']]

            #Disulfide_analysis
            params.disulfide_bool = parse_bool(splits[ppos['disulfides']])
            params.uniprot_offset = int(splits[ppos['uniprot_offset']])
            params.ss_allowbroken = int(splits[ppos['ss_allowbroken']])
            # print(f"ppos['naturally_redcys']= {ppos['naturally_redcys']}")
            # print(f"splits[ppos['naturally_redcys']]= {splits[ppos['naturally_redcys']]}")
            natredcysstr = splits[ppos['naturally_redcys']]
            params.naturally_redcys = str_to_ls(natredcysstr)


            #Parsing disulfides
            disulfides_str = splits[ppos['disulfides_ls']]
            params.disulfide_ls = parse_disulf_ls(disulfides_str)
            params_dict[params.analysisNum].append(params)

            processed_analysis = current_analysis
            # print(f"after processed: {processed_analysis}")

    return params_dict

def str_to_ls(modstr):
    """
    :param modstr:
    :return:
    """
    # Parameteres for modification permutations for disulfide breakage

    modssplit = modstr.split(';')
    mods_ls = []
    for mod in modssplit:
        mods_ls.append((mod))
    return mods_ls

def parse_bool(param_string):
    """
    Parse input strings to boolean
    :param param_string: string
    :return: bool
    """
    if param_string.lower() in ['t', 'true', 'yes', 'y']:
        return True
    elif param_string.lower() in ['f', 'false', 'no', 'n']:
        return False
    else:
        raise ValueError('Invalid boolean: {}'.format(param_string))

class Parameters(object):
    """
    Container to hold all parameters for searches to simplify method calls/etc
    """

    def __init__(self):
        """
        No parameters initialized immediately
        """
        self.params_dict = {}

        # ion prediction parameters
        self.analysisName = None
        self.analysisNum = None
        self.seq = None
        self.iontypes = None
        self.mincharge = None
        self.maxcharge = None
        self.noncysmods = None

        # Internal Fragment length
        self.min_len = None
        self.max_len = None

        # Disulfide Analysis
        self.arr = None
        self.r = None
        self.disulfide_bool = None
        self.uniprot_offset = None
        self.ss_allowbroken = None



    def set_params(self, params_dict):
        """
        Set a series of parameters given a dictionary of (parameter name, value) pairs
        :param params_dict: Dictionary, key=param name, value=param value
        :return: void
        """
        for name, value in params_dict.items():
            try:
                # only set the attribute if it is present in the object - otherwise, raise attribute error
                self.__getattribute__(name)
                self.__setattr__(name, value)
            except AttributeError:
                # no such parameter
                print('No parameter name for param: ' + name)
                continue
        self.update_dict()

    def combodict_calc(self):
        combo_dict = combination.batch_combos(self.arr, int(self.r))
        return combo_dict

    def update_dict(self):
        """
        Build (or rebuild) a dictionary of all attributes contained in this object
        :return: void
        """
        for field in vars(self):
            value = self.__getattribute__(field)
            self.params_dict[field] = value

    def __str__(self):
        """
        string
        :return: string
        """
        return '<Params> protein {}'.format(self.analysisName)
    __repr__ = __str__


def isotope_xtractor():
    """
    Function that takes a .csv file for experimental ions picked manually and a .csv file with the xy coordinate fo the mass spectrum.
    A .csv file ready to be input into the internal fragmentor is produced.
    The output contains the neutral mass, the m/z value, the charge and a list or values for the m/z and int dimensions.
    This way isotope envelopes can be compared between experimental and theoretical ions
    """

    # set output directory
    main_outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(main_outdir)

    # Input experimental ions unmatched by the terminal fragmentor
    expions = filedialog.askopenfilenames(title='expions to find', filetypes=[('CSV', '.csv')])
    print(expions)

    for expfile in expions:

        # Extract file name
        pathsplit = expfile.split('/')
        samplename = pathsplit[-1].strip(".csv")
        print(samplename)

        # Initiate dictionary  to store the coordinates for the peaks
        expion_dict = {}

        # Extracting the csv as a pandas DataFrame
        pd_file = pd.read_csv(expfile, header=0, engine='python')
        # pd_mz = pd_file['mz']
        # pd_int = pd_file['z']

        # Extract mz and z from file
        for index in range(len(pd_file)):
            ion_info = pd_file.loc[index]
            print(ion_info[0], ion_info[1])
            expion_dict[ion_info[0]] = []
            expion_dict[ion_info[0]].append(ion_info[1])
            expion_dict[ion_info[0]].append(ion_info[2])

        print(expion_dict)

        # input xy coordinates
        coordinates = filedialog.askopenfilenames(title='XY Coordinates', filetypes=[('CSV', '.csv')])
        print(coordinates)

        for file in coordinates:
            pathsplit = file.split('/')
            # proteiname = pathsplit[-1].strip("_expions.csv")
            coorname = pathsplit[-1]
            # print(coorname)

            # Extracting the csv as a pandas DataFrame
            pd_file = pd.read_csv(file, header=0, engine='python', usecols=['X(MassToCharge)', 'Y(Counts)'])
            # pd_mz = pd_file['X(MassToCharge)']
            # pd_int = pd_file['Y(Counts)']

            # Go an look for the range of mz and int tha covered the envelope of each peak
            for peak in expion_dict:
                round_peak = round(peak, 4)

                # What is the minimun mz value
                rndone_mzrange = pd_file.loc[pd_file['X(MassToCharge)'] > round_peak - 0.5]
                # print(rndone_mzrange)
                # print(type(rndone_mzrange))

                # What is the maximun mz value
                rndtwo_mzrange = rndone_mzrange.loc[rndone_mzrange['X(MassToCharge)'] < round_peak + 4]
                # print(rndtwo_mzrange)

                # extrac columns as numpy arrays
                np_coormz = rndtwo_mzrange['X(MassToCharge)'].to_numpy()
                np_coorint = rndtwo_mzrange['Y(Counts)'].to_numpy()

                # print(np_coormz)
                # print(np_coorint)

                # ";" is a separator so that input .csv file can be used by the internal fragmentor
                mz_output = ";".join(str(x) for x in np_coormz)
                int_output = ";".join(str(x) for x in np_coorint)

                # Add to dictionary
                expion_dict[peak].append((mz_output, int_output))

            # print(expion_dict)

            # parse each peak as a list in a list to create the output DataFrame
            final_df_ls = []
            for val in expion_dict:
                int = expion_dict[val][1]
                z_val = expion_dict[val][0]
                mz_val = val
                # Calculate neutral mass
                neut_val = (mz_val * z_val) - (z_val * 1.0078)
                mz_ls = expion_dict[val][2][0]
                int_ls = expion_dict[val][2][1]
                final_df_ls.append([neut_val, z_val, mz_val, int, mz_ls, int_ls])
            # print(final_df_ls)

            # Create Output DataFrame
            final_df = pd.DataFrame(final_df_ls,
                                    columns=['#neut', 'z', 'mz', 'int','isoenv_mz', 'isoenv_int'])

            print(final_df)

            # Convert FataFrame to .csv
            final_df.to_csv(f"{samplename}_expion_isoenv.csv", index=False)


class expionObj:
    """
    Contianer of experimental ions (Fragmentor)
    :param exp_neut: Neutral mass
    :param exp_mz: Experimental mass
    :param charge: ion charge
    :param mz_isoenv: list with m/z values covering the isotopic envelope of the experimental ion
    :param int_isoenv: list with intensity values covering the isotopic envelope of the experimental ion

    """
    def __init__(self, exp_neut, exp_mz, charge, pkht_cluster, mz_isoenv, int_isoenv):
        self.exp_neut = round(exp_neut, 4)
        self.exp_mz = round(exp_mz,4)
        self.charge = float(charge)
        self.pkht_cluster = pkht_cluster
        self.mz_isoenv = mz_isoenv
        self.int_isoenv = int_isoenv

    def __str__(self):
        return f"m:{self.exp_neut}\tz:{self.charge}"

    __repr__ = __str__



def expion_parser(expfile):
    """
    Function to parse a .csv file with the experimental neutral masses
    :return: a list of experimental ion objects
    """
    expion_ls = []
    pathsplit = expfile.split('/')
    proteiname = pathsplit[-1].strip("_expions.csv")


    with open(expfile, 'r') as batch:
        lines = list(batch)

        for line in lines:

            if line.startswith("#"):
                continue
            else:
                mz_isoenv_ls = []
                int_isoenv_ls = []
                # print(line)
                line = line.strip("\n")
                splits = line.split(',')
                if splits[0]:
                    neutral_mono = float(splits[0])
                    z = splits[1]
                    mz = float(splits[2])
                    int = float(splits[3])

                    mz_isoenv = splits[4]
                    if mz_isoenv:
                        mz_isoenv_spl = mz_isoenv.split(';')
                        for  mz_val in mz_isoenv_spl:

                            mz_isoenv_ls.append(float(mz_val))

                    int_isoenv = splits[5]
                    if int_isoenv:
                        int_isoenv_spl = int_isoenv.split(';')
                        for int_val in int_isoenv_spl:
                            int_isoenv_ls.append(float(int_val))

                # print(mz_isoenv_ls)
                # print(int_isoenv_ls)
                expion = expionObj(neutral_mono, mz, z, int, mz_isoenv_ls,int_isoenv_ls)
                expion_ls.append(expion)



    return expion_ls, proteiname

if __name__ == "__main__":
    pass