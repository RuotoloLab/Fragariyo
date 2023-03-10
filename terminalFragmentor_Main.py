"""
IMAnnotatorv3: PeakMatcher for Terminal Fragments
Author: Carolina Rojas Ramirez
Date: 05/22/2020
In silico fragmentation of proteins using mass.fast_mass2 from pyteomics
Pieces of code from Daniel A. Polasky's IMAnnotatorv2
Pieces of code from Carolina Rojas Ramirez' Fragmentor (internal fragment analysis + Disulfinator)
"""

import time
import combination
from PyQt5 import QtWidgets
import tkinter
import sys
import RenameIMTBXoutputs
from pyteomics import mass
from tkinter import filedialog
import os
import pickle
import Parameter_Parser_terminal
import re
from tkinter import messagebox
# from Modifications import mods_repo
from PeakMatch import matchmaker_terminal_multipass
import multiprocessing

# update Pyteomics masses to use custom ion types
mass.std_ion_comp.update({
    'M':        mass.Composition(formula=''),
    'M-H2O':    mass.Composition(formula='H-2O-1'),
    'M-NH3':    mass.Composition(formula='N-1H-3'),
    'a':        mass.Composition(formula='H-2O-1' + 'C-1O-1'),
    'a+1':        mass.Composition(formula='H-2O-1' + 'C-1O-1'+'H1'),
    'a-H2O':    mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3':    mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b':        mass.Composition(formula='H-2O-1'),
    'b-H2O':    mass.Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3':    mass.Composition(formula='H-2O-1' + 'N-1H-3'),
    'c':        mass.Composition(formula='H-2O-1' + 'NH3'),
    'c-dot':    mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1'),
    'c-1':    mass.Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+1':      mass.Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+2':      mass.Composition(formula='H-2O-1' + 'NH3' + 'H2'),
    'c-H2O':    mass.Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3':    mass.Composition(formula='H-2O-1'),
    'x':        mass.Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O':    mass.Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3':    mass.Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y':        mass.Composition(formula=''),
    'y-H2O':    mass.Composition(formula='H-2O-1'),
    'y-NH3':    mass.Composition(formula='N-1H-3'),
    'z-dot':        mass.Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z':    mass.Composition(formula='H-2O-1' + 'ON-1'),
    'z+1':      mass.Composition(formula='H-2O-1' + 'ON-1H1'),
    'z+2':      mass.Composition(formula='H-2O-1' + 'ON-1H2'),
    'z+3':      mass.Composition(formula='H-2O-1' + 'ON-1H3'),
    'z-H2O':    mass.Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3':    mass.Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
    #Internal fragments ion types
    'c-z': mass.Composition(formula='H-2O-1' + 'NH3'+'H-2O-1' + 'ON-1'),
    'c-zdot': mass.Composition(formula='H-2O-1' + 'NH3'+'H-2O-1' + 'ON-1H-1'),
    'cdot-z': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1'+'H-2O-1' + 'ON-1'),
    'c-y': mass.Composition(formula='H-2O-1' + 'NH3'+ ''),
    'cdot-y': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1'+''),
    'a-z': mass.Composition(formula='H-2O-1' + 'C-1O-1'+'H-2O-1' + 'ON-1'),
    'a-zdot': mass.Composition(formula='H-2O-1' + 'C-1O-1'+'H-2O-1' + 'ON-1H-1')
    })


class FragmentSite:
    """
    Container to hold information about a particular site along the protein backbone.
    """

    def __init__(self, sequence, terminal, protein_seq, startindex, endindex, resi_dict, ss_ls, intrabroken_ss, naturally_reducedcys_ls):
        """
        :param sequence: str, the fragment site sequence (endpoint and start point depends from which terminal the fragment is originating from.
        :param terminal: str, the fragment site origin: N or C.
        :param protein_seq: str, the full protein sequence.
        :param startindex: int, at what position in the protein does the fragment site sequence starts
        :param endindex: int, at what position in the protein does the fragment site sequence ends
        :param resi_dict: dict, keys:amino acids (str); values: set, set of indexes (int) corresponding to the indeces
        at which the amino acid is located in the fragment site sequence
        :param ss_ls: ls, a list of the disulfide bonds (set) in the protein
        :param intrabroken_ss: int, how many disulfide bonds (kept intact inside the fragment site) are meant to be broken
        :param naturally_reducedcys_ls: ls, a list of the location (int) of cysteines  that do not participate in disulfide bonds
        """
        self.seq = sequence
        self.term = terminal
        #Dictionary to store theoretical ions belonging to this FragmentSite
        self.theo_ions = {}
        #List to whihc the matched fragment ions are added
        self.hits = []
        #Obtain sequence index using protein sequence
        self.seq_index = self.get_seq_index(len(protein_seq))
        self.full_protein_seq = protein_seq
        self.startindex = startindex
        self.endindex = endindex
        self.resi_dict = resi_dict
        self.ss_ls = ss_ls
        self.intrabroken_ss = intrabroken_ss
        self.naturally_reducedcys_ls = naturally_reducedcys_ls

    def disulfide_analysis(self):
        """
        :return: unboundcys, cysloc, disulfide_counter
        unboundcys: int, the number of cystines that are naturally reduced or have lost their disulfide counterpart once the fragment site is formed.
        cysloc: set, the indeces at which Cysteines are located in the Fragmentsite and which participate in disulfides
        disulfide_counter: int, how many intact disulfides are in this fragment
        """

        #Get cysteine locations
        cysloc = self.resi_dict["C"]
        # How many cysteines are there? -Total number cysteines based on location list
        cysloclen = len(self.resi_dict["C"])

        #If there are zero cysteines there should not be any unbound cys or disulfide bonds
        if cysloclen == 0:
            unboundcys = 0
            disulfide_counter = 0
        else:

            #remove the naturally reduced cysteines
            #  from the all cysteine locations
            for cys in cysloc:
                if cys in self.naturally_reducedcys_ls:
                    # print(cys)
                    cysloclen -= 1

            disulfide_counter = 0
            # If ss_ls is obtained form uniprot, the uniprot offset is removed in the Parameter_Parser_terminal.py
            for pair in self.ss_ls:
                if pair.issubset(cysloc):
                    # print(disulfide)
                    disulfide_counter += 1

            # intrabroken_ss refers to how many disulfides inside the fragment to be broken
            if self.intrabroken_ss:
                unboundcys = cysloclen - (disulfide_counter * 2) + (self.intrabroken_ss * 2)
                if unboundcys > cysloclen:
                    unboundcys = cysloclen

            else:

                # Cysteines available for modification, usually considering disulfides reaching outside the fragment
                unboundcys = cysloclen - (disulfide_counter * 2)
                if unboundcys > cysloclen:
                    unboundcys = cysloclen

        return unboundcys, cysloc, disulfide_counter


    def get_seq_index(self, protein_len):
        """
        Determines the N-terminal relative position *starting at 0* for a site.
        E.g. a b10 ion (N-term, site #10, starting at 1) would be index 9. A y10 ion in a 100 res protein would be
        index 90 (100 - 10).
        :param protein_len: length of the protein sequence this site is part of (for computing C-term positions)
        :return: int - sequence position of this site from N-term, starting at 0
        """

        if self.term == 'N':
            # site_num = len(self.seq) - 1
            site_num = len(self.seq)
        else:
            site_num = protein_len - len(self.seq) + 1
        return site_num

    #attributes necessary for using FragmentSite as dict keys

    #Attribute to compare if the FragmentsSites are the same based not sequence
    def __eq__(self, other):
        return self.seq == other.seq

    # Attribute to make FragmentSite obj unique items
    def __hash__(self):
        # print(hash(str(self)))
        return hash((self.seq, self.seq_index, self.term))

    #To print
    def __str__(self):
        """
        print string
        :return: string
        """
        return '<Site> {}-{}/theo: {}_hits: {}'.format(self.term, len(self.seq), len(self.theo_ions),len(self.hits))

    # Representation
    __repr__ = __str__



class ThyIon:
    """
    Holds predicted/theoretical ("thy") information about a specific fragment.
    """

    def __init__(self, mz_mono, sequence, charge, ion_type, ion_type_indx, neutlosses, cys_num, mono_neutral, cysmods, thy_mods, cysloc, ss_count, terminal,
                 reverse):
        """
        :param mz_mono = int, m/z of fragment
        :param sequence = str, AA sequence of fragment
        :param charge = int, ion charge
        :param ion_type = str, ion type (e.g. 'a', 'b', 'c', etc...)
        :param ion_type_indx = str, where is the ion type coming from (e.g. '7' for a7, '345' for b345, etcc)
        :param neutlosses = bool, considered neutral losses (H2O or NH4)
        :param cys_num = int, number of cysteines
        :param mono_neutral = int, neutral mass of ion
        :param cysmods = ls, modifications for disulfide bond breakage
        :param thy_mods = ls, list of modifications (str) in the theoretical ions. for variable modifications the str is composed of {modname}x(number pf mods in the theoretical ion).
        :param cysloc: set, set of indeces (int) where the cys residues are located in the fragment site
        :param ss_count = int, how many intact disulfides are in this theoretical fragment
        :param terminal = str, 'N' or 'C'
        :param reverse = bool, is this a reverse sequence fragment site?
        """
        self.mz_mono = round(mz_mono, 8)
        self.charge = charge
        self.ion_type = ion_type
        self.ion_type_indx = ion_type_indx
        self.neutlosses = neutlosses
        self.sequence = sequence
        self.cys_num = cys_num
        self.mono_neutral = round(mono_neutral, 8)
        self.thy_mods = thy_mods
        self.cysmods = cysmods
        self.cysloc = cysloc
        self.ss_count = ss_count
        self.terminal = terminal
        self.reverse = reverse

    # Added so that a list of objects can be sorted by Thy Ion
    def __lt__(self, other):
        return self.mz_mono < other.mz_mono

    # Added so that hits can be functional dictionary keys
    def __eq__(self, other):
        return self.mz_mono == other.mz_mono

    # Added so that objects are hashable
    def __hash__(self):
        # gave issues when averaging, but good when comparing
        #return hash((self.mz, self.iontype, self.charge, self.mods))

        #Good for both averaging and comparing
        return hash((self.mz_mono, self.ion_type, self.charge))

    # Attribute to represent objects in stdout
    def __str__(self):
        return f"{self.mz_mono}\t{self.sequence}\t{self.charge}\t{self.mono_neutral}\t{self.thy_mods}\t{self.cysmods}\t{self.ion_type}\t{self.cysloc}\t{self.ss_count}\t{self.cys_num}\t{self.terminal}\t{self.reverse}"

    # Attribute to represent objects in other objects
    __repr__ = __str__

def varmod_combos(dict_of_varmods):
    """
    :param dict_of_varmods: A dictionary of the variable modifications to analyze
    :return: A dictionary of the variable modifications to analyze with the number of modifications
    """
    varmodls = []
    #how mnay mods to add?
    for varmod in dict_of_varmods:
        for varmod_num in dict_of_varmods[varmod]:
            varmodls.append(f"{varmod}_{varmod_num}")

    # print(f"varmodls = {varmodls}")

    #What are the total amount of mods in the fragment
    totvarmods = len(dict_of_varmods)

    #How can they be combined?
    varmodcombos = combination.rSubset(varmodls, totvarmods)
    # print(varmodcombos)
    return varmodcombos

def varmods_processing(FragmentSiteObj, modification, amino, var_mods_dict, mods_repo):
    """
    Adding variable modifications
    :param FragmentSiteObj: Fragment Site container
    :param modification: Modification object
    :param amino: Residue to be modified
    :param var_mods_dict: Dictionary with the modification and their amount
    :return: A dictionary of variable modification tailored to the current fragment site
    """
    # Variable mod part
    # print(f"before adding varmods = {mods_repo[modification].current_num}")

    #Current number of modifications on the ion
    varmods_num = mods_repo[modification].current_num

    #Let's place the modifications
    if len(FragmentSiteObj.resi_dict[amino]) >= mods_repo[modification].max_num:

        #Based on the maximun number of modifications allowed
        for num_mod in range(1, mods_repo[modification].max_num + 1):

            # print(f"Starting at each num_mod = {mods_repo[modification].current_num}")

            # print(f"put {num_mod} times = m/z {(mods_repo[modification].mass) * num_mod}")
            varmods_num += num_mod

            # print(f"mod current num = {varmods_num}")

            #Don't add more than the allowed per modification
            if varmods_num > mods_repo[modification].max_num:
                continue
            else:
                # var_mods_dict[f"{modification}[{num_mod}]"] = (mods_repo[modification].mass) * num_mod

                var_mods_dict[f"{modification}"].append(num_mod)
            varmods_num -= num_mod

        # mods.remove(modification)

    else:
        # print(f"The places for var_mods {FragmentSiteObj.resi_dict[amino]} is less than the allowed mod max")

        # for num_mod in range(1, len(FragmentSiteObj.resi_dict[amino]) + 1):
        for num_mod in range(1, len(FragmentSiteObj.resi_dict[amino]) + 1):
            # print(f"Available aminos = {len(FragmentSiteObj.resi_dict[amino])}")

            # print(f"put {num_mod} times = m/z {(mods_repo[modification].mass) * num_mod}")
            varmods_num += num_mod
            # print(f"mod current num = {varmods_num}")

            # Don't add more modifications than there are residues
            if varmods_num > len(FragmentSiteObj.resi_dict[amino]):
                continue
            else:
                # var_mods_dict[f"{modification}[{num_mod}]"] = (mods_repo[
                #                                                   modification].mass) * num_mod
                var_mods_dict[f"{modification}"].append(num_mod)
            varmods_num  -= num_mod
        # mods.remove(modification)

    #Reset mods for each fragment
    mods_repo[modification].current_num = 0
    return var_mods_dict

def modificator(FragmentSiteObj, mod_ls, charge, var_mods_dict, mz_mono, neutral_mono, mods_repo):
    """
    Places modifications on a Fragment Site object using all the information necessary
    :param FragmentSiteObj: container for the Fragment Site info
    :param mod_ls: A list of modifications to be considered
    :param charge: Charge state of the theoretical ion to be calculated
    :param var_mods_dict: Dictionary of the modification to be palced with their respective amounts
    :param mz_mono: float, m/z value
    :param neutral_mono: float, deconvoluted mass
    :return: mz_mono, neutral_mono, mods, var_mods_dict for a theoretical ion
    """
    # If modifications need to be considered

    # mods = mod_ls.copy()

    #Intitialized list and masses to be calculated for the theoretical ion
    mods = []
    accum_mz_fx = 0
    accum_neutral_fx = 0

    for modification in mod_ls:
        mods_repo[modification].current_num = 0
        # print(f"the modObj is {mods_repo[modification]} with type {type(mods_repo[modification])}")
        # print(f"the FragSite has the an amino acid composition of {FragmentSiteObj.resi_dict}")
        # print(f"the modification {mods_repo[modification]} can be placed on {mods_repo[modification].target_aas}")

        #Place mod if amino acid is there
        if mods_repo[modification].target_aas:
            for amino in mods_repo[modification].target_aas:
                # print(f"In the  fragment, {amino} is at positions {FragmentSiteObj.resi_dict[amino]}")

                if FragmentSiteObj.resi_dict[amino]:
                    if mods_repo[modification].fixed:
                        for res in mods_repo[modification].fixed:
                            if res in FragmentSiteObj.resi_dict[amino]:
                                # print(
                                #     f"The modification is fixed at  {mods_repo[modification].fixed} which is in this fragment!")
                                accum_mz_fx += (mods_repo[modification].mass) / charge
                                accum_neutral_fx += mods_repo[modification].mass
                                mods.append(modification)
                            # else:
                            #     # print(
                            #     #     f"The modification is fixed at  {mods_repo[modification].fixed} which is NOT in this fragment!")
                            #     mods.remove(modification)


                    elif mods_repo[modification].terminal_flag:

                        if mods_repo[modification].terminal_flag == FragmentSiteObj.term:
                            # mods_repo[modification].current_num = 0
                            # print(
                            #     f"The modification is variable with a max number of {mods_repo[modification].max_num} and it can also be terminal at the {mods_repo[modification].terminal_flag} term, so here! ")

                            mods_repo[modification].current_num += 1
                            accum_mz_fx += (mods_repo[modification].mass) / charge
                            accum_neutral_fx += mods_repo[modification].mass
                            mods.append(f"{mods_repo[modification].terminal_flag}-terminal_{mods_repo[modification].name}")

                            var_mods_dict[f"{modification}"] = []
                            var_mods_dict = varmods_processing(FragmentSiteObj, modification, amino, var_mods_dict, mods_repo)

                        else:
                            # print(
                            #     f"The modification is variable with a max number of {mods_repo[modification].max_num} and it can also be terminal at the {mods_repo[modification].terminal_flag} term, not this term")

                            var_mods_dict[f"{modification}"] = []
                            var_mods_dict = varmods_processing(FragmentSiteObj, modification, amino,
                                                                     var_mods_dict, mods_repo)
                            # mods.remove(modification)

                    else:
                        # print(f"The modification is variable with a max number of {mods_repo[modification].max_num}")
                        var_mods_dict[f"{modification}"] = []
                        var_mods_dict = varmods_processing(FragmentSiteObj, modification, amino, var_mods_dict,mods_repo)


                else:
                    # print(f"There is no {amino} to put it at!")
                    # mods.remove(modification)
                    pass

        else:
            # print(f"{modification} is only terminal!")
            if mods_repo[modification].terminal_flag == FragmentSiteObj.term:
                mods.append(
                    f"{mods_repo[modification].terminal_flag}-terminal_{mods_repo[modification].name}")
                mods_repo[modification].current_num += 1
                accum_mz_fx += (mods_repo[modification].mass) / charge
                accum_neutral_fx += mods_repo[modification].mass

            else:
                # mods.remove(modification)
                # print("This terminal modification does not belong in this terminal!")
                pass

    mz_mono += accum_mz_fx
    neutral_mono += accum_neutral_fx

    return mz_mono, neutral_mono, mods, var_mods_dict


def mass_calc(frag_counter, FragmentSiteObj, ion_types, maxcharge, neutrals = None, ss_bonds=None, cysmodmass_dict = None, cys_modls = None, cys_num=None, cysloc = None, sscount = None, mod_ls =None, reverse_flag=None, modificatorrepo = None):
    """
    Function to calculate the mass of theoretical ions based on the current fragment site and pass parameters
    :param frag_counter: int, to keep track of how many theoretical ions are being produced
    :param FragmentSiteObj: obj, container for the Fragment site that will be used to calculate theoretical ions
    :param ion_types: ls, list of ion types (str) (this are keys to the Pyteomics mass.std_ion_comp dict)
    :param maxcharge: int, the maximun charge to considered to calculate theoretical m/z values
    :param neutrals: str, either H2O or NH3
    :param ss_bonds: bool, will disulfide bonds be considered?
    :param cysmodmass_dict: dict, dictionary with number of disulfides available for modification (keys), possible modifications (values/dict)
    The latter is a dictionary of the possible modifications combinations (keys) and their total mass shift (values)
    :param cys_modls: ls, the losses or gains to be considered when disulfide bonds are broken by CID or ECD. They will
     be used to calculate the possible modifications (see cysmodmass_dict description)
    :param cys_num: int, number of cys residues available for modification
    :param cysloc: set, set of indeces (int) where the cys residues are located in the fragment site
    :param sscount: int, number of disulfide bonds
    :param mod_ls: ls, list of non-cysteine modifications (str, keys to the mods_repo dictionary)
    :param reverse_flag: bool, is this a reverse sequence fragment site?
    :return: updated frag_counter
    """

    # print(f"{FragmentSiteObj}")
    # print(f"{FragmentSiteObj.seq}")

    #For each ion type
    for ion_type in ion_types:
        #For each charge
        for charge in range(1, maxcharge + 1):

            #Use correct ion type based from which terminal the fragment site is coming from
            if FragmentSiteObj.term == 'N':
                if ion_type[0] in ['x', 'y', 'z']:
                    continue
            else:
                if ion_type[0] in ['a', 'b', 'c']:
                    continue

            # calculate m/z and neutral
            mz_mono_org = mass.fast_mass2(FragmentSiteObj.seq, ion_type=ion_type, charge=charge)
            neutral_mono_org = mass.fast_mass2(FragmentSiteObj.seq, ion_type=ion_type, charge=0)

            # print(f"\nOriginal mz = {mz_mono_org}")
            # print(f"Original mass = {neutral_mono_org}")

            #Create a new m/z and mass variable to contain the original mz and neutral mass
            mz_mono = mz_mono_org
            neutral_mono = neutral_mono_org

            #To contain the end results of going thru all the conditions below
            mods = []
            var_mods_dict = {}
            # If disulfides are considered
            if ss_bonds:

                # print(f"ss_bonds = {ss_bonds}")
                # print(f"cys_num = {cys_num}")

                # If the sequence has not free cysteines for modification
                if cys_num == 0:
                    #if besides considering disulfides, possible modifications need to be considered

                    if mod_ls:
                        mz_mono, neutral_mono, mods, var_mods_dict = modificator(FragmentSiteObj, mod_ls, charge,
                                                                                 var_mods_dict, mz_mono, neutral_mono, modificatorrepo)


                        # Removing the hydrogens due to the disulfide bonds
                        neutral_mono = neutral_mono + sscount * (-1.0078 * 2)
                        mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))


                    else:
                        # Do not modified cysteines
                        mods = []

                        # Removing the hydrogens due to the disulfide bonds
                        neutral_mono = neutral_mono + sscount * (-1.0078 * 2)
                        mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))

                #If there are disulfides to be modified because their disulfide partner is not in this fragment
                else:
                    #If modidications need to be considered
                    if mod_ls:
                        mz_mono, neutral_mono, mods, var_mods_dict = modificator(FragmentSiteObj, mod_ls, charge,
                                                                                 var_mods_dict, mz_mono, neutral_mono, modificatorrepo)

                        # except removing the hydrogens due to the disulfide bonds
                        neutral_mono = neutral_mono + sscount * (-1.0078 * 2)
                        mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))

                        # Modify neutral mass with the modifications possible
                        neutral_mono = neutral_mono + cysmodmass_dict[cys_num][cys_modls] + sscount * (-1.0078 * 2)
                        mz_mono = mz_mono + (cysmodmass_dict[cys_num][cys_modls] / charge) + (
                                sscount * ((-1.0078 / charge) * 2))

                    else:
                        # If the internal fragment contains cysteines for modifications

                            # except removing the hydrogens due to the disulfide bonds
                            neutral_mono = neutral_mono + sscount * (-1.0078 * 2)
                            mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))

                            # Modify neutral mass with the modifications possible
                            neutral_mono = neutral_mono + cysmodmass_dict[cys_num][cys_modls] + sscount * (-1.0078 * 2)
                            mz_mono = mz_mono + (cysmodmass_dict[cys_num][cys_modls] / charge) + (
                                        sscount * ((-1.0078 / charge) * 2))

            #If there is no disulfides to be considered, but only modifications
            elif mod_ls:
                mz_mono, neutral_mono, mods, var_mods_dict = modificator(FragmentSiteObj, mod_ls, charge,
                                                                         var_mods_dict, mz_mono, neutral_mono, modificatorrepo)

            # If neutrals are to be considered!
            neutloss = ''
            if neutrals == 'NH3':
                # print("Adding ammonia!")
                mz_mono = mz_mono + (-17.02655/charge)
                neutral_mono = neutral_mono + -17.02655
                neutloss += neutrals

            elif neutrals == 'H2O':
                # print("Adding water!")
                mz_mono = mz_mono  + (-18.01056 / charge)
                neutral_mono = neutral_mono + -18.01056
                neutloss += neutrals


            # Creating theoretical ion and adding it to the all_sites dict
            if FragmentSiteObj.term == 'N':
                #Add correct ion type index based on terminal
                iontype_num = FragmentSiteObj.endindex
                # print(f"For variable mods = {var_mods_dict}")

                #If there were variable modifications to considered
                if var_mods_dict:
                    # print(var_mods_dict)
                    # print(f"There are {len(var_mods_dict)} var mods!")
                    possible_varmodscombos = varmod_combos(var_mods_dict)
                    for combo in possible_varmodscombos:
                        # print(combo)
                        combols = list(combo)
                        mass_combo = 0
                        for mod in combols:
                            modsplits = mod.split("_")
                            # mass_combo += mods_repo[mod[:-2]].mass * int(mod[-1])
                            mass_combo += modificatorrepo[modsplits[0]].mass * int(modsplits[1])
                        # print(f"Mass Combo {mass_combo}")
                        mods.append(combols)
                        mz_mono += mass_combo/charge
                        neutral_mono += mass_combo
                        modscopy = mods.copy()
                        #Create theoretical ion!
                        theoretical_ion = ThyIon(mz_mono, FragmentSiteObj.seq, charge, ion_type, iontype_num, neutloss,
                                                 cys_num, neutral_mono, cys_modls, modscopy, cysloc, sscount,
                                                 FragmentSiteObj.term,
                                                 reverse_flag)
                        # print(f"ThyIon {theoretical_ion}")


                        #Add theoretical ion to its site in the main dictionary of fragment sites
                        FragmentSiteObj.theo_ions[theoretical_ion.mz_mono] = theoretical_ion
                        # print(f"ThyIon in FragSite {FragmentSiteObj.theo_ions[theoretical_ion.mz_mono]}")
                        frag_counter += 1
                        mz_mono -= mass_combo/charge
                        neutral_mono -= mass_combo

                        #Remove variable mod after being done with it
                        mods.remove(combols)

                else:

                    # print(f"Mods before creating ThyIOn {mods}")
                    # print(f"m/z before creating ThyIOn {mz_mono}")
                    # print(f"Mass before creating ThyIOn {neutral_mono}")
                    # print(f"Modifications on Cys {cys_modls}")
                    modscopy = mods.copy()
                    theoretical_ion = ThyIon(mz_mono, FragmentSiteObj.seq, charge, ion_type, iontype_num, neutloss,
                                             cys_num, neutral_mono, cys_modls, modscopy, cysloc, sscount,
                                             FragmentSiteObj.term,
                                             reverse_flag)
                    # print(f"Mods in thy ion {theoretical_ion.mods}")
                    # print(f"cys mods in thy ion {theoretical_ion.cysmods}")
                    # print(f"ThyIon {theoretical_ion}")
                    FragmentSiteObj.theo_ions[theoretical_ion.mz_mono] = theoretical_ion
                    frag_counter += 1




            else:

                # C-terminal fragments
                iontype_num = len(FragmentSiteObj.full_protein_seq) - (FragmentSiteObj.startindex -  1)

                # print(f"For variable mods = {var_mods_dict}")

                if var_mods_dict:
                    # print(var_mods_dict)
                    # print(f"There are {len(var_mods_dict)} var mods!")
                    possible_varmodscombos = varmod_combos(var_mods_dict)
                    for combo in possible_varmodscombos:
                        # print(combo)
                        combols = list(combo)
                        mass_combo = 0
                        for mod in combols:
                            modsplits = mod.split("_")
                            mass_combo += modificatorrepo[modsplits[0]].mass * int(modsplits[1])
                        # print(f"Mass Combo {mass_combo}")
                        mods.append(combols)
                        mz_mono += mass_combo / charge
                        neutral_mono += mass_combo
                        modscopy = mods.copy()
                        # Create theoretical ion!
                        theoretical_ion = ThyIon(mz_mono, FragmentSiteObj.seq, charge, ion_type, iontype_num, neutloss,
                                                 cys_num, neutral_mono, cys_modls, modscopy, cysloc, sscount,
                                                 FragmentSiteObj.term,
                                                 reverse_flag)
                        # print(f"ThyIon {theoretical_ion}")

                        # Add theoretical ion to its site in the main dictionary of fragment sites
                        FragmentSiteObj.theo_ions[theoretical_ion.mz_mono] = theoretical_ion
                        # print(f"ThyIon in FragSite {FragmentSiteObj.theo_ions[theoretical_ion.mz_mono]}")
                        frag_counter += 1
                        mz_mono -= mass_combo / charge
                        neutral_mono -= mass_combo

                        # Remove variable mod after being done with it
                        mods.remove(combols)
                else:

                    # print(f"Mods before creating ThyIOn {mods}")
                    # print(f"m/z before creating ThyIOn {mz_mono}")
                    # print(f"Mass before creating ThyIOn {neutral_mono}")
                    # print(f"Modifications on Cys {cys_modls}")
                    modscopy = mods.copy()
                    theoretical_ion = ThyIon(mz_mono, FragmentSiteObj.seq, charge, ion_type, iontype_num, neutloss,
                                             cys_num, neutral_mono, cys_modls, modscopy, cysloc, sscount,
                                             FragmentSiteObj.term,
                                             reverse_flag)
                    # print(f"Mods in thy ion {theoretical_ion.mods}")
                    # print(f"ThyIon {theoretical_ion}")
                    FragmentSiteObj.theo_ions[theoretical_ion.mz_mono] = theoretical_ion
                    frag_counter += 1

    return frag_counter


def fragments(analysis_name, sequence, types, maxcharge=1, neutral_bool = None, cystine = None, combo_dict = None, ss_ls = None, intrabroken_ss = None, natredcys = None, modbool = None, noncysmods=None,reverse_seq=None, init_tol=None, final_tol=None, cal_bool=None, libraryofmods = None):
    """
    Method to produced terminal theoretical ions for a protein sequence!
    :param analysis_name: str, name to id a pass
    :param sequence: str, full protein sequence
    :param types: ls, list of ion types (str) (this are keys to the Pyteomics mass.std_ion_comp dict)
    :param maxcharge: int, the maximin charge to considered to calculate theoretical m/z values
    :param neutral_bool: bool, 'H2O' and 'NH3' need to be removed?
    :param cystine: bool, does disulfides need to be considered?
    :param combo_dict: dict, dictionary with number of disulfides available for modification (keys), possible modifications (values/dict)
    The latter is a dictionary of the possible modifications combinations (keys) and their total mass shift (values)
    :param ss_ls: ls, a list of the disulfide bonds (set) in the protein
    :param intrabroken_ss: int, how many disulfide bonds (kept intact inside the fragment site) are meant to be broken
    :param natredcys: ls, a list of the location (int) of cysteines  that do not participate in disulfide bonds
    :param modbool: bool, do modifications need to be added?
    :param noncysmods: ls, list of non-cysteine modifications (str, keys to the mods_repo dictionary)
    :param reverse_seq: bool, is this a reverse sequence fragment site?
    :param init_tol: int, error tolerance (in ppm) to match theoretical and experimental ions (before calibrating)
    :param final_tol: int, error tolerance (in ppm) to calibrate based on median shift
    :param cal_bool: bool, is mass calibration needed?
    :return: analysis_name, all_sites, init_tol, final_tol, cal_bool
    analysis_name = str, pass name
    all_sites = dict, a dictionary with site name (str) as keys, a FragmentSite obj as values
    The objects stored the theoretical ions produced base on the parameters passed on
    init_tol, final_tol, cal_bool = passed a long to match under each pass
    """

    #To time the script
    all_start = time.time()
    print(f"The disulfide list = {ss_ls}")

    #Inititialize a counter to determine the total of theoretical ions created
    counter = 0
    print(f"Protein Length: {len(sequence)   }")

    #Create dictionaries to store sites
    all_sites = {}
    c_sites = {}

    #Produce fragments based on terminal
    for i in range(1, len(sequence)):
        n_frag = sequence[:i]
        if len(n_frag) > 2:
            n_site = fragment_processing(n_frag, sequence, 'N',ss_ls, intrabroken_ss, natredcys)
            all_sites[f"N_{len(n_frag)}"] = n_site

        c_frag = sequence[i:]
        if len(c_frag) > 2:
            c_site = fragment_processing(c_frag, sequence, 'C',ss_ls, intrabroken_ss, natredcys)
            c_sites[f"C_{len(c_frag)}"] = c_site

    #Append the c_sites here so that the dictionary has the N sites first, then the C sites
    all_sites.update(c_sites)

    #Go thru each fragment site in the all_sites dictionary and based on the parameters passed for the current pass give this info to the mass_calc function
    for site in all_sites:
        # print(site)
        #Considere modification and disulfide bonds
        if modbool and cystine:

            # print('cystein/MOO_bool!')

            #Parameters from disulfide analysis
            unboundcys, cysloc, disulfide_counter = all_sites[site].disulfide_analysis()
            # print(f"cysloc = {cysloc}")
            # print(f"disulfide_counter = {disulfide_counter}")

            #If there are cystines to be modified
            if unboundcys == 0:
                mods = []
                if neutral_bool:
                    neutralsLs = ['NH3', 'H2O']
                    for loss in neutralsLs:
                        local_counter = 0
                        counter += mass_calc(local_counter, all_sites[site], ion_types=types, maxcharge=maxcharge,
                                  neutrals=loss, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls=mods,
                                  cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter, mod_ls=noncysmods,
                                             reverse_flag=None, modificatorrepo=libraryofmods)

                else:
                    local_counter = 0
                    counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                              neutrals=neutral_bool, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls=mods,
                              cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter,
                                         mod_ls=noncysmods, reverse_flag=None, modificatorrepo=libraryofmods)

            else:
                for mods in combo_dict[unboundcys]:
                    # print(f"Cys mods {mods}")

                    #Is the removal of water or ammonia needed
                    if neutral_bool:
                        neutralsLs = ['NH3', 'H2O']
                        for loss in neutralsLs:
                            local_counter = 0
                            counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                                      neutrals=loss, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls=mods,
                                      cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter, mod_ls=noncysmods,
                                                  reverse_flag=None, modificatorrepo=libraryofmods)

                    else:
                        local_counter = 0
                        counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                                  neutrals=neutral_bool, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls=mods,
                                  cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter,
                                             mod_ls=noncysmods, reverse_flag=None, modificatorrepo=libraryofmods)

        #Only considere modifications
        elif modbool:
            # print(f'\nMoodBOOL = {noncysmods}')
            if neutral_bool:
                neutralsLs = ['NH3', 'H2O']
                for loss in neutralsLs:
                    local_counter = 0
                    counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                              neutrals=loss, ss_bonds=None, cysmodmass_dict=None, cys_modls = None, cys_num=None,
                              cysloc=None, sscount=None, mod_ls=noncysmods, reverse_flag=None, modificatorrepo=libraryofmods)



            else:
                local_counter = 0
                counter += mass_calc(local_counter,all_sites[site], ion_types = types, maxcharge= maxcharge, neutrals=neutral_bool,ss_bonds=None, cysmodmass_dict=None, cys_modls = None, cys_num=None,
                      cysloc=None, sscount=None, mod_ls=noncysmods,reverse_flag=None, modificatorrepo=libraryofmods)

        #Only do disulfide analysis
        elif cystine:

            # print('cystein_bool!')
            unboundcys, cysloc, disulfide_counter = all_sites[site].disulfide_analysis()
            # print(cysloc)
            # print(f"cysloc = {cysloc}")
            if unboundcys == 0:
                mods = []
                if neutral_bool:
                    neutralsLs = ['NH3', 'H2O']
                    for loss in neutralsLs:
                        local_counter = 0
                        counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                                  neutrals=loss, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls = mods,  cys_num=unboundcys,
                                  cysloc=cysloc, sscount=disulfide_counter, mod_ls=None, reverse_flag=None, modificatorrepo=libraryofmods)

                else:
                    local_counter = 0
                    counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                              neutrals=neutral_bool, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls = mods, cys_num=unboundcys,
                              cysloc=cysloc, sscount=disulfide_counter, mod_ls=None, reverse_flag=None, modificatorrepo=libraryofmods)

            else:
                for mods in combo_dict[unboundcys]:
                    # print(f"Cys mods {mods}")

                    if neutral_bool:
                        neutralsLs = ['NH3', 'H2O']
                        for loss in neutralsLs:
                            local_counter = 0
                            counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                                      neutrals=loss, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls = mods,  cys_num=unboundcys,
                                      cysloc=cysloc, sscount=disulfide_counter, mod_ls=None, reverse_flag=None, modificatorrepo=libraryofmods)

                    else:
                        local_counter = 0
                        counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                                  neutrals=neutral_bool, ss_bonds=True, cysmodmass_dict=combo_dict, cys_modls = mods, cys_num=unboundcys,
                                  cysloc=cysloc, sscount=disulfide_counter, mod_ls=None, reverse_flag=None, modificatorrepo=libraryofmods)

        #No modification or disulfide bond analysis needed
        else:
                # print("Nothing!")
                # print(f'neutral_bool = {neutral_bool}')

                if neutral_bool:
                    neutralsLs = ['NH3', 'H2O']
                    for loss in neutralsLs:
                        local_counter = 0
                        counter += mass_calc(local_counter,all_sites[site], ion_types=types, maxcharge=maxcharge,
                                  neutrals=loss, ss_bonds=None, cysmodmass_dict=None, cys_modls = None, cys_num=None,
                                  cysloc=None, sscount=None, mod_ls=None, reverse_flag=None, modificatorrepo=libraryofmods)



                else:
                    local_counter = 0
                    counter += mass_calc(local_counter,all_sites[site], ion_types = types, maxcharge= maxcharge, neutrals=neutral_bool,ss_bonds=None, cysmodmass_dict=None, cys_modls = None, cys_num=None,
                          cysloc=None, sscount=None, mod_ls=None,reverse_flag=None, modificatorrepo=libraryofmods)


    all_end = time.time() - all_start
    print('Total prediction time: {}'.format(round(all_end,4)))
    print(f'Total theoretical ions {counter}\n')

    # for site in all_sites:
    #     print(site)
    #     for ion in all_sites[site].theo_ions:
    #         print(all_sites[site].theo_ions[ion])


    return analysis_name, all_sites, init_tol, final_tol, cal_bool




def fragment_processing(frag, protein_sequence, terminal, ss_ls,intrabroken_ss, natredcys_ls):
    """
    Function to create fragment site objects
    :param frag: str, sequence
    :param protein_sequence: str, full protein sequence
    :param terminal: str, 'N' or 'C'
    :param ss_ls: ls, a list of the disulfide bonds (set) in the protein
    :param intrabroken_ss: int, how many disulfide bonds (kept intact inside the fragment site) are meant to be broken
    :param natredcys: ls, a list of the location (int) of cysteines  that do not participate in disulfide bonds
    :return: fragment site object
    """

    residue_dict = {"A":[],"R":[],"N":[],"D":[],"C":[],"E":[],"Q":[],"G":[],"H":[],"I":[],"L":[],"K":[],"M":[],"F":[],"P":[],"S":[],"T":[],"W":[],"Y":[],"V":[]}


    # Use regex to locate initial residue location in the whole protein sequence
    interindex = re.compile(frag, re.I)
    mo = interindex.search(protein_sequence)

    # Add one to "translate from python indexing to residue number
    seqstart = mo.start() + 1

    seqend = mo.end()

    # Analysis per Residue
    for residue in residue_dict:
        amino_acid = re.compile(fr'{residue}', re.I)
        # residue locations
        resloc = set()
        iterr = amino_acid.finditer(frag)
        for res in iterr:
            # print(cys.group(), cys.start())
            resloc.add(res.start() + seqstart)
        residue_dict[residue] = resloc
        # print(cysloc)

    # print(f"This is a {terminal}-term site = {frag}_{iseqstart}-{iseqend}. Amino Acid composition = {residue_dict}")

    return FragmentSite(frag, terminal, protein_sequence, seqstart, seqend, residue_dict, ss_ls, intrabroken_ss, natredcys_ls)


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

def write_hit_header():
    """ Returns a header string with column information for the print_hit_info method """
    header_line = 'Pass num,cal mz_mono (exp),mz_mono (thy),cal error(ppm),charge,ion type,mods,losses,neutral mass (thy), free_disulfides, cys_location, disulfide mods,' \
                  'DT mono (bins),Ht (mono),Area (mono),mz_top,DT (top),Ht (top),Area (top),#Peaks in cluster,' \
                  'top pk index,charge,mz_avg,Ht (cluster),Area (cluster),Score (averagine), noise lvl,' \
                  'mz_mono (exp) uncal, error (ppm) uncal' + '\n'
    return header_line


def print_hits(protein_seq, fragsites_input, outputpath):
    """
    Writes match information to text file (comma-delim) using HitObj's print_hit_info() method
    :param protein_seq: string, sequence of the protein in question
    :param all_sites: list of FragmtSite objects containing hit information (previously matched)
    :type all_sites: list[FragmentSite]
    :param outputpath: path to save location for the resulting file
    :return: none (void)
    """

    with open(outputpath, 'w') as outfile:
        # write protein seq and header
        # outfile.write('Protein seq:,' + protein_seq + '\n')
        header_line = write_hit_header()
        outfile.write(header_line)

        for site in fragsites_input:
            # print(f"For each site in the results_ls {site.hits}")
            if len(site.hits) > 0:
                # there is at least one hit here, so write this site to file and print hits
                siteline = site.term + ',' + str(len(site.seq)) + ',' + \
                           site.seq + '\n'
                outfile.write(siteline)
                for hit in site.hits:
                    hitline = hit.print_hit_info()
                    # hitline = str(hit)
                    outfile.write(hitline + '\n')


def print_unmatched(ls_expions, outputpath):
    """
    Writes match information to text file (comma-delim) using HitObj's print_hit_info() method
    :param all_sites: list of FragmtSite objects containing hit information (previously matched)
    :param outputpath: path to save location for the resulting file
    :return: none (void)
    """

    with open(outputpath.strip(".csv") + "_Unmatched" + '.csv', 'w') as outfile:
        # write protein seq and header
        outfile.write('mz,z,int\n')

        for ion in ls_expions:
            out_str = f"{ion.mz_mono},{ion.charge},{ion.pkht_cluster}\n"

            outfile.write(out_str)

    outfile.close()


def save_fragments(fragment_list, output_filename):
    """
    Use Pickle to save the list of matched FragmentSite objects for later retrieval in analysis methods.
    :param fragment_list: list of FragmentSite containers with Hits matched
    :param output_filename: full system path to file to save (convention is '.hits' file)
    :return: void
    """
    with open(output_filename, 'wb') as picklefile:
        pickle.dump(fragment_list, picklefile)

def fragments_to_picklefile(file_title, full_proteinseq, fragments_dict):
    """
    Save the fragment function as a binary file
    :return: void
    """
    outname = file_title + ".ions"
    out_tupsaved = [full_proteinseq, fragments_dict]
    print(f"Saving database file as {outname}")
    with open(outname, 'wb') as picklefile:
        pickle.dump(out_tupsaved, picklefile)

def unmatched_expions_outfile_writer(exp_ion_list, output_filename, csv = False):
    """
    Writes unmatched ions files. These ions will be used in following passes for search.
    :param exp_ion_list: ls, experimental ions
    :param output_filename: str, Name of the file to be written
    :param csv: file extension, if false an unmatched file will be written
    :return: void
    """
    # Test to get the unid exp-ions out and somehow, even if it manually to determine what are the neutral/radicall losses

    if csv:
        charge_dict = {}
        for ion in exp_ion_list:
            if ion.charge in charge_dict:
                charge_dict[ion.charge].append(ion)
            else:
                charge_dict[ion.charge] = [ion]

        # print(charge_dict)

        outstr = ""
        output = open(output_filename.strip(".csv") + "_Unmatched" + '.csv', 'w')
        for ion in exp_ion_list:
            outstr += f"{ion.mz_mono},{ion.charge},{ion.pkht_cluster}\n"

        output.write(outstr)
        output.close()
    else:
        print_unmatched(exp_ion_list, output_filename)
        outnamefile = output_filename.strip(".csv")
        save_fragments(exp_ion_list, outnamefile + ".unmatched")

def main_batch_multipass(main_outdir=None, modificationsrepo=None):
    """
    Main method to run several analysis with several passes for several experimental ion file
    The ions that are matched will not be considered for matching in consecutive passes
    :return: void
    """

    os.chdir(main_outdir)

    #set up Tkinter
    root = tkinter.Tk()
    root.withdraw()


    # load peaklist(s) to search
    # csv files are to expand the type of data inputed into the IMAnnotator (like in dire cases manual peak picked data) - 042120
    exp_files = filedialog.askopenfilenames(title='Choose Peaklist Files',
                                            filetypes=[('IMTBX/Grppr files', '.isotopes'), ('mMass files', '.txt'),
                                                       ('csv files', '.csv'), ('unmatched files', '.unmatched')])




    load_ions_bool = messagebox.askyesno('Load Ions File?',
                                         'Do you want to load an existing .ions file (or calculate theoretical ions fresh from a template file)? Choose YES to load .ions file or NO to open a parameter template')


    # A dict of analysis. Each one has a list of passes
    analysis_dict = {}
    protein_seq = ''
    if load_ions_bool:
        ions_file = filedialog.askopenfilename(title='Choose Ions File', filetypes=[('Ions File', '.ions')])
        with open(ions_file, 'rb') as picklefile:
            saved_database = pickle.load(picklefile)
            protein_seq = saved_database[0]
            analysis_dict = saved_database[1]

    else:

        # Load theoretical database parameters
        paramfile = filedialog.askopenfilename(title='Load Parameter File', filetypes=[('CSV', '.csv')])
        paramfilesplits = paramfile.split('/')
        paramfilesname = paramfilesplits[-1].strip('.csv')
        print(paramfilesname)
        params = Parameter_Parser_terminal.parse_param_template_batch_multipass(paramfile)
        print(f"Params = {params}")
        protein_seq = params[0]
        paramobj_dict = params[1]


        for analysis in paramobj_dict:
            print(f"Current Analysis = {analysis} of {len(paramobj_dict)}")

            #A dictionary of passes per analysis
            passes_dict = {}

            for paramobj in paramobj_dict[analysis]:
                print(f"PassName = {paramobj.analysisName}")

                print("~~~~~~Fragmentation starts~~~~~~")
                pool = multiprocessing.Pool(processes=10)
                # results = []
                argslist = [paramobj.analysisName,paramobj.seq,paramobj.iontypes,paramobj.maxcharge,paramobj.neutraloss_bool,paramobj.disulfide_bool,paramobj.combodict_calc(),
                            paramobj.disulfide_ls,paramobj.ss_allowbroken,paramobj.naturally_redcys,paramobj.mod_bool,paramobj.noncysmods,False,paramobj.init_tol,paramobj.final_tol,paramobj.cal_bool,
                            modificationsrepo]

                pool_result = pool.apply_async(fragments, args=argslist)
                # results.append(pool_result)

                pool.close()  # tell the pool we don't need it to process any more data

                passes_dict[paramobj.analysisName] = pool_result.get()


            #Add theoretical_database to the passes_dict
            analysis_dict[analysis] = passes_dict
        #
            # Create pickle file to save theoretical database


        # print(f'analysis_dict -ion_boolfalse = {analysis_dict}')
        fragments_to_picklefile(f"{paramfilesname}_theoretical_database", protein_seq, analysis_dict)

        # fragments_to_FRAGSfile(f"{analysis}", str(analysis_dict))
    #
    output_path = main_outdir
    # Search peaklists with parameters
    for file_index, exp_file in enumerate(exp_files):
        print('Searching file {} of {}...'.format(file_index + 1, len(exp_files)))

        org_expions, short_filename = Parameter_Parser_terminal.unified_exp_parser(exp_file)

        #sort exp ions so that replicates are truly replicates and go thru the search space equally
        org_expions.sort()

        print(f"File = {short_filename}")
        print(f'Expions before analysis_num loop = {len(org_expions)}')

        # Matching process begins
        for analysis_num in analysis_dict:
            output_filename = output_path + "/" + short_filename + f"-{analysis_num}_hits.csv"
            print(f'Matching currently Analysis # {analysis_num}')
            # Make a copy of the experimental ions list
            analysis_expions = org_expions.copy()
            print(f'Initial Exp ions in {analysis_num} = {len(analysis_expions)}')
            results_ls = []
            # print(f"Analysis Num = {analysis_num}")
            # print(f"Analysis Num type = {type(analysis_num)}")
            for pass_tuple in analysis_dict[analysis_num]:

                print(f"pass_tuple = {pass_tuple}")
                # print(f"analysis_dict[analysis_num][pass_tuple] = {analysis_dict[analysis_num][pass_tuple]}")


                # The matching happes here!
                matched_ions, dict_sites, median_error, average_error = matchmaker_terminal_multipass(
                    analysis_dict[analysis_num][pass_tuple],
                    analysis_expions, None)
                print(f"Matched Ions = {matched_ions}")
                # print(f"Site_dict = {dict_sites}")
                print(f"Median = {median_error}")

                # Removing matched ions from futures passes
                for ion in matched_ions:
                    if ion in analysis_expions:
                        analysis_expions.remove(ion)
                print(f'Expions after {pass_tuple}= {len(analysis_expions)}\n')

                unmatched_expions_outfile_writer(analysis_expions, output_filename)

                # Transform dictionary of FragmentSite Objs into a list of FragmentSite Objs
                if results_ls:
                    for site in results_ls:
                        for fragobj in dict_sites:
                            if site == dict_sites[fragobj]:
                                site.hits.extend(dict_sites[fragobj].hits)

                else:

                    for site in dict_sites:
                        results_ls.append(dict_sites[site])
                    # print(f'Results ls creation = {results_ls}')

                # output_filename = output_path + "/" + short_filename + f"_{pass_tuple}.csv"
                # print_hits(protein_seq, dict_sites, output_filename, input_type = 'dict')

            # Save results

            print_hits(protein_seq, results_ls, output_filename)
            save_fragments(results_ls, output_filename.rstrip('_hits.csv') + '.hits')


if __name__ == '__main__':
    # print( f"z-dot = {mass.Composition(formula='H-2O-1' + 'ON-1H-1')}  'z'= {mass.Composition(formula='H-2O-1' + 'ON-1')}")
    # print(f"c-zdot = {mass.Composition(formula='H-2O-1' + 'NH3'+'H-2O-1' + 'ON-1H-1')} ")
    # print(f"c-z = {mass.Composition(formula='H-2O-1' + 'NH3'+'H-2O-1' + 'ON-1')} ")
    # # print(f"c-y = {mass.Composition(formula='H-2O-1' + 'NH3'+ '')}")
    # print(f"a-z = {mass.Composition(formula='H-2O-1' + 'C-1O-1'+'H-2O-1' + 'ON-1')}")
    # print(f"a-zdot = {mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1' + 'ON-1H-1')}")
    # print(f"'a-y':{mass.Composition(formula='H-2O-1' + 'C-1O-1' + '')}")
    # print(f"'b-y':{mass.Composition(formula='H-2O-1' + '')}")

    main_batch_multipass()


