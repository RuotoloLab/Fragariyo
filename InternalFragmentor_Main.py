"""
Author: Carolina Rojas Ramirez
Date: 07/11/2018
In silico fragmentation of proteins using mass.fast_mass2 from pyteomics
"""

from pyteomics import mass
from pyteomics import parser
from tkinter import filedialog
from tkinter import simpledialog
import tkinter as tk
import os
import pickle
import Parameter_Parser
import re
from tkinter import messagebox
from Modifications import mods_repo
from PeakMatch import matchmaker2
from PeakMatch import matchmaker2_multipass
from Parameter_Parser import expion_parser
from Parameter_Parser import isotope_xtractor



# Setting up tkinter
root = tk.Tk()
root.withdraw()


# update Pyteomics masses
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
    'c-z': mass.Composition(formula='H-2O-1' + 'NH3'+'H-2O-1' + 'ON-1'),
    'c-zdot': mass.Composition(formula='H-2O-1' + 'NH3'+'H-2O-1' + 'ON-1H-1'),
    'cdot-z': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1'+'H-2O-1' + 'ON-1'),
    'c-y': mass.Composition(formula='H-2O-1' + 'NH3'+ ''),
    'cdot-y': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1'+''),
    'a-z': mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1' + 'ON-1'),
    'a-zdot': mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1' + 'ON-1H-1')
    })




class interfrag:
    """
    Container object for internal fragments and their properties
    :param mz_mono: int, monoisotope m/z
    :param charge: int, ion charge
    :param ion type: str
    :param sequence: str, internal fragment
    :param cys_num: int, cysteines available for modification
    :param mono_neutral: int, ion neutral mass
    :param mods: ls, mods if any
    :param cysloc: set, locations of cysteines
    :param ss_count: int, amount of disulfide bonds in the internal fragment
    :param ifragstart: Index where the fragment starts in the protein sequence
    :param ifragend: Index where the fragment ends in the protein sequence
    :param reverse: bool, is the sequence a reverse version.
    """

    def __init__(self, mz_mono, charge, ion_type, sequence, cys_num, mono_neutral, mods, cysloc, ss_count, ifragstart, ifragend, reverse):
        self.mz_mono = round(mz_mono, 8)
        self.charge = charge
        self.ion_type = ion_type
        self.sequence = sequence
        self.cys_num = cys_num
        self.mono_neutral = round(mono_neutral,8)
        self.mods = mods
        self.cysloc = cysloc
        self.ss_count = ss_count

        self.ifragstart = ifragstart
        self.ifragend = ifragend
        self.reverse = reverse


    #Attribute to be able to compare to other objects
    def __eq__(self, other):
        return self.mono_neutral == other.mono_neutral

    #Attribute to represent objects in dictionaries
    def __hash__(self):
        return hash(self.mono_neutral)

    # Attribute to represent objects in stdout
    def __str__(self):
        return f"{self.mono_neutral}\t{self.sequence}\t{self.charge}\t{self.mz_mono}\t{self.mods}\t{self.ion_type}\t{self.cysloc}\t{self.ss_count}\t{self.cys_num}\t{self.ifragstart}\t{self.ifragend}\t{self.reverse}"

    # Attribute to represent objects in other objects
    __repr__ = __str__

def mass_calc(seq, types, mincharge, maxcharge, dictionary, ss_bonds=None, cysmod_dict = None, cys_num=None, cysloc = None, sscount = None, iseqstart=None, iseqend=None, reverse_flag=None, modobj =None):
    """
    Function to calculate masses of specific internal fragments with specific disulfides configurations
    :param seq: str, interfragment sequence
    :param types: ls, ions types
    :param mincharge: Minimun charge
    :param maxcharge: Maximun charge
    :param ss_bonds: bool, disulfide bonds present
    :param cysmod_dict: Dictionary with number of disulfides available for modification (keys), possible modifications (values)
    :param cys_num: Cysteines available for modification
    param cysloc: set, location of cysteines
    :param ss_count: int, number of disulfide bonds
    :param iseqstart:  Index where the internal fragment starts at in the protein sequence
    :param iseqend:  Index where the internal fragment ends at in the protein sequence
    :param reverse_flag:  bool, is this a reverse sequence?
    :param mods = ls of noncys mods
    """

    #For each ion type
    for ion_type in types:
        #For each charge
        for charge in range(mincharge, maxcharge + 1):

            #calculate m/z and neutral
            mz_mono_org = mass.fast_mass2(seq, ion_type=ion_type, charge=charge)
            neutral_mono_org = mass.fast_mass2(seq, ion_type=ion_type, charge=0)

            #If disulfides are considered
            if ss_bonds:

                #If the sequence has not free cysteines for modification
                if cys_num == 0:
                    if modobj:
                        # print(modobj)
                        for modification in modobj:
                            mods_repo[modification].current_num = 0
                            if modification == 'oxyhemeC':
                                # print("HEMEC!!!")
                                aa_loc = set()
                                for amino in mods_repo[modification].target_aas:

                                    mod_regex = re.compile(amino, re.I)
                                    locations = mod_regex.finditer(seq)
                                    # print(locations)

                                    for loc in locations:
                                        # print(cys.group(), cys.start())
                                        aa_loc.add(loc.start() + iseqstart)
                                    print(f"Amino acid {amino} at {aa_loc} for modification by {modification}")
                                    if len(aa_loc) > 1:
                                        print("Add oxyhemeC!")

                                        # Add HemeC and remove the two hydrogen from the protein participating in binding HemeC
                                        neutral_mono = neutral_mono_org + mods_repo[modification].mass + (-1.0078 * 2)
                                        mz_mono = mz_mono_org + (mods_repo[modification].mass / charge) + (
                                                    (-1.0078 / charge) * 2)

                                        #HemeC has been added
                                        sscount -= 1

                                        if sscount < 0:
                                            sscount = 0

                                        # Removing the hydrogens due to other disulfide bonds
                                        neutral_mono = neutral_mono + (sscount * (-1.0078 * 2))
                                        mz_mono = mz_mono + (sscount * ((-1.0078 / charge) * 2))


                                        # Store internal fragment information
                                        interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono,
                                                             modification, cysloc,
                                                             sscount, iseqstart, iseqend, reverse_flag)
                                        print(interobj)
                                        dictionary[interobj.mono_neutral] = interobj

                                    else:
                                        print("Don't add oxyHemeC!")
                                        modification = ""
                                        print(f"sscount = {sscount}")
                                        # except removing the hydrogens due to the disulfide bonds
                                        neutral_mono = neutral_mono_org + (sscount * (-1.0078 * 2))
                                        mz_mono = mz_mono_org + (sscount * ((-1.0078 / charge) * 2))

                                        # Store internal fragment information
                                        interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono,
                                                             modification, cysloc,
                                                             sscount, iseqstart, iseqend, reverse_flag)
                                        print(interobj)
                                        dictionary[interobj.mono_neutral] = interobj


                            elif mods_repo[modification].target_aas:
                                for amino in mods_repo[modification].target_aas:

                                    mod_regex = re.compile(amino, re.I)
                                    mo = mod_regex.search(seq)
                                    # print(mo)
                                    print(f"In the  fragment {seq}, {amino} is at positions {mo}")
                                    if mo:
                                        # print(mo.start())
                                        # It considers uniprot
                                        proteinseqstart = iseqstart
                                        # print(proteinseqstart)

                                        print(f"Amino acid {amino} at {mo.start() + proteinseqstart} for modification by {modification}")

                                        # except removing the hydrogens due to the disulfide bonds
                                        neutral_mono = neutral_mono_org + sscount * (-1.0078 * 2)
                                        mz_mono = mz_mono_org + (sscount * ((-1.0078 / charge) * 2))

                                        # Add HemeC and remove the two hydrogen from the protein participating in binding HemeC
                                        neutral_mono = neutral_mono + modobj.mass + (-1.0078 * 2)
                                        mz_mono = mz_mono + (modobj.mass / charge) + ((-1.0078 / charge) * 2)

                                        # Store internal fragment information
                                        interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono,
                                                             modobj, cysloc,
                                                             sscount, iseqstart, iseqend, reverse_flag)
                                        dictionary[interobj.mono_neutral] = interobj

                                    else:
                                        print(f"No amino acid for modification by {modification}")
                                        # except removing the hydrogens due to the disulfide bonds
                                        neutral_mono = neutral_mono_org + sscount * (-1.0078 * 2)
                                        mz_mono = mz_mono_org + (sscount * ((-1.0078 / charge) * 2))

                                        # Store internal fragment information
                                        interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono,
                                                             modobj, cysloc,
                                                             sscount, iseqstart, iseqend, reverse_flag)
                                        dictionary[interobj.mono_neutral] = interobj

                                        print(interobj)



                    else:
                        #Do not modified cysteines
                        mods = ""

                        #except removing the hydrogens due to the disulfide bonds
                        neutral_mono = neutral_mono_org + sscount*(-1.0078*2)
                        mz_mono = mz_mono_org + (sscount * ((-1.0078/charge) * 2))

                        #Store internal fragment information
                        interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono, mods, cysloc, sscount, iseqstart, iseqend, reverse_flag)
                        dictionary[interobj.mono_neutral] = interobj
                else:

                    print("Allowing ssbond broken, handle special HemeC!")

                    if 'oxyhemeC' in modobj:
                        print(modobj)

                        #Cause HemeC makes a "special ss"
                        sscount=-1

                        if sscount < 0:
                            sscount = 0

                        for mods in cysmod_dict[cys_num]:
                            # Modify neutral mass with the modifications possible
                            neutral_mono = neutral_mono_org + cysmod_dict[cys_num][mods] + (sscount * (-1.0078 * 2))
                            mz_mono = mz_mono_org + (cysmod_dict[cys_num][mods] / charge) + (
                                    sscount * ((-1.0078 / charge) * 2))

                            # Store internal fragment information
                            interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono, mods, cysloc,
                                                 sscount, iseqstart,
                                                 iseqend, reverse_flag)
                            print(interobj)
                            dictionary[interobj.mono_neutral] = interobj


                    else:

                        for mods in cysmod_dict[cys_num]:
                            # Modify neutral mass with the modifications possible
                            neutral_mono = neutral_mono_org + cysmod_dict[cys_num][mods] + (sscount * (-1.0078 * 2))
                            mz_mono = mz_mono_org + (cysmod_dict[cys_num][mods] / charge) + (
                                    sscount * ((-1.0078 / charge) * 2))

                            # Store internal fragment information
                            interobj = interfrag(mz_mono, charge, ion_type, seq, cys_num, neutral_mono, mods, cysloc,
                                                 sscount, iseqstart,
                                                 iseqend, reverse_flag)
                            print(interobj)
                            dictionary[interobj.mono_neutral] = interobj



            else:
                print("No cysmods, but mods")
                if modobj:
                    print(modobj)
                    for modification in modobj:
                        mods_repo[modification].current_num = 0
                        if mods_repo[modification].target_aas:
                            for amino in mods_repo[modification].target_aas:
                                print(f"In the  fragment {seq}, {amino} is at positions ???")
                                mod_regex = re.compile(amino, re.I)
                                mo = mod_regex.search(seq)
                                # print(mo)
                                # print(mo.start())
                                #It considers uniprot
                                proteinseqstart = iseqstart
                                # print(proteinseqstart)
                else:
                    print("Regular seq!")
                    mods = ""
                    sscount = 0
                    # Store internal fragment information
                    interobj = interfrag(mz_mono_org, charge, ion_type, seq, cys_num, neutral_mono_org, mods, cysloc,
                                         sscount, iseqstart,
                                         iseqend, reverse_flag)
                    print(interobj)
                    dictionary[interobj.mono_neutral] = interobj




def ifragments(analysis_name, sequence, types=('b', 'y'), mincharge = 0, maxcharge=1, maxstart = 4, maxlength=10, modbool = None, max_mods=1, combo_dict = None, cystine = None, uniprot_offset=None, allow_ssbroken = None, reverse_seq=None, foo_ls=None, reduced_cys=None):
    """
    :param analysis_name: Str, the name of the current analysis or pass
    :param sequence: a string of the sequence desired to fragment.
    :param types: the ion types like b or y. To pass more than one type of ion put them in parenthesis
    :param mincharge:
    :param maxcharge: the maximun charge for the fragments.
    :param maxstart: int, the smallest internal fragment size (number of residues)
    :param maxlength: int, the largest internal fragment size (number of residues)
    :param modbool: Modifications
    :param maxmods:
    :param combo_dict: A dictionary with all the possible combinations of disulfide bon mods
    :param cystine: Bool, considered disulfie bonds
    :param uniprot_offset: The amount of residues the experimental sequence is offset from the its uniprot sequence and be able to properly locate disulfide bonds
    Important in order to match disulfide bond location
    :param allow_ssbroken: How many disulfide bonds within an internal fragment
    :param reverse_seq: If it is true it means that the protein sequence are reverse
    :param foo_ls: Disulfide ls
    :return: A dictionary with neutral masses as keys and interobj as values
    """
    print('---------Fragmentation starts--------')


    global_dict = {}
    counter = 0
    #creating internal sequences
    interfrag_list = []
    isequence = sequence[0:-1]

    for i in range(1, len(isequence)):

        #Proline effect
        if isequence[i] == "P":
            continue
        else:

            for n in range(maxstart, maxlength+1):
                interseq = isequence[i:i+n]

                #Eliminating redundancies
                if interseq not in interfrag_list:
                    interfrag_list.append(interseq)
                else:
                    continue

    # print(interfrag_list)

    #For sequences with desired length
    for iseq in interfrag_list:
        if len(iseq) > maxstart:
            print(f"\n {iseq}")
            #Use regex to locate initial residue location in the whole protein sequence
            interindex = re.compile(iseq, re.I)
            mo = interindex.search(sequence)
            # print(mo.start()+1)

            #Add one to "translate from python indexing to residue number
            #Add uniport offset to match uniport residue numbers and properly id disulfuide bonds
            iseqstart = mo.start()+ 1 +uniprot_offset
            iseqend = mo.end() + 1 + uniprot_offset
            # print(iseqindex)

            # print(f" Before modbool loop Modbool = {modbool}, lenght = {len(modbool)}")

            if cystine:
                print("Considering Cysteines!")

                unboundcys = 0

                # disulfide analysis
                Cysteines = re.compile(r'C', re.I)
                #cysteine locations
                cysloc = set()
                iterc = Cysteines.finditer(iseq)
                for cys in iterc:
                    # print(cys.group(), cys.start())
                    cysloc.add(cys.start() + iseqstart)


                # The number of cysteines in the internal fragment to be considered for disulfides and for modification
                cysloclen = len(cysloc)

                print(f"cysloc = {cysloc} and its lenght {len(cysloc)}")
                # print(f"foo_ls = {foo_ls}")

                disulfide_counter = 0
                # foos_ls is the list of disulfide bond pairs from uniprot
                for pair in foo_ls:
                    if pair.issubset(cysloc):
                        # print(disulfide)
                        disulfide_counter += 1

                print(f"disulfide count = {disulfide_counter}")

                # naturally_reducedcys_ls = [58, 476]
                # Carefull taht is a list of str
                print(f"reduced_cys = {reduced_cys}")

                naturally_reducedcys_ls = reduced_cys

                for cys in cysloc:
                    if int(cys) in naturally_reducedcys_ls:
                        # print(cys)
                        cysloclen -= 1

                # Allow disulfides inside the internal fragment to be broken
                if allow_ssbroken and cysloclen>1:
                    # Unbound cysteines are cysteines that are free, don't have their didulfide bond partner in the current fragment, or
                    # cysteines that were reduced (allow_ssbroken)
                    unboundcys = cysloclen - (disulfide_counter * 2) + (allow_ssbroken * 2)
                else:
                    # Considering disulfides that reached outside the internal fragment to be broken
                    # Cysteines available for modification
                    unboundcys = cysloclen - (disulfide_counter * 2)

                print(f"unboundcys = {unboundcys}")

                # Id the amount of unbound cys is negative it means there is no cys for modification
                if unboundcys < 0:
                    unboundcys = 0

                if len(modbool) >=1:
                    print("Considering Cysteines and modifications!")
                    mass_calc(iseq, types, mincharge, maxcharge, global_dict, ss_bonds=True, cysmod_dict=combo_dict,
                              cys_num=unboundcys, cysloc=cysloc, sscount=disulfide_counter, iseqstart=iseqstart,
                              iseqend=iseqend, reverse_flag=reverse_seq, modobj=modbool)
                    counter += 1
                else:

                    print("Considering Cysteines and no modifications!")



                    mass_calc(iseq, types, mincharge, maxcharge, global_dict, ss_bonds=True, cysmod_dict=combo_dict, cys_num=unboundcys, cysloc = cysloc, sscount=disulfide_counter, iseqstart=iseqstart, iseqend=iseqend, reverse_flag=reverse_seq)
                    counter += 1

            else:
                print("Not Considering Cysteines!")
                if len(modbool) >=1:
                    print("Considering modifications!")
                    mass_calc(iseq, types, mincharge, maxcharge, global_dict, ss_bonds=None, cysmod_dict=None,
                              cys_num=None, cysloc=None, sscount=None, iseqstart=iseqstart, iseqend=iseqend, reverse_flag=reverse_seq,
                              modobj=modbool)
                    counter += 1
                else:
                    print("Regular primary sequence!")
                    mass_calc(iseq, types, mincharge, maxcharge, global_dict, reverse_flag=reverse_seq)
                    counter += 1



    print(f"{counter} Internal fragments were produced - does not include modification isoforms\n")
    # print(global_dict)
    return (analysis_name, global_dict)




def fragments_to_picklefile(file_title, fragments_dict):
    """
    Save the fragment function as a binary file
    :return:
    """
    outname = file_title + ".theofrags"
    with open(outname, 'wb') as picklefile:
        pickle.dump(fragments_dict, picklefile)

def fragments_to_FRAGSfile(file_title, fragments):
    """
    Save the fragment function output as a .txt file (recommended only for small proteins)
    :param file_title: Str, title of the file
    :param fragments: Str, a string with the fragments
    :return: void
    """

    output = open(str(file_title) + '.txt.frags', 'w')
    output.write(fragments)
    output.close()

def main():
    """
    Main method to run an analysis on one experimental ions file
    :return: void
    """
    expions, proteinName = expion_parser()

    load_ions_bool = messagebox.askyesno('Load Theoretical Database ?',
                                         'Do you want to load an existing .theofrags file (or calculate theoretical ions fresh from a template file)? Choose YES to load .ions file or NO to open a parameter template')

    if load_ions_bool:
        print("Loading Theoretical Database")
        theofrags_file = filedialog.askopenfilename(title='Load Theoretical Database',
                                                    filetypes=[('Theoretical Fragments', '.theofrags')])
        with open(theofrags_file, 'rb') as picklefile:
            frag_dict = pickle.load(picklefile)

        print("Loading Reverse Theoretical Database")

        theofrags_file = filedialog.askopenfilename(title='Load Reverse Sequence Theoretical Database',
                                                    filetypes=[
                                                        ('Reverse Sequence Theoretical Fragments', '.theofrags')])
        with open(theofrags_file, 'rb') as picklefile:
            frag_dict_rev = pickle.load(picklefile)





    else:
        paramfiles = filedialog.askopenfilename(title='Load Parameter File', filetypes=[('CSV', '.csv')])
        paramobj = Parameter_Parser.parse_param_template(paramfiles)

        paramfiles = filedialog.askopenfilename(title='Load Reverse Sequence Parameter File',
                                                filetypes=[('CSV', '.csv')])
        paramobj_rev = Parameter_Parser.parse_param_template(paramfiles)

        frag_dict = ifragments(analysis_name=paramobj.analysisName, sequence=paramobj.seq, types=paramobj.iontypes, mincharge=paramobj.mincharge,
                               maxcharge=paramobj.maxcharge,
                               maxstart=paramobj.min_len, maxlength=paramobj.max_len, modbool=None, max_mods=0,
                               combo_dict=paramobj.combodict_calc(), cystine=paramobj.disulfide_bool,
                               uniprot_offset=paramobj.uniprot_offset, allow_ssbroken=paramobj.ss_allowbroken,
                               reverse_seq=False, foo_ls=paramobj.disulfide_ls)
        # Create pickle file
        fragments_to_picklefile(f"{paramobj.proteinName}_{paramobj.ss_allowbroken}_z{paramobj.maxcharge}", frag_dict)


        # assert paramobj_rev.seq == paramobj.seq[::-1]
        # print(paramobj.seq)
        # print(paramobj_rev.seq)
        print("Calculating Reverse Sequence Theoretical Fragments")

        # producing_fragments with a reverse sequence
        frag_dict_rev = ifragments(analysis_name=paramobj_rev.analysisName,sequence = paramobj_rev.seq, types=paramobj_rev.iontypes, mincharge=paramobj_rev.mincharge,
                                   maxcharge=paramobj_rev.maxcharge, maxstart=paramobj_rev.min_len,
                                   maxlength=paramobj_rev.max_len,
                                   modbool=None, max_mods=0, combo_dict=paramobj_rev.combodict_calc(),
                                   cystine=paramobj_rev.disulfide_bool, uniprot_offset=paramobj_rev.uniprot_offset,
                                   allow_ssbroken=paramobj_rev.ss_allowbroken, reverse_seq=True,
                                   foo_ls=paramobj_rev.disulfide_ls)
        # Create pickle file
        fragments_to_picklefile(f"{paramobj.proteinName}_Reverse_{paramobj.ss_allowbroken}_z{paramobj.maxcharge}",
                                frag_dict_rev)


    # Set reverse seqeunce theoretical fragment amount as the same as the non-reverse theoretical fragment database
    frag_dict_rev_mod = {}
    print(f"Frag_dict = {len(frag_dict)}")
    print(f"Frag_dict_rev = {len(frag_dict_rev)}")
    if len(frag_dict) != len(frag_dict_rev):
        frag_dict_rev_keys = frag_dict_rev.keys()
        counter = 0
        for x in frag_dict_rev_keys:
            if counter < len(frag_dict) + 1:
                frag_dict_rev_mod[x] = frag_dict_rev[x]
                counter += 1
    print(f"Frag_dict_rev_mod = {len(frag_dict_rev_mod)}")

    # Update the theoretical fragmentation dictionary with the reverse sequence theoretical fragmentation dictionary
    # frag_dict.update(frag_dict_rev_mod)
    # frag_dict_new = {**frag_dict, **frag_dict_rev}

    outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(outdir)

    # fragments_to_FRAGSfile(proteinName, str(frag_dict_rev))

    error_tol = simpledialog.askstring("Error_tolerance", "How much error tolerance (ppm) to use when matching?")
    analysisName = simpledialog.askstring("AnalysisName", "How do you want to name this analysis?")
    analysisSSbroken = simpledialog.askstring("SS-bonds broken", "How many SS bonds broken within an internal fragment?")
    analysischarge = simpledialog.askstring("Charge state", "What is the charge state of the experimental ions?")

    matchmaker2(frag_dict, frag_dict_rev_mod, expions, error_tol,analysisName, analysisSSbroken, analysischarge)


def main_batch():
    """
    Main method to run several analysis for one experimental ion file
    :return:
    """
    expions, proteinName = expion_parser()

    paramfile = filedialog.askopenfilename(title='Load Parameter File', filetypes=[('CSV', '.csv')])



    paramobj_dict = Parameter_Parser.parse_param_template_batch(paramfile)



    print(paramobj_dict)

    paramfile_rev = filedialog.askopenfilename(title='Load Reverse Sequence Parameter File',
                                               filetypes=[('CSV', '.csv')])
    paramobj_rev_dict = Parameter_Parser.parse_param_template_batch(paramfile_rev)

    main_outdir = filedialog.askdirectory(title='Choose Output Folder')

    analysis_counter = 1
    for paramobj in paramobj_dict:
        os.chdir(main_outdir)
        print(f"Current Analysis = {analysis_counter}")
        frag_dict = ifragments(analysis_name=paramobj_dict[paramobj].analysisName, sequence=paramobj_dict[paramobj].seq, types=paramobj_dict[paramobj].iontypes,
                               mincharge=paramobj_dict[paramobj].mincharge,
                               maxcharge=paramobj_dict[paramobj].maxcharge,
                               maxstart=paramobj_dict[paramobj].min_len, maxlength=paramobj_dict[paramobj].max_len,
                               modbool=None, max_mods=0,
                               combo_dict=paramobj_dict[paramobj].combodict_calc(),
                               cystine=paramobj_dict[paramobj].disulfide_bool,
                               uniprot_offset=paramobj_dict[paramobj].uniprot_offset,
                               allow_ssbroken=paramobj_dict[paramobj].ss_allowbroken,
                               reverse_seq=False, foo_ls=paramobj_dict[paramobj].disulfide_ls)
        # Create pickle file
        os.mkdir(f"{paramobj_dict[paramobj].analysisName}")
        os.chdir(main_outdir + f"/{paramobj_dict[paramobj].analysisName}")
        fragments_to_picklefile(
            f"{paramobj_dict[paramobj].analysisName}_{paramobj_dict[paramobj].ss_allowbroken}_z{paramobj_dict[paramobj].maxcharge}", frag_dict)

        # assert paramobj_rev.seq == paramobj.seq[::-1]
        # print(paramobj.seq)
        # print(paramobj_rev.seq)
        print("Calculating Reverse Sequence Theoretical Fragments")

        # producing_fragments with a reverse sequence
        paramobj_rev = paramobj_rev_dict[paramobj_dict[paramobj].analysisName]
        frag_dict_rev = ifragments(analysis_name=paramobj_rev.analysisName, sequence = paramobj_rev.seq, types=paramobj_rev.iontypes, mincharge=paramobj_rev.mincharge,
                                   maxcharge=paramobj_rev.maxcharge, maxstart=paramobj_rev.min_len,
                                   maxlength=paramobj_rev.max_len,
                                   modbool=None, max_mods=0, combo_dict=paramobj_rev.combodict_calc(),
                                   cystine=paramobj_rev.disulfide_bool, uniprot_offset=paramobj_rev.uniprot_offset,
                                   allow_ssbroken=paramobj_rev.ss_allowbroken, reverse_seq=True,
                                   foo_ls=paramobj_rev.disulfide_ls)
        # Create pickle file
        fragments_to_picklefile(
            f"{paramobj_rev.analysisName}_Reverse_{paramobj_rev.ss_allowbroken}_z{paramobj_rev.maxcharge}",
            frag_dict_rev)

        print(f"Current Analysis = {paramobj_dict[paramobj].analysisName}")
        # Set reverse seqeunce theoretical fragment amount as the same as the non-reverse theoretical fragment database
        frag_dict_rev_mod = {}
        print(f"Frag_dict = {len(frag_dict[1])}")
        print(f"Frag_dict_rev = {len(frag_dict_rev[1])}")
        if len(frag_dict[1]) != len(frag_dict_rev[1]):
            frag_dict_rev_keys = frag_dict_rev[1].keys()
            counter = 0
            for x in frag_dict_rev_keys:
                if counter < len(frag_dict[1]) + 1:
                    frag_dict_rev_mod[x] = frag_dict_rev[1][x]
                    counter += 1
        print(f"Frag_dict_rev_mod = {len(frag_dict_rev_mod)}")


        # error_tol = simpledialog.askstring("Error_tolerance", "How much error tolerance (ppm) to use when matching?")


        matchmaker2(frag_dict[1], frag_dict_rev_mod, expions, 1, paramobj_dict[paramobj].analysisName,
                    paramobj_dict[paramobj].ss_allowbroken, paramobj_dict[paramobj].maxcharge)
        analysis_counter += 1



def main_batch_multipass():
    """
    Main method to run several analysis with several passes for one experimental ion file
    The ions that are matched will not be considered for matching in consecutive passes
    :return: void
    """

    #Parse experimental ions
    org_expions, proteinName = expion_parser()



    #Create folder to save results for one experimental ion file
    main_outdir = filedialog.askdirectory(title='Choose Output Folder')

    #Use previously created fragments
    load_ions_bool = messagebox.askyesno('Load Ions File?',
                                         'Do you want to load an existing .theofrags file (or calculate theoretical ions fresh from a template file)? Choose YES to load .ions file or NO to open a parameter template')

    # A dict of analysis. Each one has a list of passes
    analysis_dict = {}

    # A dict of analysis. Each one hasdictionary of passes
    analysis_dict_rev = {}

    if load_ions_bool:
        theofrags_file = filedialog.askopenfilename(title='Choose Ions File', filetypes=[('Theoretical Fragments File', '.theofrags')])
        with open(theofrags_file, 'rb') as picklefile:
            analysis_dict = pickle.load(picklefile)

        revtheofrags_file = filedialog.askopenfilename(title='Choose Reverse Ions File',
                                                    filetypes=[('Theoretical Fragments File', '.theofrags')])
        with open(revtheofrags_file, 'rb') as revpicklefile:
            analysis_dict_rev = pickle.load(revpicklefile)




    else:

        # Load theoretical database parameters
        paramfile = filedialog.askopenfilename(title='Load Parameter File', filetypes=[('CSV', '.csv')])
        paramobj_dict = Parameter_Parser.parse_param_template_batch_multipass(paramfile)
        print(paramobj_dict)

        # Load theoretical database parameters - reverse sequences
        paramfile_rev = filedialog.askopenfilename(title='Load Reverse Sequence Parameter File',
                                                   filetypes=[('CSV', '.csv')])
        paramobj_rev_dict = Parameter_Parser.parse_param_template_batch_multipass(paramfile_rev)
        print(paramobj_rev_dict)

        for analysis in paramobj_dict:
            print(f"Current Analysis = {analysis} of {len(paramobj_dict)}")
            #Alwasy reset to this output path so that all passes are saved under one folder
            os.chdir(main_outdir)

            #A dictionary of passes per analysis
            passes_dict = {}
            for paramobj in paramobj_dict[analysis]:
                print(f"PassName = {paramobj.analysisName}")
                print(f"before ifragments = {paramobj.combodict_calc()}")


                # It actually returns a tuple with the analysis name in first position then the fragments dictionary
                frag_dict = ifragments(analysis_name=paramobj.analysisName, sequence=paramobj.seq, types=paramobj.iontypes,
                                       mincharge=paramobj.mincharge,
                                       maxcharge=paramobj.maxcharge,
                                       maxstart=paramobj.min_len, maxlength=paramobj.max_len,
                                       modbool=paramobj.noncysmods, max_mods=0,
                                       combo_dict=paramobj.combodict_calc(),
                                       cystine=paramobj.disulfide_bool,
                                       uniprot_offset=paramobj.uniprot_offset,
                                       allow_ssbroken=paramobj.ss_allowbroken,
                                       reverse_seq=False, foo_ls=paramobj.disulfide_ls, reduced_cys=paramobj.naturally_redcys)
                passes_dict[paramobj.analysisName] = frag_dict

            analysis_dict[analysis] = passes_dict

            # Create pickle file to save theoretical database
            out_folder_str = f"{proteinName}_{analysis}"
            os.mkdir(out_folder_str)
            os.chdir(main_outdir + f"/{out_folder_str}")
            fragments_to_picklefile(f"{out_folder_str}", analysis_dict)
            fragments_to_FRAGSfile(f"{out_folder_str}", str(analysis_dict))

        # print(f'analysis_dict = {analysis_dict}')



        print("Calculating Reverse Sequence Theoretical Fragments")

        for rev_analysis in paramobj_rev_dict:
            print(f"Current Analysis for Reverse Sequences = {rev_analysis} of {len(paramobj_dict)}")
            passes_rev_dict = {}
            for paramobj_rev in paramobj_rev_dict[rev_analysis]:
                print(f"PassName = {paramobj_rev.analysisName}")

                #It actually returns a tuple with the analysis name in first position then the fragments dictionary
                frag_dict_rev = ifragments(analysis_name= paramobj_rev.analysisName, sequence=paramobj_rev.seq, types=paramobj_rev.iontypes, mincharge=paramobj_rev.mincharge,
                                       maxcharge=paramobj_rev.maxcharge, maxstart=paramobj_rev.min_len,
                                       maxlength=paramobj_rev.max_len,
                                       modbool=paramobj_rev.noncysmods, max_mods=0, combo_dict=paramobj_rev.combodict_calc(),
                                       cystine=paramobj_rev.disulfide_bool, uniprot_offset=paramobj_rev.uniprot_offset,
                                       allow_ssbroken=paramobj_rev.ss_allowbroken, reverse_seq=True,
                                       foo_ls=paramobj_rev.disulfide_ls, reduced_cys=paramobj_rev.naturally_redcys)

                passes_rev_dict[paramobj_rev.analysisName] = frag_dict_rev

            #Analsis dictionary with passes dictionaries as values
            analysis_dict_rev[rev_analysis] = passes_rev_dict
            out_folder_str_rev = f"{proteinName}_{rev_analysis}"
            os.chdir(main_outdir + f"/{out_folder_str_rev}")
            # Create pickle file
            fragments_to_picklefile(f"{out_folder_str_rev}_Reverse", analysis_dict_rev)

        # print(f'analysis_dict_rev = {analysis_dict_rev}')
        print(f'Expions before analysis_num loop = {len(org_expions)}')






    #Matching process begins
    for analysis_num in analysis_dict:
        print(f'Matching currently Analysis # {analysis_num}')
        #Make a copy of the experimental ions list
        analysis_expions = org_expions.copy()
        print(f'Initial Exp ions in {analysis_num} = {len(analysis_expions)}')
        for pass_tuple in analysis_dict[analysis_num]:
            print(f"go back to folder = {proteinName}_{analysis_num}")
            # print(f"analysis_dict[analysis_num] = {analysis_dict[analysis_num]}")
            print(f"pass_tuple = {pass_tuple}")



            try:
                os.chdir(main_outdir + f"/{proteinName}_{analysis_num}")
            except FileNotFoundError:
                os.mkdir(main_outdir + f"/{proteinName}_{analysis_num}")
                os.chdir(main_outdir + f"/{proteinName}_{analysis_num}")


            # print(f"analysis_dict[analysis_num][pass_tuple] = {analysis_dict[analysis_num][pass_tuple]}")
            matched_expions = matchmaker2_multipass(analysis_dict[analysis_num][pass_tuple], analysis_dict_rev[analysis_num][pass_tuple] , analysis_expions, 10)
            #Removing matched ions from futures passes
            for ion in matched_expions:
                if ion in analysis_expions:
                    analysis_expions.remove(ion)
            print(f'Expions after {pass_tuple}= {len(analysis_expions)}\n')



if __name__ == '__main__':
    # isotope_xtractor()
    main_batch_multipass()
    # main_batch()
    # main()












