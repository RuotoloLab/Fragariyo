"""
Author: Carolina Rojas Ramirez
Date: May 13, 2020
Python module with all the function needed for Peak Matching using the Fragmentor and the terminalFragmentor
"""
import os
import scipy
import numpy as np
from scipy import signal
from pyteomics import mass
from pyteomics import parser
from matplotlib import pyplot as plt
from brainpy import isotopic_variants
import pythoms.molecule as pythmole
from Modifications import mods_repo

iontypedict = {
    'M': mass.Composition(formula=''),
    'M-H2O': mass.Composition(formula='H-2O-1'),
    'M-NH3': mass.Composition(formula='N-1H-3'),
    'a': mass.Composition(formula='H-2O-1' + 'C-1O-1'),
    'adot': mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H1'),
    'a-H2O': mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
    'a-NH3': mass.Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
    'b': mass.Composition(formula='H-2O-1'),
    'b-H2O': mass.Composition(formula='H-2O-1' + 'H-2O-1'),
    'b-NH3': mass.Composition(formula='H-2O-1' + 'N-1H-3'),
    'c': mass.Composition(formula='H-2O-1' + 'NH3'),
    'c-1': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1'),
    'c-dot': mass.Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+1': mass.Composition(formula='H-2O-1' + 'NH3' + 'H1'),
    'c+2': mass.Composition(formula='H-2O-1' + 'NH3' + 'H2'),
    'c-H2O': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1'),
    'c-NH3': mass.Composition(formula='H-2O-1'),
    'x': mass.Composition(formula='H-2O-1' + 'CO2'),
    'x-H2O': mass.Composition(formula='H-2O-1' + 'CO2' + 'H-2O-1'),
    'x-NH3': mass.Composition(formula='H-2O-1' + 'CO2' + 'N-1H-3'),
    'y': mass.Composition(formula=''),
    'y-H2O': mass.Composition(formula='H-2O-1'),
    'y-NH3': mass.Composition(formula='N-1H-3'),
    'z': mass.Composition(formula='H-2O-1' + 'ON-1H-1'),
    'z-dot': mass.Composition(formula='H-2O-1' + 'ON-1'),
    'z+1': mass.Composition(formula='H-2O-1' + 'ON-1H1'),
    'z+2': mass.Composition(formula='H-2O-1' + 'ON-1H2'),
    'z+3': mass.Composition(formula='H-2O-1' + 'ON-1H3'),
    'z-H2O': mass.Composition(formula='H-2O-1' + 'ON-1H-1' + 'H-2O-1'),
    'z-NH3': mass.Composition(formula='H-2O-1' + 'ON-1H-1' + 'N-1H-3'),
    'c-z':mass.Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1' + 'ON-1H-1'),
    'c-zdot': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1' + 'ON-1'),
    'c-z+1': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-2O-1' + 'ON-1H1'),
    'c-y': mass.Composition(formula='H-2O-1' + 'NH3' + ''),
    'cdot-y': mass.Composition(formula='H-2O-1' + 'NH3' + 'H-1' + ''),
    'adot-z': mass.Composition(formula='H-2O-1' + 'C-1O-1' + '1H' + 'H-2O-1' + 'ON-1H-1'),
    'adot-zdot': mass.Composition(formula='H-2O-1' + 'C-1O-1' + '1H' + 'H-2O-1' + 'ON-1'),
    'adot-z+1': mass.Composition(formula='H-2O-1' + 'C-1O-1' + '1H' + 'H-2O-1' + 'ON-1H1'),
    'a-y': mass.Composition(formula='H-2O-1' + 'C-1O-1' + ''),
    'b-y': mass.Composition(formula='H-2O-1' + ''),
    'adot-y': mass.Composition(formula='H-2O-1' + 'C-1O-1'+'H1' + ''),
    'z+1-c-1': mass.Composition(formula='H-2O-1' + 'ON-1H1' + 'H-2O-1' + 'NH3' + 'H-1'),
    }

def compodict_to_elementaldict(compodict):
    """
    Function that helps to convert a composition dictionary to an elemental one (input for correct theoretical envelope).
    :param compodict:
    :return: dict
    """

    elem_comp_dict ={}
    for s in compodict:
        # print(s)
        # print(compodict[s])

        if 'H' in s:
            elem_comp_dict['H'] = compodict[s]
        if 'C[12]' in s:
            elem_comp_dict['C'] = compodict[s]
        if 'C[13]' in s:
            elem_comp_dict['C'] += compodict[s]
        if 'N[14]' in s:
            elem_comp_dict['N'] = compodict[s]
        if 'N[15]' in s:
            elem_comp_dict['N'] += compodict[s]
        elif 'O' in s:
            elem_comp_dict['O'] = compodict[s]
        elif 'S' in s:
            elem_comp_dict['S'] = compodict[s]
        elif 'Se' in s:
            elem_comp_dict['Se'] = compodict[s]
        elif 'Fe' in s:
            elem_comp_dict['Fe'] = compodict[s]
    return elem_comp_dict


def compositionobj_to_dict(pyteo_comp_obj, ss_num, inter_mods, mods, iontype):
    """
    Function to obtain true elemental composition
    :param pyteo_comp_obj: Pyteomics Composiiton object based on the internal fragment sequence alone
    :param ss_num: Number of disulfides in the internal fragment
    :param inter_mods: List of modification possible for a cystine that is forming a disulfide with another fragment
    :param mods: noncys mods
    :param iontype: str, ion tye
    :return: dict, with the tru elemental composition
    """

    reality_bool = True

    elemcompo = pyteo_comp_obj[0]
    seqbased_dict = compodict_to_elementaldict(elemcompo)


    # print(f"Before ion type calculations = {seqbased_dict}")

    #Consider Iontype
    ionelemcompo = iontypedict[iontype]
    iontypebased_dict = compodict_to_elementaldict(ionelemcompo)


    seqion_compo = [seqbased_dict,iontypebased_dict]

    #Aggregating results
    result = {k: sum(d[k] for d in seqion_compo if k in d) for k in set(k for d in seqion_compo for k in d)}


    # print(f"Before ss calculations = {result}")
    # print(f"ss_num = {ss_num}")


    #Disulfide bonds
    disulfide_dict = {}
    if ss_num > 0:
        disulfide_dict = {"H": -2*ss_num}

    seqiondisul_compo = [seqbased_dict, iontypebased_dict,disulfide_dict]

    result = {k: sum(d[k] for d in seqiondisul_compo if k in d) for k in set(k for d in seqiondisul_compo for k in d)}

    # print(f"Before cys mods calculations = {result}")

    modcompo_dict = {'H': 0, 'C': 0, 'S': 0, 'O': 0, 'N': 0, 'Fe':0, "Se":0}
    #If there are disulfide bonds modifications to be considered
    if inter_mods:

        # print(inter_mods)
        # print(type(inter_mods))
        inter_mods = inter_mods.lstrip('(')
        inter_mods = inter_mods.strip(')')
        # print(inter_mods)
        inter_mods_spl = inter_mods.split(',')
        # print(inter_mods)


        for modss in inter_mods_spl:
            # print(f'Loss = {loss}')
            if modss:

                # print(f"{loss + '' == 'shl'}, {loss == 'sshl'}, {loss == 'chhsshl'}, {loss == 'hl'}")
                # print(f"{len(loss)}")
                # print(len('chhsshl'))

                # Must strip losses of quotation marks, otherwise they are 2 chars longer!!! Due to the combination fucntion
                modss = modss.strip(" ")
                modss = modss.lstrip(" ")
                modss = modss.strip("'")
                modss = modss.lstrip("'")

                # print(f"{len(loss)}")
                # print(f"{loss + '' == 'shl'}, {loss == 'sshl'}, {loss == 'chhsshl'}, {loss == 'hl'}")

                if modss == 'shl':
                    # print(loss == 'shl')
                    modcompo_dict["H"] += -1
                    modcompo_dict["S"] += -1
                elif modss == 'sshl':
                    # print(loss == "sshl")
                    modcompo_dict["H"] += -1
                    modcompo_dict["S"] += -2
                elif modss == 'chhsshl':
                    # print(loss == 'chhsshl')
                    modcompo_dict["H"] += -3
                    modcompo_dict["S"] += -2
                    modcompo_dict["C"] += -1
                elif modss == 'hl':
                    # print(loss == 'hl')
                    modcompo_dict["H"] += -1
                elif modss == 'h':
                    # print(loss == 'hl')
                    modcompo_dict["H"] += 1
                elif modss == 'sh':
                    # print(loss == 'hl')
                    modcompo_dict["H"] += 1
                    modcompo_dict["S"] += 1
                elif modss == "oxyhemeChl":
                    modcompo_dict["C"] += 34
                    modcompo_dict["H"] += 33
                    modcompo_dict["Fe"] += 1
                    modcompo_dict["N"] += 4
                    modcompo_dict["O"] += 4
                    modcompo_dict["S"] += 2
                elif modss == "semioxyhemeChl" or modss == "oxyhemeC":
                    modcompo_dict["C"] += 34
                    modcompo_dict["H"] += 34
                    modcompo_dict["Fe"] += 1
                    modcompo_dict["N"] += 4
                    modcompo_dict["O"] += 4
                    modcompo_dict["S"] += 2

            else:
                continue
    # print(f"Before mods calculations = {modcompo_dict}")
    # If there are modifications to be considered
    if len(mods) > 0:
        for mod in mods:
            print(f"mod = {mod}")
            if type(mod) == list:
                mod_split = mod[0].split("_")
                print(f"mod_split = {mod_split}")
                mods_elemdict = mods_repo[mod_split[0]].elemcomp
                for elem in mods_elemdict:
                    # print(elem)
                    # print(type(elem))
                    # print(mods_elemdict [elem])
                    # print(type(mods_elemdict[elem]))
                    modcompo_dict[elem] += mods_elemdict[elem] * int(mod_split[1])
            else:
                # print(f"{mod} composition is: {mods_repo[mod].elemcomp}")
                mods_elemdict = mods_repo[mod].elemcomp
                for elem in mods_elemdict:
                    # print(elem)
                    # print(type(elem))
                    # print(mods_elemdict [elem])
                    # print(type(mods_elemdict[elem]))
                    modcompo_dict[elem] += mods_elemdict[elem]

    # print(f"Inside compositionobj function = {modcompo_dict}")

    finalls_elemcompo = [seqbased_dict,iontypebased_dict,disulfide_dict,modcompo_dict]
    finalresult= {k: sum(d[k] for d in finalls_elemcompo  if k in d) for k in set(k for d in finalls_elemcompo  for k in d)}

    # Handeling negative element amounts
    for elem in finalresult:
        print(f"finalresult= {finalresult[elem]}")
        if finalresult[elem] < 0:
            reality_bool = False
            continue
    print(f"finalls_elemcompo  = {finalresult} and the reality bool is {reality_bool}")


    return finalresult, reality_bool

def elem_dict_to_isotopic_env(elem_dict, isotopologues_num= 15, charge = 4, norm_int = None, error_offset=None, resolution_param=17000):
    """
    :param elem_dict: dict, elemental composition of the internal fragment considering losses, ion types, modifications and disulfides
    :param isotopologues_num: The number of theoretical isotopologues to be calculated
    :param charge: int, charge of the experimental ion
    :param norm_int: to produced a theoretical isotope envelope with normalized intensity
    :param error_offset: The error offset between the experimental and the theoretical
    :return: m/z array and intensity array of the theoretical isotopic envelope
    """

    #Adding hydrogens to chemical composition produced by pyteomics,
    # in order to obtained right mz values when calculating envelope with pythoms
    #masses will be 0.0005 off still. Reference value of the mass of an electron: 0.000548 Da
    elem_dict["H"] += int(1*charge)

    print(f"elem_dict inside elem_dict_to_isotopic_env function = {elem_dict}\nCreating theoretical isotop envelopes with {resolution_param} mass resolving power\n")

    #From elem_dict to string
    outstr = ''
    for key in elem_dict.keys():
        if elem_dict[key] == 0:
            pass
        else:
            outstr += f"{key}{elem_dict[key]}"

    # print(f"outstr = {outstr}/ type outstr = {type(outstr)}/charge = {charge}")

    #Using pythoms isotope envelope
    mol = pythmole.IPMolecule(outstr, charge= int(charge), resolution = float(resolution_param))
    pythmole_mz = mol.bar_isotope_pattern[0]
    pythmole_int = mol.bar_isotope_pattern[1]

    # print(f"pythmole_mz: {pythmole_mz} ")

    mz_array = np.asarray(pythmole_mz)
    int_array = np.asarray(pythmole_int)

    # print(f"error_offset = {error_offset}")
    if error_offset:
        mz_array += error_offset
        #Correcting for adding Hydrogens and not protons
        mz_array -= charge*0.0005

    if norm_int:
        int_array = (int_array / int_array.max())

    print(f"Theoretical mz array after offset: {mz_array} ")

    return mz_array, int_array

def compare_isoenv(expmz ,expint, theomz, theoint, title, charge, error=None):
    """
    All are arrays. Function to compare between the experimental and the theoretical isotope envelopes
    :param theomz: m/z array of the theoretical ion
    :param theoint: intensity array of the theoretical ion
    :param expmz: m/z array of the experimental ion
    :param expint:  intensity array of the experimental ion
    :return:
    """

    # print(f"expmz = {expmz}")
    # print(f"theomz = {theomz}")
    mzscore = 0
    intscore = 0
    compound_score = 0

    #Find exp_mz maxima
    maxmz_ls = []
    maxint_ls = []

    # Normalized exp data
    try:
        maxexpint = expint.max()
    except ValueError:
        maxexpint = 1

    # print(f"maxexpint = {maxexpint}")

    exp_int_array_norm = (expint / maxexpint)
    # print(f"Max absolute int: {maxexpint}")


    #Oct 20, 2021 - Great parameters for first pass
    h = 0.4
    p = 0.3
    wl = 10

    foo = scipy.signal.find_peaks(exp_int_array_norm, height=h, prominence=p, wlen = wl)
    # print(f"foo[0] = {foo[0]}")
    # print(f"len of theo = {len(theomz)}")

    #If signal was not so good, less stringent parameters are needed
    if len(foo[0]) < 4:
        optifoo = scipy.signal.find_peaks(exp_int_array_norm, height=h - 0.2, prominence=p - 0.2, wlen= wl + 5)
        # print(f"optifoo[0] = {optifoo[0]}")
        # print(f"length optifoo[0] = {len(optifoo[0])}")
        # if len(optifoo[0]) <= len(theomz):
        #     foo = optifoo
        foo = optifoo


    #Create a list to hold exp_ion maxima
    for index in foo[0]:
        maxmz_ls.append(expmz[index])
        maxint_ls.append(exp_int_array_norm[index])

    orgmzls = maxmz_ls
    orgintls = maxint_ls
    print(maxmz_ls)



    mzscore = 0
    intscore = 0
    compound_score = 0


    #Scoring by isotopologues
    print("Scoring...")


    dictvals = {}
    for val in maxmz_ls[:5]:
        for theoval in theomz[:5]:
            # print(f"val = {val}")
            diff = abs(val-theoval)
            # print(f"diff = {diff}")

            dictvals[f"{val},{theoval}"] = diff


    print(f"dictvals = {dictvals}")

    #List containing the exp ion mz value and theoretical mz value with the smallest diff
    # initialmzisotopelogue = min(dictvals, key=dictvals.get)

    initialmzisotopelogue = {key: val for key, val in dictvals.items() if val == min(dictvals.values())}



    initialmzisotopeloguestr = ""
    for val in initialmzisotopelogue:
        initialmzisotopeloguestr = val


    isosplits = initialmzisotopeloguestr.split(",")

    print(f"{isosplits} = isosplits")

    maxmz_ls = np.asarray(maxmz_ls)

    expham = np.where(maxmz_ls == float(isosplits[0]))

    theoham = np.where(theomz == float(isosplits[1]))

    print(f"{expham} and {theoham}")

    truncmaxmz_ls = maxmz_ls[expham[0][0]:]
    truncmaxint_ls = maxint_ls[expham[0][0]:]

    truncmaxmz_array = np.asarray(truncmaxmz_ls)
    print(truncmaxmz_array)
    truncmaxint_array = np.asarray(truncmaxint_ls)

    lenmzlimit = len(truncmaxmz_array)
    print(f"lenmzlimit = {lenmzlimit}")
    lenintlimit = len(truncmaxint_array)
    print(f"lenintlimit = {lenintlimit}")

    trunctheomz = theomz[theoham[0][0]:]
    trunctheoint = theoint[theoham[0][0]:]

    lenmzlimittheo = len(theomz[theoham[0][0]:])
    lenintlimittheo = len(theoint[theoham[0][0]:])

    print(f"lenmzlimittheo = {len(theomz[theoham[0][0]:])}")
    print(f"lenintlimittheo = {len(theoint[theoham[0][0]:])}")

    if lenmzlimit < lenmzlimittheo:
        trunctheomz = trunctheomz[:lenmzlimit]
        trunctheoint = trunctheoint[:lenintlimit]

    elif lenmzlimit > lenmzlimittheo:
        truncmaxmz_array = truncmaxmz_array[:lenmzlimittheo]
        truncmaxint_array = truncmaxint_array[:lenintlimittheo]



    # trunctheomz = theomz[theoham[0][0]:theoham[0][0]+lenmzlimit]
    # print(trunctheomz)
    # trunctheoint = theoint[theoham[0][0]:theoham[0][0]+lenintlimit]

    # corr_mzscore = np.corrcoef(truncmaxmz_array,trunctheomz)
    # corre_mzscore = np.correlate(truncmaxmz_array,trunctheomz)
    #
    # corr_intscore = np.corrcoef(truncmaxint_array, trunctheoint)
    # corre_intscore = np.correlate(truncmaxint_array, trunctheoint)
    #
    # print(f"the scores = {corre_mzscore,corr_mzscore, corre_intscore, corr_intscore}")

    diff_mzarr = truncmaxmz_array - trunctheomz
    mzscore = diff_mzarr.max()
    print(f"diff_mzarr={diff_mzarr}, {diff_mzarr.max()}")

    diff_intarr = truncmaxint_array - trunctheoint
    intscore = diff_intarr.max()
    print(f"diff_mzarr={diff_intarr}, {diff_intarr.max()}")

    compound_score = mzscore+intscore

    # Diagnostic Plotting
    plt.figure('overalyn', dpi=300)
    theointzero = np.zeros(len(theomz))
    expintzero = np.zeros(len(truncmaxmz_array))

    plt.clf()
    plt.title(f"{round(error, 4)}, coor = {mzscore}")
    plt.scatter(theomz, theointzero, label="Theoretical", color="blue")
    plt.scatter(truncmaxmz_array, expintzero, label="Experimental", color="green")
    plt.xlabel("m/z")
    plt.ylabel("Relative intensity")
    plt.legend(loc='best')
    plt.plot()
    plt.savefig(title + "_comparescore.png")
    plt.close()
    plt.figure('overalyn', dpi=300)
    theointzero = np.zeros(len(theomz))
    expintzero = np.zeros(len(truncmaxint_array))
    #
    plt.clf()
    plt.title(f"{round(error, 4)}, coor = {intscore}")
    plt.scatter(theointzero, theoint, label="Theoretical", color="orange")
    plt.scatter(expintzero, truncmaxint_array, label="Experimental", color="magenta")
    #
    plt.xlabel("m/z")
    plt.ylabel("Relative intensity")
    plt.legend(loc='best')
    plt.plot()
    plt.savefig(title + "_comparescore_int.png")
    plt.close()
    #

    return truncmaxmz_array, truncmaxint_array, list(diff_mzarr), list(diff_intarr), exp_int_array_norm, compound_score, orgmzls, orgintls


def matchmaker2_multipass(theo_dict_tuple, exp_ls, mass_res, ppm_error, fullprotein_seq=None):
    """
    Function to match experimental to theoretical ions, when doing multipass searching.
    :param theo_dict_tuple: tuple. First position: str, pass name. Second postiiton: dict, keys are neutral masses; values are internal fragment objects
    :param rev_theo_dict:tuple. First position: pass name. Second postiiton: dict, keys are neutral masses; values are internal fragment objects
    :param exp_ls: ls, experimental ion objects
    :param ppm_error: int, how much tolerance error to use
    :return: ls, matched experimental objects (in order to keep tabs on what does not need to be included in subsequent passes)
    """

    #Unpacking tuples
    analysisName = theo_dict_tuple[0]
    theo_dict = theo_dict_tuple[1]
    # analysisName_rev = rev_theo_dict_tuple[0]
    # rev_theo_dict = rev_theo_dict_tuple[1]
    #
    # if analysisName != analysisName_rev:
    #     print("Unreverse and reverse passes are not the same in the matchmaker!!!!")

    #Create header
    out_str = f"{fullprotein_seq}\nMass Resolution: {mass_res} \nTolerance error: +/- {ppm_error} ppm\nneutral exp_ion\tneutral theoretical ion\tseq\tcharge\tmz_mono\ttmods\tion_type\tcysteine_locations\tss_count\tcysteines-with-mods\tcysteine mods" \
              "\tStart AA\tEnd AA\treverse_bool\tIntensity\tcyclic_density\terror\tchemical_composition\tisomz_score\tisoint_score\tfragment_score\n"

    print("~~~~~~~ Matching Fragments~~~~~~~~~")

    # Set reverse sequence theoretical fragments amount as the same as the non-reverse theoretical fragment database
    # frag_dict_rev_mod = {}
    # print(f"Frag_dict = {len(theo_dict)}")
    # print(f"Frag_dict_rev = {len(rev_theo_dict)}")
    # if len(theo_dict) != len(rev_theo_dict):
    #     frag_dict_rev_keys = rev_theo_dict.keys()
    #     counter = 0
    #     for x in frag_dict_rev_keys:
    #         if counter < len(theo_dict) + 1:
    #             frag_dict_rev_mod[x] = rev_theo_dict[x]
    #             counter += 1
    # print(f"Frag_dict_rev_mod = {len(frag_dict_rev_mod)}")
    #

    #Extract theoretical keys for matching
    theo_keys = theo_dict.keys()
    # rev_keys = rev_theo_dict.keys()

    matched_ls = []

    # Using the experimental ion list because it is shorter
    for expobj in exp_ls:

        #Extract infromation from experimental ion objecta
        expmass = expobj.exp_neut
        expmz = expobj.exp_mz
        expz = expobj.charge
        expint = expobj.pkht_cluster
        exp_mz_array = np.asarray(expobj.mz_isoenv)
        exp_int_array = np.asarray(expobj.int_isoenv)

        if len(exp_mz_array) == 0:
            continue


        for x in theo_keys:
            # print(x)
            error = expmass-x
            mz_error = expmz - theo_dict[x].mz_mono
            error_ppm = (error/x)*1000000
            # print(f" error_ppm  = {error_ppm }")

            # print(f"\nTheoretical ion = {theo_dict[x]}//expmz = {expmz} with error {error_ppm}")
            #If error is within error tolerance

            if abs(error_ppm) < float(ppm_error) and expobj.charge == theo_dict[x].charge:

                # print(f"\nTheoretical ion = {theo_dict[x]}//expmz = {expmz} with error {ppm_error}")

                #Get the theoretical sequence
                sequence = theo_dict[x].sequence

                #Parameters to pass to make the correct elemental composition
                ss_num = theo_dict[x].ss_count
                interfragss_mods = theo_dict[x].cysmods
                interfrag_mods = theo_dict[x].mods
                iontype = theo_dict[x].ion_type
                # print(ss_num)
                # ss_hydrogenloss = f"H{ss_num}"

                # Cyclic density = disulfide brige regions per length of protein sequence
                cyclic_den = ss_num / len(sequence)

                #Create and modified elemental composition whne the theoretical ion is created, if negative element amounts are obtained don't crete the theoretical on

                # Calculating elemental composition only based on simple sequence
                theo_iso_env = mass.most_probable_isotopic_composition(sequence=sequence)

                #Get correct elemental composition
                elemcomp_dict, reality_bool = compositionobj_to_dict(theo_iso_env, ss_num, interfragss_mods, interfrag_mods, iontype)
                # print(elemcomp_dict)

                #If elemental composition had only positive amounts
                if reality_bool:

                    print("Inside REality Bool")

                    #Calculate isotopic envelope for theoretical fragment
                    theomz_array, theoint_array =elem_dict_to_isotopic_env(elemcomp_dict, charge = expobj.charge, norm_int=True, error_offset=mz_error,resolution_param=mass_res)

                    if theomz_array[0] - exp_mz_array[0] > 1.75:
                        continue
                    else:
                        matched_ls.append(expobj)

                        #Compare isotopic envelopes
                        expmass_str = str(round(expmass,2)).replace(".","_")
                        fig_title = f"exp-{expmass_str}"
                        sanity_checkI, sanity_checkII, sanity_checkIII, sanity_checkIV, exp_int_array_norm, sanity_checkV, expmaximamz, expmaximaint\
                            = compare_isoenv(exp_mz_array, exp_int_array,theomz_array, theoint_array,title = fig_title, error = mz_error, charge =  expz)


                        # Create comparison plot of theoretical and experimental isotopic envelopes
                        plt.figure('overalyn', dpi=300)

                        plt.clf()
                        plt.title(f"{round(error_ppm,4)} ppm")
                        plt.plot(expmaximamz, expmaximaint, color='green', label="Experimentalone")
                        plt.plot(exp_mz_array, exp_int_array_norm, color='orange', label = "Experimental")
                        plt.scatter(theomz_array, theoint_array, label = "Theoretical", marker = 'o')
                        plt.xlabel("m/z")
                        plt.ylabel("Relative intensity")
                        plt.legend(loc='best')
                        plt.plot()

                        #There are times an OS error arrises, specify full path as a solution
                        currentdir = os.getcwd()
                        fullouputname = os.path.join(currentdir, fig_title)
                        plt.savefig(fullouputname + ".png")
                        plt.close()

                        # if x in rev_keys:
                        #     out_str += f"{expmass}\t{theo_dict[x]}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"
                        #     out_str += f"{expmass}\t{rev_theo_dict[x]}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"
                        # else:
                        #     out_str += f"{expmass}\t{theo_dict[x]}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"

                        out_str += f"{expmass}\t{theo_dict[x]}\t{expint}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"
                else:
                    continue
        #Save results
        try:
            output = open(f"{analysisName}_hits" + '.tsv', 'w')
        except OSError:
            output = open(f"{analysisName}_hits" + '.tsv', 'a')

        output.write(out_str)
        output.close()

    return matched_ls


class Hit:
    """
    Wrapper object for matches found. Contains a theoretical fragment, an experimental cluster, and the error
    (in ppm) with which they were matched.
    """
    def __init__(self, thy_ion, exp_ion, error, pass_num):
        """
        :param thy_ion = theoretical Fragment object containing information about the predicted fragment
        :param exp_ion = experimental cluster detected and matched to the theoretical fragment
        :param error = error in ppm between mz of the theoretical and experimental monoisotopic peaks.
        :param pass_num = which pass this hit was found in the multipass search (typically 1, 2, or 3)
        :param cal_error = error after calibration
        """
        self.thy_ion = thy_ion
        self.exp_ion = exp_ion
        self.error = error
        self.pass_num = pass_num
        self.cal_error = 0

    def __eq__(self, other):
        return self.thy_ion.mz_mono == other.thy_ion.mz_mono


    def __hash__(self):
        # print(hash(str(self)))
        return hash(self.thy_ion)

    # def __hash__(self):
    #     print(hash(str(self)))
        # return hash((self.seq, self.seq_index, self.term))

    def __lt__(self, other):
        return self.thy_ion.mz_mono < other.thy_ion.mz_mono

    def print_hit_info(self):
        """
        Returns/prints all information from a hit to a single line in the same order as default header, comma separated
        :return:
        """
        # handle/convert information for printing, including checking for types (mods, losses) which may be None
        cys_loc_arg = ''
        if self.thy_ion.cysloc is not None:
            thy_modcl = str(self.thy_ion.cysloc)
            cys_loc_arg += thy_modcl.replace(',',';')

        cys_loc_arg = cys_loc_arg[1:-1]


        mod_arg = ''
        if self.thy_ion.thy_mods is not None:
            thy_modstr = str(self.thy_ion.thy_mods)
            mod_arg += thy_modstr.replace(',',';')

        mod_arg = mod_arg[1:-1]


        cys_mod_arg = ''

        if self.thy_ion.cysmods is not None:
            # print(f"Before processing it to printable forms: {self.thy_ion.cysmods}-{type(self.thy_ion.cysmods)}")
            cysmodstr =  str(self.thy_ion.cysmods)
            cysmodstr = cysmodstr[1:-1]
            cysmodssplit = cysmodstr.split(',')
            for cmod in cysmodssplit:
                if cmod == ' ':
                    continue
                else:
                    cys_mod_arg += f"{cmod};"
            cys_mod_arg.strip(';')
        cys_mod_arg = cys_mod_arg[:-1]


        # print(f"Disulfides mods = {cys_mod_arg}")

        ion_str = '({}){}'.format(self.thy_ion.ion_type, self.thy_ion.ion_type_indx)

        # format and write information to file
        try:

            line = f'{self.pass_num},{self.exp_ion.cal_mz_mono}, {self.thy_ion.mz_mono}, {self.cal_error}, {int(self.exp_ion.charge)}, {ion_str}, {mod_arg}, {self.thy_ion.neutlosses}, {self.thy_ion.mono_neutral}, {self.thy_ion.cys_num},{cys_loc_arg},{cys_mod_arg}'
        except AttributeError:
            # catch data that was not auto-calibrated

            line = f'{self.pass_num},{self.exp_ion.mz_mono}, {self.thy_ion.mz_mono}, {self.error}, {int(self.exp_ion.charge)}, {ion_str}, {mod_arg}, {self.thy_ion.neutlosses}, {self.thy_ion.mono_neutral},{self.thy_ion.cys_num}, {cys_loc_arg},{cys_mod_arg}'

        expline_args = ['{}'.format(x) for x in self.exp_ion.data_list[1:]]
        expline_args.extend([str(self.exp_ion.mz_mono), str(self.error)])
        expline = ','.join(expline_args)
        line = line + ',' + expline

        # print(line)
        return line

    def __str__(self):
        """
        string representation for debugging
        :return: string
        """
        if not self.cal_error == 0:
            return '<Hit> mz: {:.2f}, {:.1f} ppm'.format(self.exp_ion.mz_mono, self.cal_error)
        else:
            return '<Hit> mz: {:.2f}, {:.1f} ppm'.format(self.exp_ion.mz_mono, self.error)
    __repr__ = __str__

def clear_site_hits(siteDict):
    """
    Clear all hits from a list of FragmentSite containers in between analyses to prevent
    them from carrying over to the next analysis (but preserving the theoretical ions calculated
    for that site)
    :param site_list: list of FragmentSite containers
    :type site_list: list[FragmentSite]
    :return: updated site list with hits cleared
    """
    for site in siteDict:
        siteDict[site].hits = []
    return siteDict

def matchmaker_terminal_multipass(theo_dict_tuple, exp_ls, rev_theo_dict_tuple=None):
    """
    Function to match experimental to theoretical ions, when doing multipass searching.
    :param theo_dict_tuple: tuple. First position: str, pass name. Second postiiton: dict, sites-keys; FragSites-values
    :param rev_theo_dict:tuple. First position: pass name. Second postiiton: dict, keys are neutral masses; values are internal fragment objects
    :param exp_ls: ls, experimental ion objects
    :param ppm_error: int, how much tolerance error to use
    :return: ls, matched experimental objects (in order to keep tabs on what does not need to be included in subsequent passes)
    """

    #Unpacking tuples
    analysisName = theo_dict_tuple[0]
    site_dict = theo_dict_tuple[1]

    #To prevent carrying over hits to different passes or files
    site_dict = clear_site_hits(site_dict)

    init_tol = theo_dict_tuple[2]
    final_tol = theo_dict_tuple[3]
    cal_bool = theo_dict_tuple[4]
    # analysisName_rev = rev_theo_dict_tuple[0]
    # rev_theo_dict = rev_theo_dict_tuple[1]

    # if analysisName != analysisName_rev:
    #     print("Unreverse and reverse passes are not the same in the matchmaker!!!!")

    print("~~~~~~~ Matching Fragments~~~~~~~~~")

    # Set reverse sequence theoretical fragments amount as the same as the non-reverse theoretical fragment database
    # frag_dict_rev_mod = {}
    # print(f"Frag_dict = {len(theo_dict)}")
    # print(f"Frag_dict_rev = {len(rev_theo_dict)}")
    # if len(theo_dict) != len(rev_theo_dict):
    #     frag_dict_rev_keys = rev_theo_dict.keys()
    #     counter = 0
    #     for x in frag_dict_rev_keys:
    #         if counter < len(theo_dict) + 1:
    #             frag_dict_rev_mod[x] = rev_theo_dict[x]
    #             counter += 1
    # print(f"Frag_dict_rev_mod = {len(frag_dict_rev_mod)}")
    #

    #Extract theoretical keys for matching
    # keys = theo_dict.keys()
    # rev_keys = rev_theo_dict.keys()

    matched_ls = []

    # Using the experimental ion list because it is shorter
    for expobj in exp_ls:

        #Extract infromation from experimental ion objecta
        expmz = expobj.mz_mono
        charge = expobj.charge

        for site in site_dict:
            # print(site)
            # print(site_dict[site])

            for ion in site_dict[site].theo_ions:
                # print(ion)
                theoretical_ion = site_dict[site].theo_ions[ion]

                mz_error = ion - expmz
                error_ppm = (mz_error / ion) * 1000000

                # print(f"error_ppm = {error_ppm}")

                # Match possible if error is less than tolerance
                if abs(error_ppm) < init_tol:
                    # confirm that charge state matches theoretical
                    if charge == theoretical_ion.charge:
                        matched_ls.append(expobj)
                        matched_hit = Hit(theoretical_ion, expobj, error_ppm, analysisName)
                        # print(f"Matched Hit = {matched_hit}")
                        site_dict[site].hits.append(matched_hit)


    #If calibration is desired
    if cal_bool:
        all_sites, median_error, average_error, expion_restore, matched_ls_aftercal = calibrate_data(site_dict, final_tol, matchedls=matched_ls)

        for ion in expion_restore:
            if ion in exp_ls:
                continue
            else:
                exp_ls.append(ion)

        # for site in all_sites:
            # print(all_sites[site].hits)
    else:
        all_sites, median_error, average_error, expion_restore, matched_ls_aftercal = calibrate_data(site_dict, final_tol, matchedls=matched_ls, cal=False)


    return matched_ls_aftercal, all_sites, median_error, average_error



    # #Save results
    # output = open(f"{analysisName}_hits" + '.tsv', 'w')
    # output.write(out_str)
    # output.close()



def calibrate_data(dict_of_sites, final_tol, matchedls = None, cal = True):
    """
    Method to determine appropriate calibration values for data and apply it to the output. Assumes that
    a wide error (e.g. 100ppm) calibration has been applied initially and that a narrower calibration
    (to the final_tol parameter) should be applied around the median of the initial errors.
    :param all_sites: List of fragment site objects containing all annotated fragmt information (n-term)
    :type all_sites: list[FragmentSite]
    :param final_tol: final tolerance (ppm) to filter the initial hit lists
    :return: updated site list with errors calibrated and 'hits' outside the new (narrow) tolerance removed,
     median error, average error, exp_ions from hits that after calibtration had a large error
    """
    hits = []
    # find all sites with hits and make a list of hits
    for site in dict_of_sites:
        # standard auto-cal
        hits.extend(dict_of_sites[site].hits)
        # In order to sort needed to add __lt__  to Hit and Thy Ion objects
        hits.sort()

    #TODO: someday fit polynomials and stuff like that
    #Determine calibration as median of errors in annotated data if it hasn't been provided
    size = int(len(hits) / 3)

    errors1 = [hit.error for hit in hits[:size]]
    errors2 = [hit.error for hit in hits[size:(2 * size)]]
    errors3 = [hit.error for hit in hits[(2 * size):]]

    median_error1 = np.median(np.asarray(errors1))
    average_error1 = np.average(np.asarray(errors1))
    median_error2 = np.median(np.asarray(errors2))
    average_error2 = np.average(np.asarray(errors2))
    median_error3 = np.median(np.asarray(errors3))
    average_error3 = np.average(np.asarray(errors3))

    errors12 = errors1 + errors2
    errors = errors12+errors3
    median_error = np.median(np.asarray(errors))
    average_error = np.average(np.asarray(errors))

    # filter out any hits that don't meet the new calibration criteria and add a cal_exp_mono_m/z field to the others
    exp_ion_restorelist = []

    if cal:
        for site in dict_of_sites:
            hits_remove_list = []

            for hit in dict_of_sites[site].hits:
                if hit.error in errors1:
                    hit.cal_error = hit.error - median_error1
                elif hit.error in errors2:
                    hit.cal_error = hit.error - median_error2
                elif hit.error in errors3:
                    hit.cal_error = hit.error - median_error3

                if -final_tol < hit.cal_error < final_tol:
                    # The calibrated hit is within the narrow tolerance - update it to include a calibrated m/z
                    hit.exp_ion.cal_mz_mono = hit.exp_ion.mz_mono / (1 - median_error / 1000000)    # adjusted mz

                else:
                    # This 'hit' is not within final error - remove it from the site's hit list
                    matchedls.remove(hit.exp_ion)
                    exp_ion_restorelist.append(hit.exp_ion)
                    hits_remove_list.append(hit)

            for remove_hit in hits_remove_list:
                dict_of_sites[site].hits.remove(remove_hit)
    else:
        pass


    return dict_of_sites, median_error, average_error, exp_ion_restorelist, matchedls
