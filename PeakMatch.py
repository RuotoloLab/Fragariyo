"""
Author: Carolina Rojas Ramirez
Date: May 13, 2020
Python module with all the function needed for Peak Matching using the Fragmentor and the terminalFragmentor
"""
import scipy
import numpy as np
from scipy import signal
from pyteomics import mass
from matplotlib import pyplot as plt
from brainpy import isotopic_variants
import pythoms.molecule as pythmole
from Modifications import mods_repo


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


    elem_comp_dict = {}
    elem_ls = ['C', 'H', 'N', 'O','S','Fe']

    for element in elem_ls:
        elem_comp_dict[element] = 0

    print(f"Original composition = {pyteo_comp_obj[0]}")
    compostr = str(pyteo_comp_obj[0])
    compostr = compostr.lstrip('Composition(')
    compostr = compostr.strip(')')
    compostr = compostr.lstrip('{')
    compostr = compostr.strip('}')

    composplt = compostr.split(',')

    for x in composplt:
        # print(x)
        # print(type(x))
        x_spl = x.split(":")
        # print(x_spl[0])
        # print(x_spl[1])
        if 'H' in x_spl[0]:
            elem_comp_dict['H'] = int(x_spl[1].lstrip(' '))

        #Why did I wrote this? I don't remember, but I bet it corrects for a weird bug!
        #Nov 19th, 2021 = I remember it corrects for the different carbon and nitrogen isotopes calculated by pythonms
        elif 'C' in x_spl[0]:
        #If a second isotope is found added to the already key in the dict, don't replace the previous value
            if elem_comp_dict['C']:
                elem_comp_dict['C'] += int(x_spl[1].lstrip(' '))
            else:
                elem_comp_dict['C'] += int(x_spl[1].lstrip(' '))

        elif 'N' in x_spl[0]:
            if elem_comp_dict['N']:
                elem_comp_dict['N'] += int(x_spl[1].lstrip(' '))
            else:
                elem_comp_dict['N'] += int(x_spl[1].lstrip(' '))

        elif 'O' in x_spl[0]:
            elem_comp_dict['O'] = int(x_spl[1].lstrip(' '))
        elif 'S' in x_spl[0]:
            elem_comp_dict['S'] = int(x_spl[1].lstrip(' '))
        elif 'Fe' in x_spl[0]:
            elem_comp_dict['Fe'] = int(x_spl[1].lstrip(' '))

    print(f"Before ion type calculations = {elem_comp_dict}")
    #Consider Iontype
    if iontype == 'c-z':
        elem_comp_dict["H"] += -1
        elem_comp_dict["O"] += -1
    elif iontype == 'c-zdot':
        elem_comp_dict["H"] += -2
        elem_comp_dict["O"] += -1
    elif iontype == 'c-y':
        elem_comp_dict["H"] += 1
        elem_comp_dict["O"] += -1
        elem_comp_dict["N"] += 1
    elif iontype == 'cdot-y':
        elem_comp_dict["O"] += -1
        elem_comp_dict["N"] += 1
    elif iontype == 'a-z':
        elem_comp_dict["H"] += -4
        elem_comp_dict["O"] += -2
        elem_comp_dict["C"] += -1
        elem_comp_dict["N"] += -1
    elif iontype == 'a-zdot':
        elem_comp_dict["H"] += -5
        elem_comp_dict["O"] += -2
        elem_comp_dict["C"] += -1
        elem_comp_dict["N"] += -1
    elif iontype == 'b-y':
        elem_comp_dict["H"] += -2
        elem_comp_dict["O"] += -1

    elif iontype == 'a-y':
        elem_comp_dict["H"] += -2
        elem_comp_dict["O"] += -2
        elem_comp_dict["C"] += -1
    elif iontype == 'x-c':
        elem_comp_dict["H"] += -1
        elem_comp_dict["N"] += 1
        elem_comp_dict["C"] += 1
    elif iontype == 'x-cdot':
        elem_comp_dict["H"] += -2
        elem_comp_dict["N"] += 1
        elem_comp_dict["C"] += 1





    print(f"ss_num = {ss_num}")
    print(f"Before ss calculations = {elem_comp_dict}")
    #Disulfide bonds
    if ss_num > 0:
        elem_comp_dict["H"] -= 2*ss_num

    print(f"Before cys mods calculations = {elem_comp_dict}")
    #If there are disulfide bonds modifications to be considered
    if inter_mods:

        print(inter_mods)
        print(type(inter_mods))
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
                    elem_comp_dict["H"] += -1
                    elem_comp_dict["S"] += -1
                elif modss == 'sshl':
                    # print(loss == "sshl")
                    elem_comp_dict["H"] += -1
                    elem_comp_dict["S"] += -2
                elif modss == 'chhsshl':
                    # print(loss == 'chhsshl')
                    elem_comp_dict["H"] += -3
                    elem_comp_dict["S"] += -2
                    elem_comp_dict["C"] += -1
                elif modss == 'hl':
                    # print(loss == 'hl')
                    elem_comp_dict["H"] += -1
                elif modss == 'h':
                    # print(loss == 'hl')
                    elem_comp_dict["H"] += 1
                elif modss == 'sh':
                    # print(loss == 'hl')
                    elem_comp_dict["H"] += 1
                    elem_comp_dict["S"] += 1
                elif modss == "oxyhemeChl":
                    elem_comp_dict["C"] += 34
                    elem_comp_dict["H"] += 33
                    elem_comp_dict["Fe"] += 1
                    elem_comp_dict["N"] += 4
                    elem_comp_dict["O"] += 4
                    elem_comp_dict["S"] += 2
                elif modss == "semioxyhemeChl" or modss == "oxyhemeC":
                    elem_comp_dict["C"] += 34
                    elem_comp_dict["H"] += 34
                    elem_comp_dict["Fe"] += 1
                    elem_comp_dict["N"] += 4
                    elem_comp_dict["O"] += 4
                    elem_comp_dict["S"] += 2

            else:
                continue
    print(f"Before mods calculations = {elem_comp_dict}")
    # If there are modifications to be considered
    if len(mods) > 0:
        for mod in mods:
            print(mod)
            if type(mod) == list:
                mod_split = mod[0].split("_")
                print(mod_split)
                mods_elemdict = mods_repo[mod_split[0]].elemcomp
                for elem in mods_elemdict:
                    # print(elem)
                    # print(type(elem))
                    # print(mods_elemdict [elem])
                    # print(type(mods_elemdict[elem]))
                    elem_comp_dict[elem] += mods_elemdict[elem] * int(mod_split[1])
            else:
                # print(f"{mod} composition is: {mods_repo[mod].elemcomp}")
                mods_elemdict = mods_repo[mod].elemcomp
                for elem in mods_elemdict:
                    # print(elem)
                    # print(type(elem))
                    # print(mods_elemdict [elem])
                    # print(type(mods_elemdict[elem]))
                    elem_comp_dict[elem] += mods_elemdict[elem]

    print(f"Inside compositionobj function = {elem_comp_dict}")
    #Handeling negative element amounts
    for elem in elem_comp_dict:
        print(f"elem_comp_dict[elem] = {elem_comp_dict[elem]}")
        if elem_comp_dict[elem] < 0:
            elem_comp_dict[elem] = 0

    return elem_comp_dict

def elem_dict_to_isotopic_env(elem_dict, isotopologues_num= 15, charge = 4, norm_int = None, error_offset=None):
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
    #masses will be 0.0005 off still
    elem_dict["H"] += int(1*charge)

    print(f"elem_dict inside elem_dict_to_isotopic_env function = {elem_dict}")

    #From elem_dict to string
    outstr = ''
    for key in elem_dict.keys():
        if elem_dict[key] == 0:
            pass
        else:
            outstr += f"{key}{elem_dict[key]}"

    # print(f"outstr = {outstr}/ type outstr = {type(outstr)}/charge = {charge}")

    #Using pythoms isotope envelope
    mol = pythmole.IPMolecule(outstr, charge= int(charge), resolution = 17000)
    pythmole_mz = mol.bar_isotope_pattern[0]
    pythmole_int = mol.bar_isotope_pattern[1]

    # print(f"pythmole_mz: {pythmole_mz} ")

    mz_array = np.asarray(pythmole_mz)
    int_array = np.asarray(pythmole_int)

    # print(f"error_offset = {error_offset}")
    if error_offset:
        mz_array += error_offset

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
    h = 0.5
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


    # print(maxmz_ls)
    # print(f"theomz = {theomz}")
    # print(maxmz_ls)

    # indx = 0
    # if theomz[indx] - maxmz_ls[indx] > 0.1:
    #     print("Move to next value in maxmz_ls")


    mzscore = 0
    intscore = 0
    compound_score = 0


    #Scoring by isotopologues
    print("Scoring...")
    for val in maxmz_ls:
        diff = abs(val-theomz[0])
        print(f"diff = {diff}")
        # print(f"val = {val}")
        if diff < 0.4:
            ham = np.where(maxmz_ls == val)
            # print(f"Where to begin: {ham}")
            # print(ham[0])

            maxmz_ls = maxmz_ls[ham[0][0]:]
            maxint_ls = maxint_ls[ham[0][0]:]

            mz_int_zip = zip(maxmz_ls, maxint_ls)
            mz_int_ls = list(mz_int_zip)
            print(f"mz_int_ls = {mz_int_ls}")
            new_ls = []
            prev = 0
            for item in mz_int_ls:
                mzval = item[0]
                # print(item[0])
                if prev == 0:
                    prev = mzval
                    new_ls.append(item)
                else:
                    diff = 0
                    diff = round(abs(prev - mzval),1)
                    print(f"{prev} - {mzval} = {diff}")

                    if diff == (1/charge):
                        new_ls.append(item)
                        prev = mzval
                    else:
                        continue
            print(f"charge = {charge}")
            print(f"after charge correction = {new_ls}")

            new_ls_zip = zip(*new_ls)
            two_ls = list(new_ls_zip)
            # print(f"two_ls = {two_ls}")
            # print(two_ls[0])
            # print(type(two_ls[0]))

            maxmz_array = np.asarray(two_ls[0])
            maxint_array = np.asarray(two_ls[1])

            print(f"len(maxmz_array) = {len(maxmz_array)}")
            print(f"len(maxint_array) = {len(maxint_array)}")

            if len(maxmz_array) >= 2:

                trunc_mzarr = maxmz_array
                trunc_intarr = maxint_array

                print(f"trunc_mzarr = {trunc_mzarr}")
                print(f"trunc_intarr = {trunc_intarr}")

                if len(trunc_mzarr) > len(theomz):

                    print("mz > int")

                    lendiff = len(theomz) - len(trunc_mzarr)
                    trunc_mzarr = trunc_mzarr[:lendiff]
                    trunc_intarr = trunc_intarr[:lendiff]

                    print(len(trunc_mzarr), len(theomz))


                elif len(trunc_mzarr) < len(theomz):

                    print("mz < int")

                    lendiff = len(trunc_mzarr) - len(theomz)
                    theomz = theomz[:lendiff]
                    theoint = theoint[:lendiff]

                    print(len(trunc_mzarr), len(theomz))

                elif len(trunc_mzarr) == len(theomz):

                    pass

                # corr_mzscore = np.corrcoef(trunc_mzarr, theomz)
                # corre_mzscore = np.correlate(trunc_mzarr, theomz)
                # cov_mzscore = np.cov(trunc_mzarr, theomz)
                diff_mzarr = trunc_mzarr - theomz
                print(f"Diff = {trunc_mzarr - theomz}")
                notmached = 0
                for item in diff_mzarr:
                    if abs(item) > 0.05:
                        notmached += 1
                matchedpoints = len(diff_mzarr) - notmached

                try:
                    mzscore = (matchedpoints / (len(diff_mzarr))) * 100
                except ZeroDivisionError:
                    mzscore = 0

                    # print(f"Let's try numpy's Pearson Product correlation  (mz): {corr_mzscore}")
                # # print(f"Let's try numpy's correlate (mz): {corre_mzscore}")
                # # print(f"Let's try numpy's correlate (mz): {cov_mzscore}")
                # print(f"Std. Dev (mz): {np.std(trunc_mzarr)}")
                # mzscore = corr_mzscore[0][1]

                # After determining where to began comparison based on mz dimension, let's use numpy's correalte to compare intensity
                # corr_intscore = np.corrcoef(trunc_intarr, theoint)
                # print(f"Let's try numpy's Pearson Product correlation: {corr_intscore[0][1]}")
                # print(f"Std. Dev (mz): {np.std(trunc_intarr)}")
                # intscore = corr_intscore[0][1]

                diff_intarr = trunc_intarr - theoint
                print(f"Int Diff = {trunc_intarr - theoint}")
                notmached_int = 0
                for item in diff_intarr:
                    if abs(item) > 0.01:
                        notmached_int += 1
                matchedpoints_int = len(diff_intarr) - notmached_int
                try:
                    intscore = (matchedpoints_int / (len(diff_intarr))) * 100
                except ZeroDivisionError:
                    intscore = 0

                if mzscore == 0:

                    compound_score = 0

                else:
                    compound_score = round((mzscore + abs(intscore)) / 2, 3)

                # Diagnostic Plotting
                plt.figure('overalyn', dpi=300)
                theointzero = np.zeros(len(theomz))
                expintzero = np.zeros(len(trunc_mzarr))

                plt.clf()
                plt.title(f"{round(error, 4)}, coor = {mzscore}")
                plt.scatter(theomz, theointzero, label="Theoretical", color="blue")
                plt.scatter(trunc_mzarr, expintzero, label="Experimental", color="green")

                plt.xlabel("m/z")
                plt.ylabel("Relative intensity")
                plt.legend(loc='best')
                plt.plot()
                plt.savefig(title + "_comparescore.png")
                plt.close()

                plt.figure('overalyn', dpi=300)
                theointzero = np.zeros(len(theomz))
                expintzero = np.zeros(len(trunc_intarr))

                plt.clf()
                plt.title(f"{round(error, 4)}, coor = {intscore}")
                plt.scatter(theointzero, theoint, label="Theoretical", color="orange")
                plt.scatter(expintzero, trunc_intarr, label="Experimental", color="magenta")

                plt.xlabel("m/z")
                plt.ylabel("Relative intensity")
                plt.legend(loc='best')
                plt.plot()
                plt.savefig(title + "_comparescore_int.png")
                plt.close()

            else:
                continue



        else:
            # print("No possible to compare isotopologues")
            # mzscore = 0
            # intscore = 0
            # compound_score = 0
            continue






    return maxmz_ls, maxint_ls, mzscore, intscore, exp_int_array_norm, compound_score


def matchmaker2_multipass(theo_dict_tuple, exp_ls, ppm_error, fullprotein_seq=None):
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
    out_str = f"{fullprotein_seq}\nneutral exp_ion\tneutral theoretical ion\tseq\tcharge\tmz_mono\ttmods\tion_type\tcysteine locations\tss_count\tcysteines-with-mods\tcysteine mods" \
              "\tStart AA\tEnd AA\treverse_bool\tIntensity\tcyclic density\terror\tchemical_composition\tisomz_score\tisoint_score\tfragment_score\n"

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
            # print(f" error = {error}")


            #If error is within error tolerance
            if abs(error_ppm) < float(ppm_error) and expobj.charge == theo_dict[x].charge:

                print(f"\nTheoretical ion = {theo_dict[x]}//expmz = {expmz} wit error {ppm_error}")

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

                #TODO: Create and modified elemental composition whne the theoretical ion is created, if negative element amounts are obtained don't crete the theoretical on

                # Calculating elemental composition only based on simple sequence
                theo_iso_env = mass.most_probable_isotopic_composition(sequence=sequence)

                #Get correct elemental composition
                elemcomp_dict = compositionobj_to_dict(theo_iso_env, ss_num, interfragss_mods, interfrag_mods, iontype)
                # print(elemcomp_dict)

                #Calculate isotopic envelope for theoretical fragment
                theomz_array, theoint_array =elem_dict_to_isotopic_env(elemcomp_dict, charge = expobj.charge, norm_int=True, error_offset=mz_error)

                if theomz_array[0] - exp_mz_array[0] > 1.75:
                    continue
                else:
                    matched_ls.append(expobj)

                    #Compare isotopic envelopes
                    expmass_str = str(round(expmass,2)).replace(".","-")
                    fig_title = f"exp-{expmass_str}"
                    sanity_checkI, sanity_checkII, sanity_checkIII, sanity_checkIV, exp_int_array_norm, sanity_checkV = compare_isoenv(exp_mz_array, exp_int_array,
                                                                                    theomz_array, theoint_array,title = fig_title, error = mz_error, charge =  expz)


                    # Create comparison plot of theoretical and experimental isotopic envelopes
                    plt.figure('overalyn', dpi=300)

                    plt.clf()
                    plt.title(f"{round(error_ppm,4)} ppm")
                    plt.plot(exp_mz_array, exp_int_array_norm, color='orange', label = "Experimental")
                    plt.scatter(theomz_array, theoint_array, label = "Theoretical", marker = 'o')
                    plt.xlabel("m/z")
                    plt.ylabel("Relative intensity")
                    plt.legend(loc='best')
                    plt.plot()
                    plt.savefig(fig_title +".png")
                    plt.close()

                    # if x in rev_keys:
                    #     out_str += f"{expmass}\t{theo_dict[x]}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"
                    #     out_str += f"{expmass}\t{rev_theo_dict[x]}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"
                    # else:
                    #     out_str += f"{expmass}\t{theo_dict[x]}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"

                    out_str += f"{expmass}\t{theo_dict[x]}\t{expint}\t{cyclic_den}\t{error_ppm}\t{elemcomp_dict}\t{sanity_checkIII}\t{sanity_checkIV}\t{sanity_checkV}\n"

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
