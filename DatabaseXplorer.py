"""
Author: Carolina Rojas Ramirez
Date: 11/09/2021
"""
from tkinter import filedialog
import pickle
from terminalFragmentor_Main import FragmentSite
from terminalFragmentor_Main import ThyIon
from InternalFragmentor_Main import interfrag


def internal_database():
    """
    Function to explore .ions files that contain a terminal fragment database
    :return:
    """
    ions_file = filedialog.askopenfilename(title='Choose Ions File', filetypes=[('Theoretical Fragments', '.theofrags')])

    with open(ions_file, 'rb') as picklefile:
        saved_database = pickle.load(picklefile)

        for analysis_num in saved_database:
            count = 0
            print(f"The analysis is {analysis_num} and the protein sequence is {saved_database[analysis_num][0]}")
            print(f"The amount of passes are :{len(saved_database[analysis_num][1])}")

            passes = saved_database[analysis_num][1]

            for pass_label in passes:
                print(pass_label)
                for item in passes[pass_label]:
                    # print(item)
                    # print(len(item))
                    count+=len(item)

                    # for ion in item:
                    #     print(ion)

            print(count)




def terminal_database():
    """
    Function to explore .ions files that contain a terminal fragment database
    :return:
    """
    ions_file = filedialog.askopenfilename(title='Choose Ions File', filetypes=[('Ions File', '.ions')])
    with open(ions_file, 'rb') as picklefile:
        saved_database = pickle.load(picklefile)
        protein_seq = saved_database[0]
        analysis_dict = saved_database[1]

        # print(analysis_dict)

        count_ions = 0
        #Analysis
        for analysis in analysis_dict:
            # print(analysis)
            # print(analysis_dict[analysis])
            for analysis_pass in analysis_dict[analysis]:
                # print(analysis_pass)
                # print(analysis_dict[analysis][analysis_pass])
                for site in analysis_dict[analysis][analysis_pass][1]:
                    # print(site)
                    fragsite = analysis_dict[analysis][analysis_pass][1][site]
                    # print(fragsite)
                    # print(fragsite.theo_ions)
                    for ion in fragsite.theo_ions:

                        # print(f"{fragsite.theo_ions[ion].mz_mono}_{fragsite.theo_ions[ion].charge}_{fragsite.theo_ions[ion].ion_type}{fragsite.theo_ions[ion].ion_type_indx}_{fragsite.theo_ions[ion].thy_mods}")
                        count_ions += 1

        print(count_ions)


if __name__ == '__main__':

    # terminal_database()
    internal_database()