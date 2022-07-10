"""
Author: Carolina Rojas Ramirez
Date: 11/09/2021
"""
from tkinter import filedialog
import pickle
from terminalFragmentor_Main import FragmentSite
from terminalFragmentor_Main import ThyIon

ions_file = filedialog.askopenfilename(title='Choose Ions File', filetypes=[('Ions File', '.ions')])
with open(ions_file, 'rb') as picklefile:
    saved_database = pickle.load(picklefile)
    protein_seq = saved_database[0]
    analysis_dict = saved_database[1]

    # print(analysis_dict)

    #Analysis
    for analysis in analysis_dict:
        print(analysis)
        # print(analysis_dict[analysis])
        for analysis_pass in analysis_dict[analysis]:
            print(analysis_pass)
            # print(analysis_dict[analysis][analysis_pass])
            for site in analysis_dict[analysis][analysis_pass][1]:
                print(site)
                fragsite = analysis_dict[analysis][analysis_pass][1][site]
                # print(fragsite)
                # print(fragsite.theo_ions)
                for ion in fragsite.theo_ions:

                    print(f"{fragsite.theo_ions[ion].mz_mono}_{fragsite.theo_ions[ion].charge}_{fragsite.theo_ions[ion].ion_type}{fragsite.theo_ions[ion].ion_type_indx}_{fragsite.theo_ions[ion].thy_mods}")