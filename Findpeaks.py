import scipy
import pandas as pd
from tkinter import filedialog
from scipy import signal
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import os


def rms_calc(array):
    """
    Calculates Root Mean Square of an array
    :param array: an array or array slice
    :return: RMS
    """
    y = array

    # Square all values in the array, take the arithmetic mean and square it
    rms = np.sqrt(np.mean(y ** 2))

    return rms

def noise_xtractor():
    """
    Function to parse a .csv file with the experimental neutral masses
    :return: a list of experimental masses
    """
    out_str = ""
    expfile = filedialog.askopenfilename(title='XY Coordinates for noise', filetypes=[('CSV', '.csv')])

    pathsplit = expfile.split('/')
    # proteiname = pathsplit[-1].strip("_expions.csv")
    proteiname = pathsplit[-1]
    print(proteiname)

    pd_file = pd.read_csv(expfile, header=0, engine='python')
    pd_mz = pd_file['mz']
    pd_int = pd_file['int']

    # print(pd_file.head())
    # print(pd_file['int'])

    np_mz = pd_mz.to_numpy()
    np_int = pd_int.to_numpy()

    # print(np_int)
    # return  np_mz, np_int

    rms = rms_calc(np_int)
    print(f"RMS of noise = {rms}")

    return rms

def spectra_xtractor():
    """
    Function to parse a .csv file with the experimental neutral masses
    :return: a list of experimental masses
    """
    main_outdir = filedialog.askdirectory(title='Choose Output Folder')
    os.chdir(main_outdir)

    out_str = ""
    expfiles = filedialog.askopenfilenames(title='XY Coordinates', filetypes=[('CSV', '.csv')])
    print(expfiles)

    noise = noise_xtractor()
    for file in expfiles:
        pathsplit = file.split('/')
        # proteiname = pathsplit[-1].strip("_expions.csv")
        proteiname = pathsplit[-1]
        print(proteiname)

        pd_file = pd.read_csv(file, header=0, engine='python')
        pd_mz = pd_file['mz']
        pd_int = pd_file['int']

        # print(pd_file.head())
        # print(pd_file['int'])

        np_mz = pd_mz.to_numpy()
        np_int = pd_int.to_numpy()

        # print(np_int)
        # return  np_mz, np_int

        rms = rms_calc(np_int)
        print(f"RMS = {rms}")


        sn = rms/noise
        print(f"S/N = {sn}")

        # For maxima only
        # peaks = scipy.signal.find_peaks(np_int)
        # print(peaks)
        # print(len(peaks[0]))
        # prominances = scipy.signal.peak_prominences(np_int, peaks[0])
        # print(f"Prominances = {prominances[0]}")
        #
        # output_intstr = ""
        # output_mzstr = ""
        # for index in peaks[0]:
        #     output_intstr += f"{np_int[index]};"
        #     output_mzstr += f"{np_mz[index]};"

        #Non-maxima
        output_intstr = ""
        output_mzstr = ""
        for mz in np_mz:
            output_mzstr += f"{mz};"
        for int in np_int:
            output_intstr += f"{int};"

        output_mzstr = output_mzstr.strip(";")
        output_intstr = output_intstr.strip(";")

        out_str += f"{proteiname}\n{output_mzstr}\t{output_intstr}\t{sn}\n"

        print("\n"+output_intstr)
        print(output_mzstr+"\n")

    output = open("Batch_XY-extractor_z03" + '.tsv', 'w')
    output.write(out_str)
    output.close()

def isotope_xtractor(main_outdir):

    """
    Script to get the isotope envelopes for the experimental ions. Input to internal fragments

    """
    os.chdir(main_outdir)
    expions = filedialog.askopenfilenames(title='expions to find (e.g unmatched peaks)', filetypes=[('CSV', '.csv')])
    print(expions)

    for expfile in expions:

        pathsplit = expfile.split('/')
        samplename = pathsplit[-1].strip(".csv")
        print(samplename)

        #Initiate dictiornay to store the coordinates for the peaks
        expion_dict = {}

        pd_file = pd.read_csv(expfile, header=0, engine='python')
        # pd_mz = pd_file['mz']
        # pd_int = pd_file['z']


        for index in range(len(pd_file)):
            ion_info = pd_file.loc[index]
            # print(ion_info[0], ion_info[1], ion_info[2])
            expion_dict[ion_info[0]] = []
            expion_dict[ion_info[0]].append(ion_info[1])
            expion_dict[ion_info[0]].append(ion_info[2])

        # print(expion_dict)

        coordinates = filedialog.askopenfilenames(title='XY Coordinates', filetypes=[('CSV', '.csv'),('XY', '.xy')])

        for file in coordinates:
            pathsplit = file.split('/')
            # proteiname = pathsplit[-1].strip("_expions.csv")
            coorname = pathsplit[-1]
            print(f"coorname = {coorname}")
            coornamespl = coorname.split('.')
            coorfiletype = coornamespl[-1]

            # print(f"coorname type = {coorfiletype}")


            if coorfiletype == "csv":
                # Agilent
                # pd_file = pd.read_csv(file, header=0, engine='python', usecols=['X(MassToCharge)', 'Y(Counts)'])
                pd_file = pd.read_csv(file, header=None, engine='python')
                # print(pd_file.head())
                pd_file.columns = ['X(MassToCharge)', 'Y(Counts)']
                # print(pd_file.head())
            else:
                # Breuker
                pd_file = pd.read_csv(file, header=0, engine='python', usecols=[0,1], sep="\t")
                pd_file.columns = ['X(MassToCharge)', 'Y(Counts)']


            for peak in expion_dict:
                round_peak = round(peak,4)
                rndone_mzrange = pd_file.loc[pd_file['X(MassToCharge)'] > round_peak - 1 ]
                # print(rndone_mzrange)
                # print(type(rndone_mzrange))

                rndtwo_mzrange = rndone_mzrange.loc[rndone_mzrange['X(MassToCharge)'] < round_peak + 3]
                # print(rndtwo_mzrange)

                np_coormz = rndtwo_mzrange['X(MassToCharge)'].to_numpy()
                np_coorint = rndtwo_mzrange['Y(Counts)'].to_numpy()

                # print(np_coormz)
                # print(np_coorint)

                #";" is a separator so that .csv file can be used
                mz_output = ";".join(str(x) for x in np_coormz)
                int_output = ";".join(str(x) for x in np_coorint)

                expion_dict[peak].append((mz_output, int_output))

            print(expion_dict)

            final_df_ls = []
            for val in expion_dict:
                charge = expion_dict[val][0]
                intensity = expion_dict[val][1]
                neutral_val = (val*charge) - (charge*1.0078)
                mz_val = expion_dict[val][2][0]
                int_val = expion_dict[val][2][1]

                final_df_ls.append([neutral_val,charge,val,intensity,mz_val, int_val])
            print(final_df_ls)

            final_df = pd.DataFrame(final_df_ls,
                              columns=['#neutral_mono','z','mz','int', 'isoenv_mz','isoenv_int'])

            print(final_df)

            final_df.to_csv(f"{samplename}_expion_isoenv.csv", index=False)




            # # For maxima only
            # peaks = scipy.signal.find_peaks(np_int)
            # print(peaks)
            # print(len(peaks[0]))
            # prominances = scipy.signal.peak_prominences(np_int, peaks[0])
            # print(f"Prominances = {prominances[0]}")





if __name__ == "__main__":

    # spectra_xtractor()
    isotope_xtractor()


    #
    # # Define some test data which is close to Gaussian
    # data = np.random.normal(size=10000)
    #
    # hist, bin_edges = np.histogram(data, density=True)
    # bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    #
    #
    # # Define model function to be used to fit to the data above:
    # def gauss(x, *p):
    #     A, mu, sigma = p
    #     return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))
    #
    #
    # # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    # p0 = [1., 0., 1.]
    #
    # coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    #
    # # Get the fitted curve
    # hist_fit = gauss(bin_centres, *coeff)
    # print(hist_fit)
    #
    # plt.plot(bin_centres, hist, label='Test data')
    # plt.plot(bin_centres, hist_fit, label='Fitted data')
    #
    # # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
    # print('Fitted mean = ', coeff[1])
    # print('Fitted standard deviation = ', coeff[2])
    #
    # plt.show()