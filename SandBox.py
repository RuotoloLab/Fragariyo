"""Testing my own peak picking script
https://medium.com/@chrisjpulliam/quickly-finding-peaks-in-mass-spectrometry-data-using-scipy-fcf3999c5057
"""
from tkinter import filedialog
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks
import ms_deisotope
import time


coordinatesfile = filedialog.askopenfilename(title='Select Coordinate files', filetypes=[('CSV', '.csv')])

ms_data = pandas.read_csv(coordinatesfile, header=1)
ms_data.columns = ["mz", "int"]
print(ms_data)
# default kind is line
ms_data.plot(x='mz', y='int', color='blue')
plt.title('MS spectrum')
plt.xlabel('m/z')
plt.ylabel('Intensity')


peak_idx, _= find_peaks(ms_data["int"],
                          prominence = 75,
                          height = 150,
                          distance = None)

peak_data = ms_data.iloc[peak_idx]

#
# print(peak_data)
# print(peak_data.index)
print(peak_data.values.tolist())


sns.scatterplot(x = peak_data["mz"], y = peak_data["int"],
                color = 'red', alpha = 0.5)
plt.show()
# It seems I don't need to peak picked...it will do it for you...never mind is taking froever! (more than 12 hrs)
# need to prepare first into PeakSet Object: https://github.com/mobiusklein/ms_deisotope?tab=readme-ov-file
start_time = time.time()
preparedpeaks = ms_deisotope.deconvolution.utils.prepare_peaklist(peak_data.values.tolist())
print(preparedpeaks)

print("--- %s seconds preparing ---" % (time.time() - start_time))

start_time = time.time()
deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(preparedpeaks, averagine=ms_deisotope.peptide,
                                                       scorer=ms_deisotope.MSDeconVFitter(10.), use_quick_charge=True)
print(deconvoluted_peaks)

print("--- %s seconds deconvoluting ---" % (time.time() - start_time))

for peak in deconvoluted_peaks:
    print(peak)


