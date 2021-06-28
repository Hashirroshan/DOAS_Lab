
import glob
import numpy as np

import pandas as pd
import os


from scipy.interpolate import interp1d

from scipy import linalg
from scipy import stats


import matplotlib.pyplot as plt

import numpy as np

# Data measurements with cell and without cell
no_cell_wave = glob.glob('without_cell_dat/NC*.dat')
sorted_no_cw=sorted(no_cell_wave)
cell_wave = glob.glob('with_cell_dat/WC*.dat')
sorted_cw=sorted(cell_wave)
dark = np.loadtxt('dark/dark_1.dat', comments = '*')

len(cell_wave)

Wavelength = []
I0 = []
I_c = []
no2_conc = np.zeros(len(cell_wave))
time = np.zeros(len(cell_wave))

# substract dark current from I0 and I

for i in range(len(cell_wave)):
    I_cell = np.loadtxt(sorted_cw[i], comments = '*')
    I_without_dark = I_cell - dark
    I_c = pd.DataFrame(I_without_dark)

    I0_no_cell = np.loadtxt(sorted_no_cw[i], comments='*')
    I0_without_dark = I0_no_cell - dark
    I0 = pd.DataFrame(I0_without_dark)


    # 3. Wavelength Calibration


    """
    lambdai is the wavelength
    N is the number of CCD pixels: for NO2: N = 1024
    a0 .. a2 are the calibration coefficients which have the values

    """
N = 1024
lambda_i = np.zeros(N)
for i in range(1024):
    a0 = 429.494
    a1 = 93.112
    a2 = -6.050
    lambda_i[i] = a0 + a1 * ((i - 1) / (N - 1)) + a2 * ((i - 1) / (N - 1)) * ((i - 1) / (N - 1))
    # lambda_i = a0+a1*((i-1)/(N-1))+a2*((i-1)/(N-1))*((i-1)/(N-1))
    Wavelength = pd.DataFrame(lambda_i)

# Task:4. choose the values of I and I0 in the fitting window: NO2: 432.5 – 465 nm

I0.insert(0, "wavelength", Wavelength[0], True)
I0.insert(0, "I", I_c[1], True)
I0 = I0[I0.wavelength <= 465]
I0 = I0[I0.wavelength >= 432.5]
I0.reset_index(drop=True, inplace=True)

# Task 5: calculate ln(I0/I)

calc_log = (I0[1] / I0['I'])
opt_depth = pd.DataFrame(calc_log)
depth_value = np.log(opt_depth)
ln_I0_I = pd.DataFrame(depth_value)

#Task 6: fit a polynomial of order 3 to ln(I0/I)

x = I0['wavelength']
fit_a_polynomial=np.polyfit(x, ln_I0_I, 3)

#Task 7: subtract the fitted polynomial from ln(I0/I) to get the differential ln(I0/I)

fitted_polynomial = np.polyval(fit_a_polynomial, x)
fitted_polynomial_df = pd.DataFrame(fitted_polynomial)
diff = ln_I0_I - fitted_polynomial_df

#Task 8: interpolate the differential cross section to the measurement wavelengths

crosssection = np.loadtxt(os.path.join("help_files/NO2_DIFFXSECTION.DAT"), comments='*')
cr = pd.DataFrame(crosssection)

cr = cr[cr[0] <= 465]
cr = cr[cr[0] >= 432.5]
cr.reset_index(drop=True, inplace=True)
X = np.array(cr[0])
Y = np.array(cr[1])
f = interp1d(X, Y)

#Task 9: DOAS fit, fit the differential cross-section to the differential ln(I0/I) :
"""
    Note: use column vectors for the cross-sections and the differential OD. The function
    also returns the standard error of the fitted parameter.
"""

cr.insert(2, "differential", diff[0], True)
D_C = f(x)  # take differential absorption cross-section
Diff_Cross = np.array(D_C)
Diff_lnI_I = np.array(diff[0])

#Task 10: Calculate the concentration of NO2 in the cell

m, c, r_value, p_value, std_err = stats.linregress(Diff_Cross, Diff_lnI_I)
y = m * Diff_Cross + c
no2_slant_column = m
no2concentration = m * 1e+18

no2_conc=abs(no2concentration/10) #dividing by cell height

#Task11: Plot results
plt.figure(1)
plt.plot(I0['wavelength'], I0[1])
plt.plot(I0['wavelength'], I0['I'])
plt.title('Intensity with cell and Intensity without cell for Last data set')
plt.legend(['Intensity without NO\u2082 cell', 'Intensity with NO\u2082 cell'])
plt.ylabel('Intensity without dark current ')
plt.xlabel('Wavelength [nm]')
plt.show()
plt.savefig('Task_1_2.png')


plt.figure(2)
plt.plot(I0['wavelength'], ln_I0_I[0])
plt.plot(I0['wavelength'], diff[0])
plt.title('Optical Density and Differential Optical Density for Last data set')
plt.legend(['ln (I0/I) vs wavelength', 'differential optical density vs wavelength'])
plt.ylabel('Optical density')
plt.xlabel('Wavelength [nm]')
plt.show()
plt.savefig('Task_1_3.png')

plt.figure(3)
plt.plot(X, f(X), '-')
plt.scatter(X, Y)
plt.title('Differential absorption Cross-section of NO\u2082 Interpolation',pad = 20)
plt.ylabel('Differential absorption Cross-section of NO\u2082 [cm\u00b2] ')
plt.xlabel('Wavelength[nm]')
plt.show()
plt.savefig('Task_1_4.png')

plt.figure(4)
plt.scatter(Diff_Cross, Diff_lnI_I)
plt.plot(Diff_Cross, y, 'r')
plt.title('NO\u2082 Slant Column')
plt.ylabel('Differential Optical Depth ')
plt.xlabel('NO\u2082 Differential Cross-section [cm\u00b2/molec] 10^-18')
plt.show()
plt.savefig('Task_1_5.png')
