from pathlib import Path
import glob
import numpy as np

import pandas as pd
import os

from scipy import stats
from scipy.interpolate import interp1d

from scipy import linalg
from scipy import stats
from scipy.optimize import curve_fit

from uncertainties import ufloat


# path = Path('with_cell_1')
# for file in path.glob('*.dat'):
#     with open(file) as f:
#         with open(str(file)+'1.csv', "w") as f1:
#             for Counts in f:
#                 f1.write(Counts)


import matplotlib.pyplot as plt

import numpy as np

# Data measurements with cell (wavelengths and Counts)
no_cell_wave = glob.glob('without_cell_dat/NC*.dat')
sorted_no_cw=sorted(no_cell_wave)
cell_wave = glob.glob('with_cell_dat/WC*.dat')
sorted_cw=sorted(cell_wave)
dark = np.loadtxt('dark/dark_1.dat', comments = '*')

len(cell_wave)

# cell_dir = os.path.dirname("with_cell_dat/WC__*.DAT")
# no_cell_dir = os.path.dirname("without_cell_dat/NC__*.dat")
# dark_dir = os.path.dirname("dark/")
# crs_sec_dir = os.path.dirname("help_files/")
#
# cell = glob.glob(os.path.join(cell_dir, "WC__*.DAT"))
# no_cell = glob.glob(os.path.join(no_cell_dir, "NC__*.dat"))
# dark = np.loadtxt(os.path.join(dark_dir, "dark_1.dat"), comments='*')

Wavelength = []
I0 = []
I_c = []
no2_conc = np.zeros(len(cell_wave))
time = np.zeros(len(cell_wave))

value_1=[]
value_2=[]
value_3=[]
value_4=[]
value_6=[]
value_7=[]
value_8=[]
value_5=[]
# value_12=[]
# value_13=[]
# value_14=[]

# substract dark current from I0 and I

for a in range(len(cell_wave)):
    I_cell = np.loadtxt(sorted_cw[a], comments = '*')
    I_without_dark = I_cell - dark
    I_c = pd.DataFrame(I_without_dark)

    I0_no_cell = np.loadtxt(sorted_no_cw[a], comments='*')
    I0_without_dark = I0_no_cell - dark
    I0 = pd.DataFrame(I0_without_dark)


# Wavelength = []
# I0 = []
# I_x = []
# no2con = np.zeros(157)
# time = np.zeros(157)
#
# # 2. subtract dark current from I0 and I :
#
# for a in range(157):
#     I_dt = np.loadtxt(cell[a], comments='*')
#     I_withoutdark = I_dt - dark
#     I_x = pd.DataFrame(I_withoutdark)
#
#     I0_dt = np.loadtxt(no_cell[a], comments='*')
#     I0_withoutdark = I0_dt - dark
#     I0 = pd.DataFrame(I0_withoutdark)

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
    Diff_Cross_x = np.array(D_C)
    Diff_lnI_I_y = np.array(diff[0])

#Task 10: Calculate the concentration of NO2 in the cell

    m, c, r_value, p_value, std_err = stats.linregress(Diff_Cross_x, Diff_lnI_I_y)
    print(std_err)
    y = m * Diff_Cross_x + c
    # no2_slant_column = m
    # no2concentration = m
    #
    # no2_conc=abs(no2concentration/10) #dividing by cell height
    no2_slant_column = m
    no2_conc=abs(no2_slant_column/10) #dividing by cell height
    no2concentration = no2_conc

    error=std_err
    no2_conc=no2concentration
    value_1.append(no2_slant_column)
    value_2.append(error)
    value_3.append(no2_conc)
    value_4.append(a)

Er_Rep = ufloat(no2_slant_column, std_err) #sc value = sc +/- sc_uncertainty
Hg1_cell = ufloat(10,0.01) #cell height = 10 +/- 0.1
Hg2_cell = ufloat(10,0)  #cell height = 10 +/- 0
value_er_1 = Er_Rep/Hg1_cell   # computation of VC assuming 1mm uncertainty in cell height: VC value = VC +/- VC_uncertainty (including fit error and cell uncertainty)
er_1=value_er_1.std_dev
value_5.append(value_er_1)
value_6.append(er_1)
value_er_2 = Er_Rep/Hg2_cell  # computation of VC assuming no uncertainty in cell height: VC value = VC +/- VC_uncertainty (including just fit error)
er_2=value_er_2.std_dev
value_7.append(value_er_2)
value_8.append(er_2)


NO2_SC_std_err_concentration=pd.DataFrame(value_1)
NO2_SC_std_err_concentration.columns = ['NO2_SC_molec/cm^2']
value_2_df=pd.DataFrame(value_2)
value_3_df=pd.DataFrame(value_3)
value_4_df=pd.DataFrame(value_4)
value_5_df=pd.DataFrame(value_5)
value_6_df=pd.DataFrame(value_6)
value_7_df=pd.DataFrame(value_7)
value_8_df=pd.DataFrame(value_8)
value_6_df.columns = ['N02con_fitted _&_1mm_NO2cell_uncertainty']
value_6_df.insert(1,"N02con_just_fitted _error", value_8_df[0], True)


NO2_SC_std_err_concentration.insert(1, "NO2_Con_molec/cm^3", value_3_df[0], True)
NO2_SC_std_err_concentration.insert(2, "Data_set_NO", value_4_df[0], True)
NO2_SC_std_err_concentration.insert(3, "NO2_std_err", value_2_df[0], True)
pd.set_option("display.max_rows", None, "display.max_columns", None)
result = pd.DataFrame(NO2_SC_std_err_concentration)
j = result.to_csv('error_analysis.csv')
# print(value_6_df)
print(result)


#Task11: Plot results

plt.figure(1)
plt.scatter(time, no2_conc)
plt.title('NO2 Concentration in Cell vs Time')
plt.ylim(ymin=0)
plt.ylabel('NO\u2082 Concentration in Cell [molec/cm\u00b3]')
plt.xlabel('TIME (from 8:34 am to 10:26 am)')
plt.show()
plt.savefig('Task_1_1.png')

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
plt.scatter(Diff_Cross_x, Diff_lnI_I_y)
plt.plot(Diff_Cross_x, y, 'r')
plt.title('NO\u2082 Slant Column')
plt.ylabel('Differential Optical Depth ')
plt.xlabel('NO\u2082 Differential Cross-section [cm\u00b2/molec] 10^-18')
plt.show()
plt.savefig('Task_1_5.png')
