import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d
import glob
import matplotlib.pyplot as plt
import os

# 1. Load data files :

hrz_dir = os.path.dirname("horizon_measurements/")

cell_H = glob.glob("horizontal/horiz_*.DAT")
I0_dt_H = np.loadtxt('horizontal/zenith.DAT', comments='*')
dark_H = np.loadtxt('horizontal/dark.DAT', comments='*')
# print(len(cell_H))
no2con_H = np.zeros(20)
time_H = np.zeros(20)

# 2. subtract dark current from I0 and I :

I0_without_dark_H = I0_dt_H - dark_H
I0_H = pd.DataFrame(I0_without_dark_H)
Wavelength_H = []
I_H = []
I0_H_fix = I0_H.copy()


for a in range(15):
    I0_H = I0_H_fix.copy()
    I_dt_H = np.loadtxt(cell_H[a], comments='*')
    I_without_dark_H = I_dt_H - dark_H
    I_H = pd.DataFrame(I_without_dark_H)

    # 3. Wavelength Calibration :

    wl_H = np.zeros(1024)
    for j in range(1024):
        a0 = 429.494
        a1 = 93.112
        a2 = -6.050
        N = 1024

        wl_H[j] = a0 + a1 * ((j - 1) / (N - 1)) + a2 * ((j - 1) / (N - 1)) * ((j - 1) / (N - 1))
        Wavelength_H = pd.DataFrame(wl_H)

    # 4. choose the values of I and I0 in the fitting window: NO2: 432.5 – 465 nm

    I0_H.insert(0, "wavelength", Wavelength_H[0], True)
    I0_H.insert(0, "I", I_H[1], True)
    I0_H = I0_H[I0_H.wavelength <= 465]
    I0_H = I0_H[I0_H.wavelength >= 432.5]
    I0_H.reset_index(drop=True, inplace=True)

    # 5. calculate ln(I0/I)  :

    OD_H = I0_H[1] / I0_H['I']
    depth_value_H = np.log(OD_H)
    ln_I0_I_H = pd.DataFrame(depth_value_H)

    # 6. fit a polynomial of order 3 to ln(I0/I) :

    x_H = I0_H['wavelength']
    fit_a_polynomial_H = np.polyfit(x_H, ln_I0_I_H, 3)

    # 7. subtract the fitted polynomial from ln(I0/I) to get the differential ln(I0/I) :

    fitted_polynomial_H = np.polyval(fit_a_polynomial_H, x_H)
    fitted_polynomial_df_H = pd.DataFrame(fitted_polynomial_H)
    diff_H = ln_I0_I_H - fitted_polynomial_df_H

    # 8. interpolate the differential cross section to the measurement wavelengths :

    crosssection_H = np.loadtxt(os.path.abspath("help_files/NO2_DIFFXSECTION.DAT"), comments='*')
    R = pd.DataFrame(crosssection_H)

    R = R[R[0] <= 465]
    R = R[R[0] >= 432.5]
    R.reset_index(drop=True, inplace=True)

    X_H = np.array(R[0])
    Y_H = np.array(R[1])
    f_H = interp1d(X_H, Y_H)

    # 9. DOAS fit, fit the differential cross-section to the differential ln(I0/I) :

    R.insert(2, "differential", diff_H[0], True)
    D_C_H = f_H(x_H)  # take differential absorption cross-section
    Diff_Cross_H = np.array(D_C_H * 1e+18)
    Diff_lnI_I_H = np.array(diff_H[0])

    # 10. Calculate the concentration of NO2 in the cell :

    m1, c1, r_value1, p_value1, std_err1 = stats.linregress(Diff_Cross_H, Diff_lnI_I_H)
    y_H = m1 * Diff_Cross_H + c1
    no2concentration_H = m1 * 1e+18

    no2con_H[a] = no2concentration_H
    time_H[a] = a


# 11. Plot results :

plt.figure(1)
plt.scatter(time_H,no2con_H)
plt.title('NO\u2082 Slant Column vs Angle ')
plt.ylim(ymin=0)
plt.ylabel('NO\u2082 Slant Column [molec/cm\u00b2] ')
plt.xlabel('Elevation Angle [°]')
plt.show()
plt.savefig('Task_2_1.png')


plt.figure(2)
plt.plot(I0_H['wavelength'], I0_H[1])
plt.plot(I0_H['wavelength'], I0_H['I'])
plt.title('Intensity at 20° elevation angle')
plt.legend(['Intensity at zenith angle', 'Intensity at 20° elevation angle'])
plt.ylabel('Intensity without dark current ')
plt.xlabel('Wavelength [nm]')
plt.show()
plt.savefig('Task_2_2.png')

plt.figure(3)
plt.plot(I0_H['wavelength'], ln_I0_I_H[0])
plt.plot(I0_H['wavelength'], diff_H[0])
plt.title('Optical Density and Differential Optical Density vs wavelength')
plt.legend(['ln (I0/I)_H vs wavelength', 'differential optical density vs wavelength'])
plt.ylabel('Optical density')
plt.xlabel('Wavelength [nm]')
plt.show()
plt.savefig('Task_2_3.png')

plt.figure(4)
plt.plot(X_H, f_H(X_H), '-')
plt.scatter(X_H, Y_H)
plt.title('Differential absorption Cross-section of NO\u2082 Interpolation',pad = 20)
plt.ylabel('Differential absorption Cross-section of NO\u2082 [cm\u00b2]')
plt.xlabel('Wavelength[nm]')
plt.show()
plt.savefig('Task_2_4.png')

plt.figure(5)
plt.scatter(Diff_Cross_H, Diff_lnI_I_H)
plt.plot(Diff_Cross_H, y_H, 'r')
plt.title('NO\u2082 Slant Column')
plt.ylabel('Differential Optical Depth ')
plt.xlabel('NO\u2082 Differential Cross-section [cm\u00b2/molec] 10^-18')
plt.show()
plt.savefig('Task_2_5.png')
