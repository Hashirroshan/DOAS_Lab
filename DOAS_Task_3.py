import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.interpolate import interp1d




slant_columns_DATA = np.loadtxt("slant_columns/NO2_SC_BREMEN_20200501.DAT", comments="*")
SC_DATA = pd.DataFrame(slant_columns_DATA)
# SC_DATA.columns = ['Time', 'Solar Zenith Angle', 'NO2 Slant Column', 'NO2 Slant Column uncertainty %']
AMF_DATA = np.loadtxt("help_files/NO2_AMF.DAT", comments='*')q
AMF = pd.DataFrame(AMF_DATA)
AMF.columns = ['angle', 'Air_Mass_Factor']

# SC_DATA['SZAAM'] = SC_DATA['NO2 Slant Column'] / AMF['Air_Mass_Factor']
#
# plt.plot(SC_DATA['Time'],SC_DATA['SZAAM'])
# plt.ylabel('NO2 VERTICAL COLUMN [10e15 molec cm-2]')
# plt.xlabel('TIME [UT]')
# plt.show()


g = interp1d(AMF.angle, AMF.Air_Mass_Factor)
SC_DATA["AMF"] = g(SC_DATA[1])
AMF_at_38 g(38)
SC_DATA["Delta_AMF"] = SC_DATA["AMF"] - AMF_at_38
SC_DATA["VC"] = SC_DATA[2] / SC_DATA["Delta_AMF"]


# calculation @ sunrise SR and sunset SS

SR = pd.DataFrame(SC_DATA)
SS = pd.DataFrame(SC_DATA)

SR = SR[SR[0] <= 8]
SS = SS[SS[0] >= 15.1]

temp = np.array(SR[1])
I_VC = np.array(SR['VC'])
i = interp1d(temp, I_VC)
VC_90_Sunrise = i(90)

temp2 = np.array(SS[1])
J_VC = np.array(SS['VC'])
j = interp1d(temp2, J_VC)
VC_90_Sunset = j(90)

plt.figure(1)
plt.plot(SC_DATA[0], SC_DATA['VC'])
plt.title(' NO\u2082 vertical columns above Bremen',pad = 20)
plt.ylabel('NO\u2082 Vertical Columns [molec/cm\u00b2]')
plt.xlabel('time [UT]')
plt.show()
plt.savefig('Task_3_1.png')

plt.figure(2)
plt.errorbar(SC_DATA[0], SC_DATA['VC'],ecolor='red', yerr = np.multiply(SC_DATA[3], SC_DATA['VC'])/100 ,fmt='o')      #yerr = 1e16*0.5*SC_DATA[3]
plt.plot(SC_DATA[0], SC_DATA['VC'])
plt.title(' NO\u2082 vertical columns above Bremen with Uncertainities',pad = 20)
plt.ylabel('NO\u2082 Vertical Columns [molec/cm\u00b2]')
plt.xlabel('time [UT]')
plt.savefig('Task_3_2.png')

print('vertical column during the sunrise when sun at 90° is %f' % VC_90_Sunrise)
plt.figure(3)
plt.errorbar(SR[0], SR['VC'],ecolor='red',yerr = np.multiply(SR[3], SR['VC'])/100 ,fmt='o')
plt.plot(SR[0], SR['VC'])
plt.title(' NO\u2082 VC above Bremen during Sunrise', pad = 20)
plt.ylabel('NO\u2082 Vertical Columns [molec/cm\u00b2]')
plt.xlabel('time [UT]')
plt.show()
plt.savefig('Task_3_3.png')

print('vertical column during the sunset when sun at 90° is %f' % VC_90_Sunset)
plt.figure(4)
plt.errorbar(SS[0], SS['VC'],ecolor='red',yerr = np.multiply(SS[3], SS['VC'])/100 ,fmt='o')
plt.plot(SS[0], SS['VC'])
plt.title(' NO\u2082 VC above Bremen during Sunset',pad = 20)
plt.ylabel('NO\u2082 Vertical Columns [molec/cm\u00b2]')
plt.xlabel('time [UT]')
plt.show()
plt.savefig('Task_3_4.png')
