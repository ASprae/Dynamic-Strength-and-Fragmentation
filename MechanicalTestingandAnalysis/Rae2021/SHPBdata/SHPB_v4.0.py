import os as os
import csv
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import integrate

#######################
## INPUT INFORMATION ##
#######################
folder = 'SeeSst'
sample = 'SeeSst_SHPB_001_Ti25_30b_Al12b_0_001'

sample_length = 42.00 # mm
sample_diameter = 41.00 # mm

striker_length = 0.25 # m
striker_material = 'Ti' # Al, Ti, or Steel
striker_velocity = np.nan # Unnecessary!

t_sampling_f = 1.25E+6 #1.25E+6

#Pulse Width Factor - Factors used to trim data correctly, values should be edited to cut the data over the correct interval
ffactor = 6.0 # 6.0 - 7.0 usually works for the pulse-shapers in this study
bfactor = 3.5 # 3.0 - 4.5 usually works for the pulse-shapers in this study
#######################

############################
## Non-Variable Constants ##
############################
sample_diameter = sample_diameter/1000
sample_length = sample_length/1000

if striker_material == 'Al':
    striker_E = 70000000000 # Pa
    striker_rho = 2810 # kg/m^3
elif striker_material == 'Ti':
    striker_E = 110000000000 # Pa
    striker_rho = 4430 # kg/m^3
elif striker_material == 'Steel':
    striker_E = 190000000000 # Pa
    striker_rho = 8000 # kg/m^3
else:
    print ('ERROR - FAILED TO SELECT VALID STRIKER MATERIAL')

C0 = np.sqrt(striker_E/striker_rho) # m/s

striker_diameter = 50.0 # mm
incident_length = 1.25000 #1.2525 #m
#incident_length = 1.50000 #1.2525 #m
transmitted_length = 1.0000 #1.0025 #m

t_sampling_t = 1/t_sampling_f

striker_diameter = striker_diameter/1000

U0 = 5 # V
V = 400 # Amplification factor
k = 1.99 # Strain gauge k-factor

if not os.path.exists(folder + '/' + sample):
    os.mkdir(folder + '/' + sample)

Dlamda_int, CnC0_int = np.genfromtxt('Pochhammer-Chree.csv', skip_header=1, delimiter = ',', unpack=True)
############################

print ("-------------------")
print (sample)
print ("-------------------")

# Load Data
t, Inc_1, Inc_2, Tra_1, Tra_2 = np.genfromtxt('{}/{}.txt'.format(folder,sample),skip_header=11, delimiter = ';', unpack=True, dtype=None)

for data in [t, Inc_1, Inc_2, Tra_1, Tra_2]:
    data = data.tolist()

time = []
I1 = []
I2 = []
T1 = []
T2 = []

# Users may need to uncomment the L87-91, and comment out L92-96 depending on the nature of the data file (Seeberger Sandstone inputs have commas as periods but the remaining inputs just have periods!)

for i in np.arange(0, len(t),1):
    #time.append((float(t[i])))
    #I1.append((float(Inc_1[i])))
    #I2.append((float(Inc_2[i])))
    #T1.append((float(Tra_1[i])))
    #T2.append((float(Tra_2[i])))
    time.append((float(t[i].replace(',', '.'))))
    I1.append((float(Inc_1[i].replace(',', '.'))))
    I2.append((float(Inc_2[i].replace(',', '.'))))
    T1.append((float(Tra_1[i].replace(',', '.'))))
    T2.append((float(Tra_2[i].replace(',', '.'))))

I1 = np.array(I1)
I2 = np.array(I2)
T1 = np.array(T1)
T2 = np.array(T2)
time = np.array(time)


# Plot Raw Data
fig01a = plt.figure(figsize=(8,4))
ax01 = fig01a.add_subplot(111)

ax01.set_xlabel('Time (s)')
ax01.set_ylabel('Strain Gauge Voltage (V)')

ax01.plot(time, I1, label = "I1")
ax01.plot(time, I2, label = "I2")
ax01.plot(time, T1, label = "T1")
ax01.plot(time, T2, label = "T2")

ax01.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax01.legend()

fig01a.tight_layout()
fig01a.savefig('{}/{}/F1a_RawData.png'.format(folder,sample),dpi=300)


########## FILTER DATA - USER CAN MAKE CHANGES HERE IF REQUIRED ##########
print ("-------------------")
print ("Did any of the strain-gauges fall off?")
choice = float(input("Yes (1) or No (0)? - "))

if choice == 1:
    #I1 = ma.masked_where(time > 0, I1)
    #I2 = ma.masked_where(time > 0.00, I2)
    #T1 = ma.masked_where(time > 0.0, T1)
    #T2 = ma.masked_where(time > 0.000, T2)

    # Plot Raw Data
    fig01b = plt.figure(figsize=(6,3))
    ax01 = fig01b.add_subplot(111)

    ax01.set_xlabel('Time (s)')
    ax01.set_ylabel('Strain Gauge Voltage (V)')

    ax01.plot(time, I1, label = "I1")
    ax01.plot(time, I2, label = "I2")
    ax01.plot(time, T1, label = "T1")
    ax01.plot(time, T2, label = "T2")

    ax01.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    ax01.legend()

    fig01b.tight_layout()
    fig01b.savefig('{}/{}/F1b_FilteredData.png'.format(folder,sample),dpi=300)

    I = ma.array((I1, I2)).mean(axis=0)
    T = ma.array((T1, T2)).mean(axis=0)

elif choice == 0:
    I_data = np.array([I1, I2])
    T_data = np.array([T1, T2])
    
    I = np.ma.mean(I_data, axis = 0)
    T = np.average(T_data, axis = 0)
    print ("Great!")
else:
    print ("*************************************")
    print ("ERROR - FAILED TO SELECT VALID OPTION")
    print ("*************************************")

# Convert Voltages to Strains
Ei = (4 * I) / (U0 * V * k)
Et = (4 * T) / (U0 * V * k)

#Calculate Time corrections
inc_Tcorr = incident_length/C0
ref_Tcorr = -incident_length/C0
tra_Tcorr = -transmitted_length/C0

time_incident = time + inc_Tcorr

time_inc = time_incident - np.amin(time_incident)
time_ref = time + ref_Tcorr - np.amin(time_incident)
time_tra = time + tra_Tcorr - np.amin(time_incident)

t_pulse = 2 * (striker_length/C0)
t_pulse_indexwidth = int(t_pulse/t_sampling_t)

t_peak = time_inc[(Ei.tolist()).index(np.amax(Ei))]

index_i = np.argmax(Ei) #Determines the time-index of the peak strain in the incident bar
index_r = (np.abs(time_ref - t_peak)).argmin() #Determines time-index of equivalent to "index_i" in time-adjusted incident bar signal
index_t = (np.abs(time_tra - t_peak)).argmin() #Determines time-index of equivalent to "index_i" in time-adjusted transmitted bar signal

inc_arg_start = int(index_i - (ffactor/2 * t_pulse_indexwidth)) #Index of Start of Incident Wave signal
inc_arg_end = int(index_i + (bfactor/2 * t_pulse_indexwidth)) + 1 #Index of End of Incident Wave signal

ref_arg_start = int(index_r - (ffactor/2 * t_pulse_indexwidth)) #Index of Start of Reflected Wave signal
ref_arg_end = int(index_r + (bfactor/2 * t_pulse_indexwidth)) + 1 #Index of End of Reflected Wave signal

tra_arg_start = int(index_t - (ffactor/2 * t_pulse_indexwidth)) #Index of Start of Transmitted Wave signal
tra_arg_end = int(index_t + (bfactor/2 * t_pulse_indexwidth)) + 1 #Index of End of Transmitted Wave signal

if inc_arg_start < 0:
    print ("*******************************************************")
    print ("ERROR - WINDOW IS TOO LARGE, CALLING DATA OUTSIDE FRAME")
    print ("********** Try reducing the value of 'factor' *********")
    print ("*******************************************************")

########## FINAL DATA ARRAYS ##########
TIME = np.array(time_inc.tolist()[inc_arg_start: inc_arg_end] - np.amin(time_inc.tolist()[inc_arg_start: inc_arg_end]))
STRAIN_I = np.array(Ei.tolist()[inc_arg_start: inc_arg_end])
STRAIN_R = np.array(Ei.tolist()[ref_arg_start: ref_arg_end])
STRAIN_T = np.array(Et.tolist()[tra_arg_start: tra_arg_end])
#######################################

### Pochhammer-Chree Filtering ###
tmax = np.amax(TIME)
omega = (2 * np.pi) / tmax
N = 30
n = np.arange(1,N+1,1)
iter_n = 20

areas_I = []
for i in np.arange(0,len(STRAIN_I)-1,1):
    areas_I.append(((STRAIN_I[i] + STRAIN_I[i+1])/2) * t_sampling_t)
areas_I = np.array(areas_I)
A0_I = np.sum(areas_I) * (2/tmax)

areas_R = []
for i in np.arange(0,len(STRAIN_R)-1,1):
    areas_R.append(((STRAIN_R[i] + STRAIN_R[i+1])/2) * t_sampling_t)
areas_R = np.array(areas_R)
A0_R = np.sum(areas_R) * (2/tmax)

areas_T = []
for i in np.arange(0,len(STRAIN_T)-1,1):
    areas_T.append(((STRAIN_T[i] + STRAIN_T[i+1])/2) * t_sampling_t)
areas_T = np.array(areas_T)
A0_T = np.sum(areas_T) * (2/tmax)

a_n_I = []
b_n_I = []
a_n_R = []
b_n_R = []
a_n_T = []
b_n_T = []
for m in n:
    a_I = []
    b_I = []
    a_R = []
    b_R = []
    a_T = []
    b_T = []
    for i in np.arange(0,len(STRAIN_I)-1,1):
        a_I.append(areas_I[i] * np.cos(m * omega * (TIME[i] + TIME[i+1])/2))
        b_I.append(areas_I[i] * np.sin(m * omega * (TIME[i] + TIME[i+1])/2))
        a_R.append(areas_R[i] * np.cos(m * omega * (TIME[i] + TIME[i+1])/2))
        b_R.append(areas_R[i] * np.sin(m * omega * (TIME[i] + TIME[i+1])/2))
        a_T.append(areas_T[i] * np.cos(m * omega * (TIME[i] + TIME[i+1])/2))
        b_T.append(areas_T[i] * np.sin(m * omega * (TIME[i] + TIME[i+1])/2))
    a_n_I.append(np.sum(a_I) * (2/tmax))
    b_n_I.append(np.sum(b_I) * (2/tmax))
    a_n_R.append(np.sum(a_R) * (2/tmax))
    b_n_R.append(np.sum(b_R) * (2/tmax))
    a_n_T.append(np.sum(a_T) * (2/tmax))
    b_n_T.append(np.sum(b_T) * (2/tmax))

a_n_I = np.array(a_n_I)
b_n_I = np.array(b_n_I)
a_n_R = np.array(a_n_R)
b_n_R = np.array(b_n_R)
a_n_T = np.array(a_n_T)
b_n_T = np.array(b_n_T)

D_n_I = np.sqrt(a_n_I**2 + b_n_I**2)
D_n_R = np.sqrt(a_n_R**2 + b_n_R**2)
D_n_T = np.sqrt(a_n_T**2 + b_n_T**2)

###################################################
fig02a = plt.figure(figsize=(6,3))
ax01 = fig02a.add_subplot(111)

ax01.plot(n,D_n_I)
ax01.plot(n,D_n_R)
ax01.plot(n,D_n_T)

ax01.set_xlabel('Wavenumber')
ax01.set_ylabel('Strain')

fig02a.tight_layout()
fig02a.savefig('{}/{}/F2a_FourierAnalysis.png'.format(folder, sample),dpi=300)
###################################################

phi_n_I = []
phi_n_R = []
phi_n_T = []
for m in np.arange(0,N,1):
    if a_n_I[m] < 0:
        phi_n_I.append(np.arctan(b_n_I[m]/a_n_I[m]) + np.pi)
    else:
        phi_n_I.append(np.arctan(b_n_I[m]/a_n_I[m]))

    if a_n_R[m] < 0:
        phi_n_R.append(np.arctan(b_n_R[m]/a_n_R[m]) + np.pi)
    else:
        phi_n_R.append(np.arctan(b_n_R[m]/a_n_R[m]))

    if a_n_T[m] < 0:
        phi_n_T.append(np.arctan(b_n_T[m]/a_n_T[m]) + np.pi)
    else:
        phi_n_T.append(np.arctan(b_n_T[m]/a_n_T[m]))

lamda_i = t_pulse * C0 / n

###################################################
fig02b = plt.figure(figsize=(6,3))
ax01 = fig02b.add_subplot(111)

for k in np.arange(0,iter_n,1):
    ci = (n * omega) / (2 * np.pi) * lamda_i
    x = striker_diameter / lamda_i
    
    CnC0_int_new = np.interp(x, Dlamda_int, CnC0_int)
    Cj = CnC0_int_new * C0
    
    lamda_j = Cj * 2 * np.pi/(n * omega)
    
    CnC0 = Cj / C0
    Dla = striker_diameter / lamda_j
    
    ax01.plot(Dla, CnC0)
    
    lamda_i = lamda_j

ax01.plot(Dlamda_int, CnC0_int, ls = '--', lw = 3, c='k', label = 'P-C')
ax01.plot(Dla, CnC0, ls = ':', lw = 1, c='0.5', label = 'Result')

ax01.set_xlabel('D / $\lambda_n$')
ax01.set_ylabel('$C_n$ \ $C_0$')

ax01.legend()

fig02b.tight_layout()
fig02b.savefig('{}/{}/F2b_PC-Interp.png'.format(folder,sample),dpi=300)
###################################################

Cn = Cj

phi_dn_I = (n * omega) * ((incident_length/Cn) - (incident_length/C0))
phi_dn_R = (n * omega) * ((-incident_length/Cn) - (-incident_length/C0))
phi_dn_T = (n * omega) * ((-transmitted_length/Cn) - (-transmitted_length/C0))

phi_new_I = phi_n_I + phi_dn_I
phi_new_R = phi_n_R + phi_dn_R
phi_new_T = phi_n_T + phi_dn_T

sc_I = []
sc_R = []
sc_T = []
for i in TIME:
    scv_I = D_n_I * np.cos((n * omega * i) - phi_new_I)
    sc_I.append((A0_I/2) + np.sum(scv_I))

    scv_R = D_n_R * np.cos((n * omega * i) - phi_new_R)
    sc_R.append((A0_R/2) + np.sum(scv_R))

    scv_T = D_n_T * np.cos((n * omega * i) - phi_new_T)
    sc_T.append((A0_T/2) + np.sum(scv_T))

###################################################
fig02c = plt.figure(figsize=(6,3))
ax01 = fig02c.add_subplot(111)

ax01.plot(TIME, STRAIN_I, c = 'r', ls = '--', lw = 0.5)
ax01.plot(TIME, STRAIN_R, c = 'g', ls = '--', lw = 0.5)
ax01.plot(TIME, STRAIN_T, c = 'b', ls = '--', lw = 0.5)

STRAIN_I = np.array(sc_I)
STRAIN_R = np.array(sc_R)
STRAIN_T = np.array(sc_T)

ax01.plot(TIME, STRAIN_I, c = 'r', label='$\epsilon_{I}$')
ax01.plot(TIME, STRAIN_R, c = 'g', label='$\epsilon_{R}$')
ax01.plot(TIME, STRAIN_T, c = 'b', label='$\epsilon_{T}$')
ax01.plot(TIME, STRAIN_T-STRAIN_R, c = 'k', label='$\epsilon_{T} - \epsilon_{R}$')

ax01.legend()

ax01.set_xlabel('Time (s)')
ax01.set_ylabel('Strain')

ax01.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

fig02c.tight_layout()
fig02c.savefig('{}/{}/F2c_PCcorrectedI.png'.format(folder,sample),dpi=300)
###################################################

# Plot Time-Adjusted Strain Data
fig03a = plt.figure(figsize=(6,3))
ax01 = fig03a.add_subplot(111)

ax01.set_xlabel('Time (s)')
ax01.set_ylabel('Strain')

ax01.plot(TIME, STRAIN_I, c = 'r', label='$\epsilon_{I}$')
ax01.plot(TIME, STRAIN_R, c = 'g', label='$\epsilon_{R}$')
ax01.plot(TIME, STRAIN_T, c = 'b', label='$\epsilon_{T}$')
ax01.plot(TIME, STRAIN_T-STRAIN_R, c = 'k', lw=0.5, ls = '--', label='$\epsilon_{T} - \epsilon_{R}$')

ax01.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax01.legend()

fig03a.tight_layout()
fig03a.savefig('{}/{}/F3a_IRTwaves.png'.format(folder,sample),dpi=300)



# Calculate Stress, Strain Rate, and Strain
sigma_A = ((striker_diameter**2)/(sample_diameter**2)) * striker_E * (STRAIN_I + STRAIN_R)
sigma_B = ((striker_diameter**2)/(sample_diameter**2)) * striker_E * (STRAIN_T)

sigma_av = ((striker_diameter**2)/(sample_diameter**2)) * striker_E * ((STRAIN_I + STRAIN_R + STRAIN_T)/2)

strain_rate = (C0/sample_length) * (STRAIN_I - STRAIN_R - STRAIN_T)

strain = (C0/sample_length) * integrate.cumtrapz((STRAIN_I - STRAIN_R - STRAIN_T), dx = t_sampling_t, initial=0)

# Calculate Energies
A_bar = np.pi * (striker_diameter/2)**2

Const =  A_bar * C0 * striker_E

E_I_cum = Const * integrate.cumtrapz(STRAIN_I**2, x = TIME)
E_R_cum = Const * integrate.cumtrapz(STRAIN_R**2, x = TIME)
E_T_cum = Const * integrate.cumtrapz(STRAIN_T**2, x = TIME)

E_del_cum = E_I_cum - E_R_cum
E_abs_cum = E_I_cum - E_R_cum - E_T_cum

E_del = E_del_cum[-1]
E_abs = E_abs_cum[-1]

# Plot Stress and Strain-rate Data
fig04a = plt.figure(figsize=(14,8))
ax01 = fig04a.add_subplot(221)
ax02 = fig04a.add_subplot(222)
ax03 = fig04a.add_subplot(223)
ax04 = fig04a.add_subplot(224)

ax01.set_xlabel('Time (s)')
ax01.set_ylabel('Stress (Pa)')

ax01.fill_between(TIME, sigma_A, sigma_B, alpha = 0.25, color ='k')
ax01.plot(TIME, sigma_A, c='k', ls = ':', lw = 0.5, label = '$\sigma_A$')
ax01.plot(TIME, sigma_B, c='k', ls = '--', lw = 0.5, label = '$\sigma_B$')
ax01.plot(TIME, sigma_av, c='r', ls = '-', lw = 1, label = '$\sigma_{mean}$')
ax01.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax01.legend()


ax02.set_xlabel('Time (s)')
ax02.set_ylabel('Strain Rate ($s^{-1}$)')

ax02.plot(TIME, strain_rate, c = 'r')
ax02.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax03.set_xlabel('Time')
ax03.set_ylabel('Strain')

ax03.plot(TIME, strain, c = 'r')
ax03.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax03.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax04.set_xlabel('Strain')
ax04.set_ylabel('Stress (Pa)')

ax04.fill_between(strain, sigma_A, sigma_B, alpha = 0.25, color ='k')
ax04.plot(strain, sigma_A, c='k', ls = ':', lw = 0.5, label = '$\sigma_A$')
ax04.plot(strain, sigma_B, c='k', ls = '--', lw = 0.5, label = '$\sigma_B$')
ax04.plot(strain, sigma_av, c='r', ls = '-', lw = 1, label = '$\sigma_{mean}$')

ax04.ticklabel_format(style='sci', axis='x', scilimits=(0,0))



# Print Results
Y = np.amax(sigma_av) / 1000000
a = np.abs(sigma_A[np.argmax(sigma_av)] - np.amax(sigma_av)) / 1000000
b = np.abs(sigma_B[np.argmax(sigma_av)] - np.amax(sigma_av)) / 1000000
Y_err = np.mean([a, b])

e_max = strain[np.argmax(sigma_av)]

e_r = np.ma.masked_where(sigma_av < (0.5 * Y * 1000000), strain_rate)
e_r = np.ma.masked_where(TIME > TIME[np.argmax(sigma_av)], e_r)
e_rate = np.ma.mean(e_r)
e_rate_err = np.ma.std(e_r)

a = np.ma.masked_where(TIME > TIME[np.argmax(sigma_av)], sigma_av)
b = np.ma.masked_where(sigma_av < (0.25 * Y * 1000000), a)
c = np.ma.masked_where(sigma_av > (0.75 * Y * 1000000), b)
mask = np.ma.getmask(c)
masked_e = np.ma.masked_array(strain, mask = mask)
masked_sigma = np.ma.masked_array(sigma_av, mask = mask)
masked_e = masked_e.compressed()
masked_sigma = masked_sigma.compressed()
a, V = np.polyfit(masked_e, masked_sigma, 1, cov = True)
E = a[0]
E_err = np.sqrt(V[0][0])

print ("----------- RESULTS -------------")
print ("Y = {:.1f} +/- {:.1f} MPa".format(Y, Y_err))
print ("e_max = {:.4f}".format(e_max))
print ("e_r = {:.1f} +/- {:.1f} /s".format(e_rate, e_rate_err))
print ("E = {:.1f} +/- {:.1f} GPa".format(E * 1.E-9, E_err * 1.E-9))
print ("E_abs/del = {:.1f} J / {:.1f} J".format(E_abs, E_del))
print ("---------------------------------")

ax04.text(0, np.amax(sigma_av)*0.975, 'Y = {:.1f} +/- {:.1f} MPa'.format(Y, Y_err), horizontalalignment='left', verticalalignment='center')
ax04.text(0, np.amax(sigma_av)*0.9, r'$\epsilon_m$ = {:.4f}'.format(e_max), horizontalalignment='left', verticalalignment='center')
ax04.text(0, np.amax(sigma_av)*0.825, r'$\dot \epsilon$ = {:.1f} +/- {:.1f} /s'.format(e_rate, e_rate_err), horizontalalignment='left', verticalalignment='center')
ax04.text(0, np.amax(sigma_av)*0.75, 'E = {:.1f} +/- {:.1f} GPa'.format(E * 1.E-9, E_err * 1.E-9), horizontalalignment='left', verticalalignment='center')

ax04.legend(loc = 6)

fig04a.tight_layout()
fig04a.savefig('{}/{}/F4a_DeformationCharacterisation.png'.format(folder, sample),dpi=300)


# Plot Energy Information
fig04b = plt.figure(figsize=(15,5))
ax01 = fig04b.add_subplot(131, aspect = 'equal')
ax02 = fig04b.add_subplot(132)
ax03 = fig04b.add_subplot(133)

x = np.linspace(0, np.amax(E_del_cum), num = 10)
y = x

ax01.plot(E_del_cum, E_T_cum, c = 'r')
ax01.plot(x, y, c = 'k', ls=':', label = 'Full Transmission')

ax01.legend()

strain_NEW = np.delete(strain, -1)
sigma_av_NEW = np.delete(sigma_av, -1)

ax02.plot(E_abs_cum, strain_NEW , c = 'r')
ax03.plot(E_abs_cum, sigma_av_NEW , c = 'r')

ax01.set_xlabel('Energy Delivered (J)')
ax01.set_ylabel('Energy Transmitted (Pa)')

ax02.set_xlabel('Energy Absorbed (J)')
ax02.set_ylabel('Strain')

ax03.set_xlabel('Energy Absorbed (J)')
ax03.set_ylabel('Stress (Pa)')

ax02.text(np.amax(E_abs_cum)*0.975, np.amax(strain)*0.00, 'Total Energy Delivered = {:.1f} J'.format(E_del), horizontalalignment='right', verticalalignment='center')
ax02.text(np.amax(E_abs_cum)*0.975, np.amax(strain)*0.05, 'Final Energy Absorbed = {:.1f} J'.format(E_abs), horizontalalignment='right', verticalalignment='center')
ax02.text(np.amax(E_abs_cum)*0.975, np.amax(strain)*0.10, r'$\epsilon_m$ = {:.4f}'.format(e_max), horizontalalignment='right', verticalalignment='center')
ax02.text(np.amax(E_abs_cum)*0.975, np.amax(strain)*0.15, r'$\dot \epsilon$ = {:.1f} +/- {:.1f} /s'.format(e_rate, e_rate_err), horizontalalignment='right', verticalalignment='center')

fig04b.tight_layout()
fig04b.savefig('{}/{}/F4b_EnergyCharacterisation.png'.format(folder, sample),dpi=300)



DATA = zip(TIME, sigma_av, sigma_A, sigma_B, strain, strain_rate)
np.savetxt('{}/{}/output.txt'.format(folder, sample), np.array(DATA), delimiter = ',', header = "Time,Stress,Stress_A,Stress_B,Strain,Strain_Rate")
