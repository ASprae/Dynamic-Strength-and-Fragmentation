import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

fname = 'Analysis'

samples = ['SaLi_002', 'SaLi_008', 'SaLi_014', 'CaMa_004', 'CaMa_009', 'CaMa_016', 'TaQu_006', 'TaQu_008', 'TaQu_019', 'SeeSst_014', 'SeeSst_022', 'SeeSst_028']

pixpercm = [113.39, 112.23, 156.91, 206.25, 165.85, 156.77, 115.66, 115.11, 137.95, 135.10, 126.11, 164.15]

material = ['SaLi', 'SaLi', 'SaLi', 'CaMa', 'CaMa', 'CaMa', 'TaQu', 'TaQu', 'TaQu', 'SeeSst', 'SeeSst', 'SeeSst']
erate = [58.7, 141.8, 345.5, 91.9, 126.7, 68.5, 183.9, 98.5, 40.1, 145.9, 48.1, 280.4]
erate_sig = [3.5, 24.4, 83.2, 20.1, 19.8, 15.1, 5.8, 10.2, 0.1, 8.8, 5.1, 42.7]

color_dict = { 'SaLi':'turquoise', 'CaMa':'teal', 'TaQu':'gold', 'SeeSst':'goldenrod' }

unique_M = set(material)
M_number = len(unique_M)

n = -1

C_data = []

Mod_C = []
Mean_C = []
SD_C = []
Skew_C = []
r_C = []

AR_data = []

Mod_AR = []
Mean_AR = []
SD_AR = []
Skew_AR = []
r_AR = []

val = float(input("Enter Minimum Sieve Size (cm): "))

for s in samples:

    print ("================= {} =================".format(s))

    if not os.path.exists('{}/{}'.format(s, fname)):
        os.makedirs('{}/{}'.format(s, fname))

    n = n+1

    res = pixpercm[n]

    # 1 - Area
    # 2 - X
    # 3 - Y
    # 4 - Perim
    # 5 - Major
    # 6 - Minor
    # 7 - Angle
    # 8 - Circularity
    # 9 - Max Feret Diameter
    # 10 - FeretX
    # 11 - Feret Y
    # 12 - Feret Angle
    # 13 - Min Feret Diameter
    # 14 - Aspect Ratio
    # 15 - Roundness
    # 16 - Solidity

    area, perim, maj, min, C, AR = np.genfromtxt('{}/{}.txt'.format(s, s), skip_header = 1, usecols = (1, 4, 5, 6, 8, 14), unpack = True)
    
    d_equi = 2 * np.sqrt(area/np.pi)
    Ax_ratio = 1 / AR

    # PLOT UNFILTERED SIZE DATA
    binnum = 50
    
    d_equivalent = d_equi[d_equi < 2 * val]
    
    fig = plt.figure(figsize = (6, 5))
    ax = fig.add_subplot(111)
    
    ax.set_ylabel('Frequency Density')
    ax.set_xlabel('Diameter (cm)')
    
    kde_min = stats.gaussian_kde(d_equivalent)
    span_min = np.linspace(0, np.amax(d_equivalent) * 1.1, num = 1000)
    kdepdf_min = kde_min.evaluate(span_min)

    ax.hist(d_equivalent, density = True, bins = binnum, color = 'darkblue', alpha = 0.25)
    ax.plot(span_min, kdepdf_min, c = 'darkblue')

    dips = np.where((kdepdf_min[1:-1] < kdepdf_min[0:-2]) * (kdepdf_min[1:-1] < kdepdf_min[2:]))[0] + 1
    
    if (len(dips) == 0):
        print ("ERROR - No minima in equivalent diameter")
    elif (len(dips) == 1):
        print ("Minimum diameter = {:.3f} cm".format(span_min[dips[-1]]))
        minD = span_min[dips[-1]]
        ax.axvline(x = minD, c = 'k')
    else:
        print ("More than 1 minima in equivalent diameter")
        print ("Picking the smallest diameter...")
        print ("Minimum diameter = {:.3f} cm".format(span_min[dips[0]]))
        minD = span_min[dips[0]]
        ax.axvline(x = minD, c = 'k')
        for d in dips:
            ax.axvline(x = span_min[d], c = 'k', ls = ':')
    
    samplesize = len(d_equi[d_equi > minD])

    print ("Eliminating {} Fragments...".format(len(d_equi[d_equi < minD])))
    print ("Sample Size = {} Fragments".format(samplesize))
    
    fig.savefig('{}/{}/F1_{}_sizefiltering.png'.format(s, fname, s))
    plt.close(fig)
    
    
    area = area[d_equi > minD]
    pix = area * res**2

    C = C[d_equi > minD]
    Ax_ratio = Ax_ratio[d_equi > minD]
    d_equi = d_equi[d_equi > minD]
    
    C_data.append(C)
    AR_data.append(Ax_ratio)
    
    print ("Smallest Fragment Area = {:.3f} cm^2 ({:.0f} pixels)".format(np.amin(area), np.amin(pix)))
    print ("--------------------------------------------")
    
    
    # Diameter Histogram
    fig = plt.figure(figsize = (6, 5))
    ax = fig.add_subplot(111)
    
    ax.set_ylabel('Frequency Density')
    ax.set_xlabel('Diameter (cm)')
    
    kde_min = stats.gaussian_kde(d_equi)
    span_min = np.linspace(0, np.amax(d_equi) * 1.1, num = 1000)
    kdepdf_min = kde_min.evaluate(span_min)

    ax.set_xlim([minD, np.amax(d_equi) * 1.1])

    ax.hist(d_equi, density = True, bins = binnum, color = 'darkblue', alpha = 0.25)
    ax.plot(span_min, kdepdf_min, c = 'darkblue')
    
    fig.savefig('{}/{}/F1_{}_size.png'.format(s, fname, s))
    plt.close(fig)
    
    
    # Shape Parameters Histogram
    fig = plt.figure(figsize = (12, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    for ax in [ax1, ax2]:
        ax.set_ylabel('Frequency Density')
        
    ax1.set_xlabel('Circularity')
    ax2.set_xlabel('Axial Ratio')
    
    ax1.set_xlim([0, 1])
    ax2.set_xlim([0, 1])
    
    # ax1 - Circularity
    hist, bin_edges = np.histogram(C, density = True, bins = binnum)
    kde = stats.gaussian_kde(C)
    
    span = np.linspace(0, 1, num = 1000)
    kdepdf = kde.evaluate(span)

    ax1.hist(C, density = True, bins = binnum, color = 'darkblue', alpha = 0.25)
    ax1.plot(span, kdepdf, c = 'darkblue')
    
    mode_C = span[np.argmax(kdepdf)]
    mean_C = np.mean(C)
    sd_C = np.std(C)
    
    Mod_C.append(mode_C)
    Mean_C.append(mean_C)
    SD_C.append(sd_C)
    Skew_C.append((mean_C - mode_C)/sd_C)
    
    print ("Modal Circularity = {:.3f}".format(mode_C))
    print ("Mean Circularity = {:.3f} +/- {:.3f}".format(mean_C, sd_C))
    print ("Circularity Skewness = {:.3f}".format((mean_C - mode_C)/sd_C))
    print ("--------------------------------------------")
    
    # ax2 - Axial ratio
    hist, bin_edges = np.histogram(Ax_ratio, density = True, bins = binnum)
    kde = stats.gaussian_kde(Ax_ratio)
    
    span = np.linspace(0, 1, num = 1000)
    kdepdf = kde.evaluate(span)

    ax2.hist(Ax_ratio, density = True, bins = binnum, color = 'darkblue', alpha = 0.25)
    ax2.plot(span, kdepdf, c = 'darkblue')
    
    mode_AR = span[np.argmax(kdepdf)]
    mean_AR = np.mean(Ax_ratio)
    sd_AR = np.std(Ax_ratio)
    
    Mod_AR.append(mode_AR)
    Mean_AR.append(mean_AR)
    SD_AR.append(sd_AR)
    Skew_AR.append((mean_AR - mode_AR)/sd_AR)
    
    print ("Modal Axial Ratio = {:.3f}".format(mode_AR))
    print ("Mean Axial Ratio = {:.3f} +/- {:.3f}".format(mean_AR, sd_AR))
    print ("Axial Ratio Skewness = {:.3f}".format((mean_AR - mode_AR)/sd_AR))
    print ("--------------------------------------------")
    
    
    fig.tight_layout()
    fig.savefig('{}/{}/{}_shape.png'.format(s, fname, s))
    
    plt.close(fig)


    # Shape vs Size
    fig = plt.figure(figsize = (12, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_ylabel('Circularity')
    ax2.set_ylabel('Axial Ratio')

    for ax in [ax1, ax2]:
        ax.set_xlabel('Diameter (cm)')
        
        ax.set_xlim([minD, np.amax(d_equi) * 1.1])
        ax.set_ylim([0, 1])
        
        
    # ax1 - Circularity
    k = stats.gaussian_kde(np.vstack([d_equi, C]), weights = area)
    xi, yi = np.mgrid[minD: d_equi.max()*1.1 :100j, 0 : 1 : 100j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    ax1.contourf(xi, yi, zi.reshape(xi.shape), alpha=0.5, cmap= 'Blues')
    ax1.scatter(d_equi, C, c = 'k', s = 0.1)
    
    r = stats.pearsonr(d_equi, C)
    print ("C-r = {:.3f}".format(r[0]))
    r_C.append(r[0])
    
    # ax2 - Axial ratio
    k = stats.gaussian_kde(np.vstack([d_equi, Ax_ratio]), weights = area)
    xi, yi = np.mgrid[minD: d_equi.max()*1.1 :100j, 0 : 1 : 100j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    ax2.contourf(xi, yi, zi.reshape(xi.shape), alpha=0.5, cmap= 'Blues')
    ax2.scatter(d_equi, Ax_ratio, c = 'k', s = 0.1)
    
    r = stats.pearsonr(d_equi, Ax_ratio)
    print ("AR-r = {:.3f}".format(r[0]))
    r_AR.append(r[0])
    
    fig.tight_layout()
    fig.savefig('{}/{}/{}_shape_area.png'.format(s, fname, s))
    plt.close(fig)


    # Shape Parameter Correlation
    fig = plt.figure(figsize = (6, 5))
    ax = fig.add_subplot(111)
    
    ax.set_ylabel('Circularity')
    ax.set_xlabel('Axial Ratio')
    
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

    p = ax.scatter(Ax_ratio, C, s = 10, c = d_equi, vmin = minD, vmax = np.amax(d_equi))
    
    b = fig.colorbar(p,orientation='vertical')
    b.set_label('Diameter (cm)')
    
    fig.savefig('{}/{}/F1_{}_shape_corr.png'.format(s, fname, s))
    plt.close(fig)


Sfig01 = plt.figure(figsize = (16, 10))
ax1 = Sfig01.add_subplot(2,3,1)
ax2 = Sfig01.add_subplot(2,3,4)

ax1b = Sfig01.add_subplot(2,3,2)
ax2b = Sfig01.add_subplot(2,3,5)

ax1c = Sfig01.add_subplot(2,3,3)
ax2c = Sfig01.add_subplot(2,3,6)

for ax in [ax1, ax1b, ax1c, ax2, ax2b, ax2c]:
    ax.set_xlabel('Strain Rate (/s)')
    ax.set_xlim([0, 450])
    
ax1.set_ylabel('Circularity')
ax2.set_ylabel('Axial Ratio')

ax1b.set_ylabel('Circularity Skewness')
ax2b.set_ylabel('Axial Ratio Skewness')

ax1c.set_ylabel('Circularity r')
ax2c.set_ylabel('Axial Ratio r')

ax1.set_ylim([0, 1])
ax2.set_ylim([0, 1])

ax1b.set_ylim([-1, 1])
ax2b.set_ylim([-1, 1])

ax1c.set_ylim([-1, 1])
ax2c.set_ylim([-1, 1])

for i in unique_M:
    filter = (np.array(material) == i)
    
    rate = np.array(erate)[filter]
    rate_sig = np.array(erate_sig)[filter]
    
    Cmode = np.array(Mod_C)[filter]
    Cmean = np.array(Mean_C)[filter]
    Csd = np.array(SD_C)[filter]
    Cskew = np.array(Skew_C)[filter]
    Cr = np.array(r_C)[filter]

    ARmode = np.array(Mod_AR)[filter]
    ARmean = np.array(Mean_AR)[filter]
    ARsd = np.array(SD_AR)[filter]
    ARskew = np.array(Skew_AR)[filter]
    ARr = np.array(r_AR)[filter]
    
    ax1.errorbar(rate, Cmean, xerr = rate_sig, yerr = Csd, fmt = 's', c = color_dict[i], label = i)
    ax2.errorbar(rate, ARmean, xerr = rate_sig, yerr = ARsd, fmt = 's', c = color_dict[i], label = i)
    
    ax1b.errorbar(rate, Cskew, xerr = rate_sig, fmt = 's', c = color_dict[i], label = i)
    ax2b.errorbar(rate, ARskew, xerr = rate_sig, fmt = 's', c = color_dict[i], label = i)
    
    ax1c.errorbar(rate, Cr, xerr = rate_sig, fmt = 's', c = color_dict[i], label = i)
    ax2c.errorbar(rate, ARr, xerr = rate_sig, fmt = 's', c = color_dict[i], label = i)
    
ax1.legend(loc = 4)

Sfig01.tight_layout()
Sfig01.savefig('Summary.png')



Sfig02 = plt.figure(figsize = (10, 8))
ax1 = Sfig02.add_subplot(2,1,1)
ax2 = Sfig02.add_subplot(2,1,2)

for ax in [ax1, ax2]:
    ax.set_xlabel('Strain Rate (/s)')
    ax.set_xlim([0, 450])
    
ax1.set_ylabel('Circularity')
ax2.set_ylabel('Axial Ratio')

ax1.set_ylim([0, 1])
ax2.set_ylim([0, 1])

vp1 = ax1.violinplot(C_data, erate, showmeans = False, showextrema = False, widths = 20)
vp2 = ax2.violinplot(AR_data, erate, showmeans = False, showextrema = False, widths = 20)

n = -1
for pc in vp1['bodies']:
    n = n+1
    if n < 3:
        pc.set_facecolor(color_dict['SaLi'])
        pc.set_edgecolor('k')
    elif n < 6:
        pc.set_facecolor(color_dict['CaMa'])
        pc.set_edgecolor('k')
    elif n < 9:
        pc.set_facecolor(color_dict['TaQu'])
        pc.set_edgecolor('k')
    elif n < 12:
        pc.set_facecolor(color_dict['SeeSst'])
        pc.set_edgecolor('k')
        
n = -1
for pc in vp2['bodies']:
    n = n+1
    if n < 3:
        pc.set_facecolor(color_dict['SaLi'])
        pc.set_edgecolor('k')
    elif n < 6:
        pc.set_facecolor(color_dict['CaMa'])
        pc.set_edgecolor('k')
    elif n < 9:
        pc.set_facecolor(color_dict['TaQu'])
        pc.set_edgecolor('k')
    elif n < 12:
        pc.set_facecolor(color_dict['SeeSst'])
        pc.set_edgecolor('k')

for i in unique_M:
    filter = (np.array(material) == i)

    rate = np.array(erate)[filter]
    rate_sig = np.array(erate_sig)[filter]
    
    Cmean = np.array(Mean_C)[filter]
    ARmean = np.array(Mean_AR)[filter]

    ax1.errorbar(rate, Cmean, xerr = rate_sig, fmt = 'o', c = color_dict[i], label = i, mec = 'k')
    ax2.errorbar(rate, ARmean, xerr = rate_sig, fmt = 'o', c = color_dict[i], label = i, mec = 'k')


ax1.legend()

Sfig02.tight_layout()
Sfig02.savefig('Violin_Summary.png')
