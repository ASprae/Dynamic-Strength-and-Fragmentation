import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as stats

from skimage import io
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb

seeds = np.linspace(25, 250, num = 10)

binnum = 50

fig1 = plt.figure(figsize = (6,5))
ax1 = fig1.add_subplot(111)

ax1.set_ylabel('Frequency Density')
ax1.set_xlabel('Diameter (pixels)')


fig2 = plt.figure(figsize = (12,5))
ax2 = fig2.add_subplot(121)
ax3 = fig2.add_subplot(122)

for ax in [ax2, ax3]:
    ax.set_ylabel('Frequency Density')
    
ax2.set_xlabel('Circularity')
ax3.set_xlabel('Axial Ratio')

ax2.set_xlim([0, 1])
ax3.set_xlim([0, 1])


fig3 = plt.figure(figsize = (6, 2.5))
ax4 = fig3.add_subplot(121)
ax5 = fig3.add_subplot(122)

for ax in [ax4, ax5]:
    ax.set_xlabel('Number of Lines')
    ax.set_xlim([0, np.amax(seeds)*1.1])
    
ax4.set_ylabel('Circularity')
ax5.set_ylabel('Axial Ratio')

ax4.set_ylim([0, 1])
ax5.set_ylim([0, 1])

C_data = []
AR_data = []


fig4 = plt.figure(figsize = (12,5))
ax6 = fig4.add_subplot(121)
ax7 = fig4.add_subplot(122)

for ax in [ax6, ax7]:
    ax.set_xlabel('Diameter (pixels)')
    
ax6.set_ylabel('Circularity')
ax7.set_ylabel('Axial Ratio')

ax6.set_ylim([0, 1])
ax7.set_ylim([0, 1])


for s in seeds:

    s = int(s)

    points1 = np.random.rand(s, 2)
    points2 = np.random.rand(s, 2)
    
    fig = plt.figure(figsize = (20,20))
    ax = fig.add_subplot(111, aspect = 'equal')
    
    for i in np.arange(0, s, 1):
    
        x1 = points1[i, 0]
        y1 = points1[i, 1]
    
        x2 = points2[i, 0]
        y2 = points2[i, 1]
    
        x = np.linspace(0, 1, num = 2)
        y = (((y2-y1)/(x2-x1)) * (x - x1)) + y1
    
        ax.plot(x, y)

    ax.set_xlim([0,1])
    ax.set_ylim([0,1])

    ax.axis('off')

    fig.tight_layout()
    fig.savefig('RandLines_{:05d}.png'.format(s))
    plt.close(fig)



    image = io.imread('RandLines_{:05d}.png'.format(s), as_gray=True)

    thresh = threshold_otsu(image)
    bw = image > thresh

    cleared = clear_border(bw)

    label_image = label(cleared)

    image_label_overlay = label2rgb(label_image, image=image, bg_label=0)

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(image_label_overlay)

    regions = regionprops(label_image)

    d_eq = []
    C = []
    AR = []

    for props in regions:
        width = props.minor_axis_length
        length = props.major_axis_length
    
        perimeter = props.perimeter
        area = props.area
        
        count =  np.count_nonzero([width, length, perimeter, area])
    
        if count == 4:
            d_equivalent = 2 * np.sqrt(area/np.pi)
            Circ = 4 * np.pi * (area / (perimeter**2))
            AxR = width/length
    
            d_eq.append(d_equivalent)
            C.append(Circ)
            AR.append(AxR)
    
    C_data.append(C)
    AR_data.append(AR)
    
    ax.set_title('{} Lines'.format(s))
    
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig('Binary_{:05d}.png'.format(s))


    kde_min = stats.gaussian_kde(d_eq)
    span_min = np.linspace(0, np.amax(d_eq) * 1.1, num = 1000)
    kdepdf_min = kde_min.evaluate(span_min)

    ax1.plot(span_min, kdepdf_min, label = s)
    
    kde = stats.gaussian_kde(C)
    span = np.linspace(0, 1, num = 1000)
    kdepdf = kde.evaluate(span)

    ax2.plot(span, kdepdf, label = s)


    kde = stats.gaussian_kde(AR)
    span = np.linspace(0, 1, num = 1000)
    kdepdf = kde.evaluate(span)

    ax3.plot(span, kdepdf, label = s)
    
    ax6.scatter(d_eq, C, label = s, s = 5)
    ax7.scatter(d_eq, AR, label = s, s = 5)


ax1.legend()
ax2.legend()
ax6.legend()

fig1.savefig('F1_size.png')
plt.close(fig)

fig2.savefig('F2_shape.png')
plt.close(fig)

fig4.savefig('F4_size_shape.png')
plt.close(fig)


C_med = []
for Circul in C_data:
    C_med.append(np.median(Circul))
    
AR_med = []
for AxialR in AR_data:
    AR_med.append(np.median(AxialR))
    
vp1 = ax4.violinplot(C_data, seeds, showmeans = False, showmedians = False, showextrema = False, widths = 16)
vp2 = ax5.violinplot(AR_data, seeds, showmeans = False, showmedians = False, showextrema = False, widths = 16)

for pc in vp1['bodies']:
    pc.set_facecolor('k')
    pc.set_alpha(0.25)
    pc.set_edgecolor('k')
    
for pc in vp2['bodies']:
    pc.set_facecolor('k')
    pc.set_alpha(0.25)
    pc.set_edgecolor('k')

ax4.scatter(seeds, C_med, edgecolors = 'k', c = 'w')
ax5.scatter(seeds, AR_med, edgecolors = 'k', c = 'w')

fig3.tight_layout()
fig3.savefig('F3_Violin_Summary.png')
