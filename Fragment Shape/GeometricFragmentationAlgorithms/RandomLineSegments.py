import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as stats
import bisect

from skimage import io
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.color import label2rgb

def on_segment(p, q, r):
    if r[0] <= max(p[0], q[0]) and r[0] >= min(p[0], q[0]) and r[1] <= max(p[1], q[1]) and r[1] >= min(p[1], q[1]):
        return True
    return False

def orientation(p, q, r):
    val = ((q[1] - p[1]) * (r[0] - q[0])) - ((q[0] - p[0]) * (r[1] - q[1]))
    if val == 0 : return 0
    return 1 if val > 0 else -1
    
def intersects(seg1, seg2):
    p1, q1 = seg1
    p2, q2 = seg2

    o1 = orientation(p1, q1, p2)

    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 and o3 != o4:
        return True
        
    if o1 == 0 and on_segment(p1, q1, p2) : return True
    if o2 == 0 and on_segment(p1, q1, q2) : return True
    if o3 == 0 and on_segment(p2, q2, p1) : return True
    if o4 == 0 and on_segment(p2, q2, q1) : return True

    return False

def intersection_coords(L1, L2): # Lx is a line segment defined as [(x1, y1), (x2, y2)]
            
    X1 = L1[0][0]
    Y1 = L1[0][1]
                
    X2 = L1[1][0]
    Y2 = L1[1][1]
                
    X3 = L2[0][0]
    Y3 = L2[0][1]
                
    X4 = L2[1][0]
    Y4 = L2[1][1]
    
    A1 = (Y1 - Y2)
    B1 = (X2 - X1)
    C1 = (X1*Y2 - X2*Y1) * -1
    
    A2 = (Y3 - Y4)
    B2 = (X4 - X3)
    C2 = (X3*Y4 - X4*Y3) * -1
    
    D = A1 * B2 - B1 * A2
    Dx = C1 * B2 - B1 * C2
    Dy = A1 * C2 - C1 * A2
    if D !=0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False



seeds = np.linspace(100, 1000, num = 10)

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
    
    #bounding lines
    coord_1 = (0,0)
    coord_2 = (0,1)
    coord_3 = (1,1)
    coord_4 = (1,0)
    
    line_1 = [coord_1, coord_2]
    line_2 = [coord_2, coord_3]
    line_3 = [coord_3, coord_4]
    line_4 = [coord_4, coord_1]
    
    L = [line_1, line_2, line_3, line_4]
    LINES = []
    
    for i in np.arange(0, s, 1):
    
        x1 = points1[i, 0]
        y1 = points1[i, 1]
        
        A = (x1, y1)
    
        x2 = points2[i, 0]
        y2 = points2[i, 1]
        
        B = (x2, y2)
        
        x_LINE = np.linspace(0, 1, num = 2)
        y_LINE = (((y2-y1)/(x2-x1)) * (x_LINE - x1)) + y1
        
        A = (x_LINE[0], y_LINE[0])
        B = (x_LINE[1], y_LINE[1])
        
        L_new = [A, B]
        
        i_sections = []
        
        for j in np.arange(0, len(L), 1):
        
            if intersects(L_new, L[j]) is True:
                IS = intersection_coords(L_new, L[j])
                
                i_sections.append(IS)
                    
        i_sections = sorted(i_sections, key=lambda x: x[0])
        
        x_chosen = np.random.uniform(i_sections[0][0], i_sections[-1][0])
        
        x_values = []
        
        for j in np.arange(0, len(i_sections), 1):
            x_values.append(i_sections[j][0])
        
        ind = bisect.bisect_left(x_values, x_chosen) - 1
        
        Line = [(i_sections[ind][0], i_sections[ind][1]), (i_sections[ind + 1][0], i_sections[ind + 1][1])]
        
        L.append(Line)
        LINES.append(Line)
        
    for ls in LINES:
        x = [ls[0][0], ls[1][0]]
        y = [ls[0][1], ls[1][1]]
    
        ax.plot(x, y)
    
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])

    ax.axis('off')

    fig.tight_layout()
    fig.savefig('RandLineSegs_{:05d}.png'.format(s))
    plt.close(fig)
    
    
    
    
    image = io.imread('RandLineSegs_{:05d}.png'.format(s), as_gray=True)

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

vp1 = ax4.violinplot(C_data, seeds, showmeans = False, showmedians = False, showextrema = False, widths = 80)
vp2 = ax5.violinplot(AR_data, seeds, showmeans = False, showmedians = False, showextrema = False, widths = 80)

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

