from sys import displayhook
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import os
import pyscal.traj_process as ptp
import numpy as np
from matplotlib import pyplot as plt
import glob
import seaborn as sns
from scipy.integrate import simps
import pandas as pd
#from IPython.display import display

#Read File
sys = pc.System()
atom = sys.atoms
timestep = 200000000#provide the timestep for the single dump file

sys.read_inputfile('Coord1.25.02/ParticleCoordW.'+str(timestep), format='lammps-dump', compressed=False, customkeys=['radius', 'c_PerParticleStress[4]']) # HP: doesnt work - customkeys FK: Does work now!

# q 4, and 6 values, make sure the cut off distance 
#Make sure cut-off is twice of the big particle radius
sys.find_neighbors(method="cutoff", cutoff = 2.5) # filter = type works apply a filter to nearest neighbor calculation. If the filter keyword is set to type, only atoms of the same type would be included in the neighbor calculations. If type_r, only atoms of a different type will be included in the calculation. Default None.

sys.calculate_q([4,6])
q4 = sys.get_qvals(4)
q6 = sys.get_qvals(6)

q4_ave = np.mean(q4)
q4_std = np.std(q4)
q6_ave = np.mean(q6)
q6_std = np.mean(q6)
print ("q4 = ", q4_ave, "q6 = ", q6_ave)
bins = np.linspace(0, 1, 100)

q4_hist, bins = np.histogram(q4, bins=bins)
q6_hist, bins = np.histogram(q6, bins=bins)


plt.figure(figsize=(3, 2), dpi = 1200)
plt.hist(q4, bins=bins, alpha=0.5, label='q4', color='blue')
plt.hist(q6, bins=bins, alpha=0.5, label='q6', color='green')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of q4, q6, and q8 Values')
plt.legend()
plt.savefig("Histogram_q4q6q8.svg", format='svg')

# q_average 4, 6, and 8 values 
sys.find_neighbors(method="cutoff", cutoff = 3.4) # Make sure cut-off distance is set contact distance

sys.calculate_q([4,6], averaged = True)
q4_avg = sys.get_qvals([4], averaged = True)
q6_avg = sys.get_qvals([6], averaged = True)

q4_avg_hist, bins = np.histogram(q4_avg, bins=bins)
q6_avg_hist, bins = np.histogram(q6_avg, bins=bins)

plt.figure(figsize=(4, 3), dpi = 1200)
plt.hist(q4_avg, bins=bins, alpha=0.5, label='q4_avg', color='blue')
plt.hist(q6_avg, bins=bins, alpha=0.5, label='q6_avg', color='green')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of Average q4, q6, and q8 Values')
plt.legend()
plt.savefig("HistogramAve_q4q6q8.svg", format='svg')

# plot probabilty density for q_i
plt.figure(figsize=(4, 3), dpi = 1200)  # Fix the width and height here ...
plt.rcParams["font.family"] = "Arial"
plt.tick_params(axis='both', which='both', top=True, right=True)
plt.gca().spines['top'].set_linewidth(1.5)  # Top edge thickness
plt.gca().spines['bottom'].set_linewidth(1.5)  # Bottom edge thickness
plt.gca().spines['left'].set_linewidth(1.5)  # Left edge thickness
plt.gca().spines['right'].set_linewidth(1.5)  # Right edge thickness
plt.xlim([0, 1])
plt.xticks(np.linspace(0,1,5))
plt.tick_params(axis='both', which='minor',labelsize=10)
sns.kdeplot(q4, color='blue', label='q4')
sns.kdeplot(q6, color='green', label='q6')
plt.xlabel(r'$q_l$',fontsize=14)
plt.ylabel('P($q_l$)',fontsize=14)
plt.tick_params(axis='both', labelsize=10,which='both', direction='in', length=5, width=1)
plt.legend([r'$q_4$', r'$q_6$'])
plt.tight_layout()
plt.savefig("Probability_Density_P_q_layered.svg", format='svg')
#plt.show()
