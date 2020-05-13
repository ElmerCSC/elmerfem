#
# Make a report from output files of the inverse method
# make basic plots of cost function and norm of the gradient 
#  and report last lines of M1QN3_<RUN_NAME>.out
# Require files: M1QN3_<RUN_NAME>.out
#                Cost_<RUN_NAME>.dat
#                CostReg_<RUN_NAME>.dat
#                GradientNormAdjoint_<RUN_NAME>.dat
from collections import deque
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import genfromtxt;
import numpy as np
import sys

try:
	sys.argv[1]
except:
	print('error: try ')
	print('	python MakeReport.py RUN_NAME')
	exit()

# Output from M1QN3
fname='M1QN3_%s.out'%(sys.argv[1])
try:
	with open(fname) as file:
		pass
except IOError as e:
	print('Unable to open file %s'%(fname))
	exit()

with open(fname) as fin:
    last = deque(fin, 7)
output=''.join(list(last))

# Cost file
fname='Cost_%s.dat'%(sys.argv[1])
try:
	with open(fname) as file:
		pass
except IOError as e:
	print('Unable to open file %s'%(fname))
	exit()
cost = genfromtxt(fname,skip_header=3);

# Gradient norm
fname='GradientNormAdjoint_%s.dat'%(sys.argv[1])
try:
	with open(fname) as file:
		pass
except IOError as e:
	print('Unable to open file %s'%(fname))
	exit()
grad = genfromtxt(fname,skip_header=1);

# Regularisation
fname='CostReg_%s.dat'%(sys.argv[1])
try:
	with open(fname) as file:
		pass
except IOError as e:
	print('Unable to open file %s'%(fname))
	exit()
reg = genfromtxt(fname,skip_header=3);
with open(fname, 'r') as f:
    first_line = f.readline().strip()
print(first_line.split((',')))

# get regularisation patrameter
line = open(fname, "r").readlines()[1]
lreg=float(line.split((','))[1])
print(lreg)

# Make the report
with PdfPages('Report_%s.pdf'%(sys.argv[1])) as pdf:
	fig = plt.figure()
	fig.suptitle("Convergence plots", fontsize=16)
	plt.subplot(3,1,1)
	plt.semilogx(cost[:,0],cost[:,2])
	plt.ylabel('rms (m/a)')

	plt.subplot(3,1,2)
	plt.loglog(cost[:,0],cost[:,1],label='$J_0$')
	plt.loglog(reg[:,0],lreg*reg[:,1],label=('%1.2e$ * J_{reg}$'%(lreg)))
	plt.legend(fontsize='10')
	plt.ylim((0.01*np.amin(cost[:,1]),5*np.amax(cost[:,1])))
	plt.ylabel('$J$')

	plt.subplot(3,1,3)
	plt.loglog(grad[:,0],grad[:,1]/grad[0,1])
	plt.xlabel('iter. number')
	plt.ylabel(r'$|\!|g|\!|/|\!|g_0|\!|$')
	plt.tight_layout()
	fig.subplots_adjust(top=0.88)
	pdf.savefig()
	plt.close()

	fig = plt.figure()
	fig.suptitle("M1QN3 last 7 lines", fontsize=16)
	fig.text(.1,.5,output)
	pdf.savefig()  # saves the current figure into a pdf page
	plt.close()
