#!/usr/bin/python
#========================================== 
#
# FILE:     makemoulin.py
# USAGE:    python makemoulin.py --meshdir mesh_dir --moulin moulin_file --partition number_of_partition 
# DESCRIPTION:  Add Moulins to partition mesh  
#
# BUGS: ---
#       20/10/2019: Correct indexes that were read as float
#       10/12/2019: Use panda to read the mesh files instead of numpy 
#                   Columns do not have the same size in the mesh.boundary file 
#                   This leads to an error 
#       22/12/2019: Add more info whether the Moulins has distributed or not
#                    between partitions 
#
# AUTHOR:   mchekki
# ORGANIZATION: CNRS 
#
# VERSION: V1 
# CREATED:  2017-09-05 16:00:35
# MODIFIED: 2019-12-23 12:40:21
#
#========================================== 
import numpy as np
import sys
import os
import getopt
import pandas as pd

def get_opts(options):
	try:
		opts, args = getopt.getopt(options,"m:o:p:hv",["meshdir=","moulin=","partition=","help","verbose"])
	except getopt.GetoptError as err:
		print (err)
		usage()
		sys.exit(2)
	nbpartition =1	
	for o, a in opts:
		if o in ("-m", "--meshdir"):
			mesh= a
		elif o in ("-o", "--moulin"):
			moulin_file= a
		elif o in ("-p", "--partition"):
			nbpartition= np.int(a)
		elif o in ("-h", "--help"):
			usage()
			sys.exit(2)
		else:
			assert False, "unhandled option"
	
	return mesh, moulin_file, nbpartition 
def usage():
        print ('USAGE : python makempoulin.py --meshdir mesh_dir --moulin moulin_file --partition number_of_partition')


def exit_error(message):
        print('\nERROR: ',message)
        sys.exit(1)

def file_exists(filename):
        if os.path.exists(filename):
                return True
        else:
                return False

def file_donot_exists(filename):
        if os.path.exists(filename) == False:
                return True
        else:
                return False

if __name__=='__main__':

	options=sys.argv[1:]
	mesh, moulin_file, nbpartition = get_opts(options)

	Err = 0.01 
	# Define the mesh directory 
	if nbpartition >1:
		mesh_dir = '%s/partitioning.%s'%(mesh,np.str(nbpartition))
		if file_donot_exists(mesh_dir):
			exit_error("Directory %s does not exit: Please use the  ElmerGrid command for mesh partitionning"%(mesh_dir))
	else:
		mesh_dir = '%s'%(mesh)
		if file_donot_exists(mesh_dir):
			exit_error("Directory %s does not exit"%(mesh_dir))
	if file_donot_exists(moulin_file):
		exit_error("Moulin file %s does not exit"%(moulin_file))

	# Read global node file and get the number of all nodes
	# Depending on the mesh, columns do not have the same size
	# In order to fix this , fill columns with -999999 value

	nodeAll=pd.read_table('%s/mesh.nodes'%(mesh), dtype=float , header=None, sep='\s+').fillna(-999999).values
	NnodeAll=np.size(nodeAll,0)
	MoulinAssign=np.zeros(NnodeAll)
	# Read the moulin coordinates
	Moulin=pd.read_table(moulin_file, dtype=float , header=None, sep='\s+').fillna(-999999).values
	if Moulin.size == 2:
	    Moulin=Moulin.reshape(1,2)	
	ms=np.size(Moulin,0)
	NodeMoulin=np.zeros(ms) 
	gbc_file='%s/mesh.boundary'%(mesh)
	gbc=pd.read_table(gbc_file, dtype=float , header=None, sep='\s+').fillna(-999999).values
	MaxBC = np.max(gbc[:,1]) 

	for kk in np.arange(nbpartition): 
		if nbpartition >1:
			nodes_file='%s/part.%s.nodes'%(mesh_dir,np.str(kk+1))
			elements_file='%s/part.%s.elements'%(mesh_dir,np.str(kk+1))
			bc_file='%s/part.%s.boundary'%(mesh_dir,np.str(kk+1))
			header_file='%s/part.%s.header'%(mesh_dir,np.str(kk+1))
		else:
			nodes_file='%s/mesh.nodes'%(mesh_dir)
			elements_file='%s/mesh.elements'%(mesh_dir)
			bc_file='%s/mesh.boundary'%(mesh_dir)
			header_file='%s/mesh.header'%(mesh_dir)
		# Open the mesh.nodes file to find the nodes number
		nodes=pd.read_table(nodes_file, dtype=float , header=None, sep='\s+').fillna(-999999).values
		# Open the mesh.elements
		elements=pd.read_table(elements_file, dtype=float , header=None, sep='\s+').fillna(-999999).values
		# Read the BC file
		bc=pd.read_table(bc_file, dtype=float , header=None, sep='\s+').fillna(-999999).values
		# Read the header file
		head=open(header_file)
		header=[]
		for row in head:
		   header.append(np.fromstring(row, sep=' '))

		if np.int(header[2][0]) == 101 :
			exit_error("Found 101 elements in the mesh header. No need to proceed...")
		
		Nnode=np.size(nodes,0)
		nBC=np.size(bc,0)
		nheader=np.size(header,0)

		NodeMoulinFound=0
		NodeMoulin=np.zeros(ms) 
		for ii in np.arange(ms):
			xm = Moulin[ii,0] 
			ym = Moulin[ii,1]
			for jj in  np.arange(Nnode):
				xn = nodes[jj,2] 
				yn = nodes[jj,3] 
				if (np.abs((xm-xn)*(xm-xn))+np.abs((ym-yn)*(ym-yn))) < Err:
					idx = np.int(nodes[jj,0]-1)
					if MoulinAssign[idx] == 0:
						# Save Moulin nodes 
						if NodeMoulin[ii] == 0:
							NodeMoulin[ii] = nodes[jj,0] 
						NodeMoulinFound=NodeMoulinFound+1
						MoulinAssign[idx] = kk+1 
					else:
						print('WARNING part.%d: Remove Moulin node %d already assigned to partition %d '%(kk+1,idx+1,MoulinAssign[idx]))
			 # Test if each Moulin has been associated a mesh node
		if np.count_nonzero(NodeMoulin)==0:
			print('WARNING: No moulin nodes found on partition: %d'%(kk+1))
		else:
			
			ms_partition=np.count_nonzero(NodeMoulin)
			print('%d moulin nodes found on partition %d'%(ms_partition, kk+1))
			
			# Write the 101 BC at the end of the mesh.boundary file
			# get the minimum element index to add for 101 BC
			MinEIndex= np.min(elements[:,0])
			
			# Rewrite the file and add the 101 elements
			fid1=open(bc_file,'w');
			for ii in np.arange(nBC):
			    boundline = bc[ii,:]
			    # Remove columns with -999999 before writing
			    indexes = np.where( boundline != -999999)
			    fid1.write(" ".join(["%g"%(v) for v in boundline[indexes]]))
			    fid1.write("\n")
			jj=0
			for ii in np.arange(ms):
			    if (NodeMoulin[ii] >0):
			    	jj = jj + 1
			    	fid1.write('%g %g %g %g %g %g \n'%(nBC+jj,MaxBC+1+ii,MinEIndex,0,101,NodeMoulin[ii]))
			fid1.close()
			
			# Change the header file
			fid1=open(header_file,'w')
			if nbpartition >1:
				fid1.write('%g %g %g \n'%(header[0][0],header[0][1],header[0][2]+ms_partition))
			else:
				fid1.write('%g %g %g \n'%(header[0][0],header[0][1],header[0][2]+ms))
			fid1.write('%g \n'%(header[1][0]+1))
			if nbpartition >1:
				fid1.write('%g %g \n'%(101,ms_partition))
			else:
				fid1.write('%g %g \n'%(101,ms))
			for ii in np.arange(2,nheader):
			    fid1.write('%g %g \n'%(header[ii][0],header[ii][1]))
			
			fid1.close()
	if NodeMoulinFound == 0:
		exit_error(' ERROR: No moulin node found on all partitions ')
