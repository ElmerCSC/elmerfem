#!/usr/bin/python
#========================================== 
#
# FILE:     makemoulin.py
# USAGE:    python makemoulin.py --meshdir mesh_dir --moulin moulin_file --partition number_of_partition 
# DESCRIPTION:  Add Moulins to partition mesh  
# BUGS: ---
#
# AUTHOR:   mchekki
# ORGANIZATION: CNRS 
#
# VERSION: V0 
# CREATED:  2017-09-05 16:00:35
# MODIFIED: 2017-11-19 08:49:03
#
#========================================== 
import numpy as np
import sys
import os
import getopt

def get_opts(options):
	try:
	    opts, args = getopt.getopt(options,"m:o:p:hv",["meshdir=","moulin=","partition=","help","verbose"])
	except getopt.GetoptError, err:
	    print str(err)
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
        print "USAGE : python makempoulin.py --meshdir mesh_dir --moulin moulin_file --partition number_of_partition"


def exit_error(message):
        print "\nERROR: ", message
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
			exit_error("Directory does not exit: Please use the  ElmerGrid command for mesh partitionning")
	else:
		mesh_dir = '%s'%(mesh)
		if file_donot_exists(mesh_dir):
			exit_error("Directory does not exit")

	NodeMoulinFound=0
	# Read global node file and get the number of all nodes
	nodeAll=np.loadtxt('%s/mesh.nodes'%(mesh))
	NnodeAll=np.size(nodeAll,0)
	MoulinAssign=np.zeros(NnodeAll)

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
	
	    # Read the moulin coordinates
	    Moulin=np.loadtxt(moulin_file)
	    # Open the mesh.nodes file to find the nodes number
	    nodes=np.loadtxt(nodes_file)
	    # Open the mesh.elements
	    elements=np.loadtxt(elements_file)
	    # Read the BC file
	    bc=np.loadtxt(bc_file)
	    # Read the header file
	    head=open(header_file)
	    header=[]
	    for row in head:
	       header.append(np.fromstring(row, sep=' '))

	    if Moulin.size == 2:
		Moulin=Moulin.reshape(1,2)	
	    ms=np.size(Moulin,0)
	    Nnode=np.size(nodes,0)
	    nBC=np.size(bc,0)
	    nheader=np.size(header,0)
	
	    NodeMoulin=np.zeros(ms) 


	    for ii in np.arange(ms):
	        xm = Moulin[ii,0] 
	        ym = Moulin[ii,1]
	        for jj in  np.arange(Nnode):
	           xn = nodes[jj,2] 
	           yn = nodes[jj,3] 
	           if (np.abs((xm-xn)*(xm-xn))+np.abs((ym-yn)*(ym-yn))) < Err:
		       idx = nodes[jj,0]-1	
		       if MoulinAssign[idx] == 0:
	               	  NodeMoulin[ii] = nodes[jj,0] 
		       	  MoulinAssign[idx] = kk+1 
		          NodeMoulinFound=NodeMoulinFound+1
		       else:
	             	  print('WARNING part.%d: Remove Moulin node %d already assigned to partition %d '%(kk+1,idx+1,MoulinAssign[idx]))
			
	    # Test if each Moulin has been associated a mesh node
	    if np.count_nonzero(NodeMoulin)==0:
	             print('WARNING: No moulin node found on partition: %d'%(kk))
	    else:
		    #print("Partition: %d"%(kk))
		    #print NodeMoulin
	
		    if nbpartition >1:
			ms_partition=np.count_nonzero(NodeMoulin)
		
		    # Write the 101 BC at the end of the mesh.boundary file
		    MaxBC = np.max(bc[:,1]) 
		    # get the minimun elment index to add fro 101 BC
		    MinEIndex= np.min(elements[:,0])
		
		    # Rewrite the file and add the 101 elements
		    fid1=open(bc_file,'w');
		    for ii in np.arange(nBC):
		        fid1.write(" ".join(["%g"%(v) for v in bc[ii,:]]))
		        fid1.write("\n")
		    jj=0
		    for ii in np.arange(ms):
			if (NodeMoulin[ii] >0):
		        	jj = jj + 1 
		        	fid1.write('%g %g %g %g %g %g \n'%(nBC+jj,NodeMoulin[ii],MinEIndex,0,101,NodeMoulin[ii]))
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
		exit_error(' No moulin node found on all partitions ')
