#!/usr/bin/python
#
# create preprocessing instructions
# from sources.lst
#
import re

sources="amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults amd_order amd_control amd_info amd_valid amd_preprocess"

amd_sources=[]

for f in re.split(" ",sources):
    nf=re.sub('amd_','amd_i_',f) + ".c"
    amd_sources.append(nf)
    print nf + ":\n\t" + "$(CPP) $(INCLUDES) -DINT " + f + ".c > " + nf

for f in re.split(" ",sources):
    nf=re.sub('amd_','amd_l_',f) + ".c"
    amd_sources.append(nf)
    print nf + ":\n\t" + "$(CPP) $(INCLUDES) -DLONG " + f + ".c > " + nf

print 
print "AMD_CPP_SOURCES = \\"
for f in amd_sources:
    if f == amd_sources[len(amd_sources)-1]:
        print "\t"+f 
    else:
        print "\t"+f + " \\"
