#!/usr/bin/python
#
# create preprocessing instructions
# from sources.lst
#
import re,sys

input=file("cflags.dep")

origSrcs={}
ppSrcs=[]

for l in input:
    d=re.split(" ",l)
    d[0]=re.sub("[ \t]","",d[0])

    flags = ""
    for i in d[2:]:
        i=re.sub("\n","",i)
        flags = flags + " " + i

    print d[1] + ":\n\t" + "$(CPP) $(INCLUDES) " + flags + " " + d[0] + " > " + d[1]
    
    ppSrcs.append(d[1])
    origSrcs[d[0]]=1

print
print "UMFPACK_CPP_SOURCES = \\"
for f in ppSrcs:
    if f == ppSrcs[len(ppSrcs)-1]:
        print "\t"+f 
    else:
        print "\t"+f + " \\"

print
print "UMFPACK_ORIG_SOURCES = \\"
ok=origSrcs.keys()
for key in ok:
    if key == ok[len(ok)-1]:
        print "\t"+key 
    else:
        print "\t"+key + " \\"


