#!/usr/bin/python
"""
 slow and simple  parser/hashtable hack for computing f90 makefile dependencies
 usage:
 makedepf95 *.f90 (and cross you fingers it'll work)

 or with elmer:
 makedepf95 *.src
 
 juha.vierinen@csc.fi
"""
import sys
import re
import os

ELMER=True

class SimpleF90Parser:
    """ A simple parser that known only use and module """
    def __init__(self):
        self.fileToObject = re.compile('^([a-zA-Z0-9_-]+).[a-zA-Z0-9]+$')
        
    def parseFile(self,file):
        """ return a tuple (define_modules,used_modules) """
        self.inComment=False
        self.inUseStatement=False
        self.inModuleStatement=False
        self.history=" "
        self.separator = re.compile('[\t, ;]')
        self.comment = re.compile('!')
        self.newline = re.compile('\n')
        self.usest = re.compile('use[ \t]+([a-zA-Z0-9_]+)[ \t,\n;:]',re.I)
        self.modulest = re.compile('module[ \t]+([a-zA-Z0-9_]+)[ \t]*\n$',re.I)
        self.endst = re.compile('end',re.I)
        self.prevWord = re.compile('[\t, ;\n]+([a-zA-Z0-9_]+)[\t, ;\n]+$')
        
        self.usedModules={}
        self.definedModules={}
        input = open(file, 'r')
        self.reset()
        
        a=input.read(1)
        while a:
            self.history=self.history+a

            if self.separator.search(a):
                self.token()
                    
            if self.comment.search(a):
                self.inComment=True
                
            if self.newline.search(a):
                self.token()
                self.reset()

            a=input.read(1)
        
        input.close()
        return (self.usedModules,self.definedModules)


    def token(self):
        if len(self.history) > 3 and not self.inComment:
            p=self.prevWord.search(self.history)
            if p:
                token=p.group(1).lower()
#                print token
                if self.endst.search(token):
                    self.inComment = True

            if self.modulest.search(self.history):
                p=self.modulest.search(self.history)
                if p:
                    self.definedModules[token]=1
                else:
                    sys.stderr.write("Warning: no module\n")

            if self.usest.search(self.history):
                p=self.usest.search(self.history)
                token=p.group(1).lower()
                if p:
                    self.usedModules[token]=1
                else:
                    sys.stderr.write("Warning: no use\n")

            #print "no token found"


    def reset(self):
        self.inComment=False
        self.inUseStatement=False
        self.inModuleStatement=False
        self.history=""

    def getObjectName(self,fname):
        """ simple regexp to create object name from blah.f90 """
        # get object name 
        om=self.fileToObject.search(fname)
        if om:
            oname = om.group(1) + ".$(OBJEXT)"
            return oname
        else:
            print "fatal, name couldn't be created " , fname
            return None

if __name__ == "__main__":
#    depfinder = F90DepFinder(sys.argv[1:])
    p = SimpleF90Parser()


    modtofile = {}
    fileuses = {}
    files=[]
    
    for f in sys.argv[1:]:
        sys.stderr.write(f + "...\n")   
        (um,dm)=p.parseFile(f)

        for module in dm:
            modtofile[module]=f
        
        uses=[]
        fileuses[f]=uses
        
        files.append(f)

        for module in um:
            uses.append(module)
            
    for f in files:
        object=p.getObjectName(f)

        if ELMER:
            line= object + ": " + re.sub('.src','.f90',f)
        else:
            line= object + ": " + f
            
        for use in fileuses[f]:
            if modtofile.has_key(use):
                depobj=p.getObjectName(modtofile[use])
                if object != depobj:
                    line=line + " " + depobj
                else:
                    sys.stderr.write("Warning, removing circular dependency from: "+object+"\n")
            else:
                sys.stderr.write("Warning, module "+use+" undefined\n")
                
        print line
    

