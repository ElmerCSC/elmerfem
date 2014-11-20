#!/usr/bin/python
import os, fnmatch, re, types, subprocess, sys

DEFAULT_BUILD_DIR="/home/jkataja/src/elmer/build"

##################################

ignore_file_glob = ['*.swp','Makefile','*.cmake',\
        'CMakeLists.txt', '*.o','*.sif']
source_file_glob = ['*.f90', '*.F90']

CMakeLists_Template = """\
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
%s

file(COPY%s DESTINATION \"${CMAKE_CURRENT_BINARY_DIR}/\")

ADD_ELMER_TEST(%s)
"""

runtests_Template = """\
include(${TEST_SOURCE}/../test_macros.cmake)
%s
RUN_ELMER_TEST()
"""

#failed_tests = []
#with open("/home/jkataja/Desktop/failed",'r') as failedfile:
    #failed_tests = map(str.rstrip, failedfile.readlines())
    #failedfile.close()

##################################

def contains(Elem, Set, Match):
    for set_elem in Set:
        if(Match(Elem,set_elem)):
            return True
    return False

def find_elems(Elem, List, Match):
    keys = []
    for m in range(0,len(List)):
        if Match(List[m], Elem):
            keys.append(m)
    return keys

def yesnomatch(string, patterns):
    Match = fnmatch.fnmatchcase
    if Match(string,patterns[0])==True and Match(string,patterns[1]) == False:
        return True
    return False

def parse_makefile(makefilename):
    dict_pattern = re.compile("^([a-zA-Z]+):$")
    rule_pattern = re.compile("^\t(.+)$")
    makefiledict = {}
    dict_key = ""
    with open(makefilename, 'r') as makefile:
        for line in makefile: 

            dict_key_ob = re.match(dict_pattern, line)
            rule_key_ob = re.match(rule_pattern, line)

            if not(dict_key_ob is None):
                dict_key = dict_key_ob.groups()[0]
                if not(makefiledict.has_key('dict_key')):
                    makefiledict[dict_key] = []

            if not(rule_key_ob is None):
                rule = rule_key_ob.groups()[0]
                key_rules = makefiledict[dict_key]
                key_rules.append(rule)
                makefiledict[dict_key] = key_rules
    return makefiledict


def increase_solver_n(testname, decrease=False):
    norm_re = re.compile('^(Solver )(\d)+( :: Reference Norm = Real )(\S*)')
    tol_re = re.compile("^(Solver )(\d)+( :: Reference Norm Tolerance = Real )(\S*)")
    filename = testname + "/" + get_sif_name(testname)
    with open(filename,'r') as infile:
        lines = infile.readlines()
        infile.close() 
        for m in range(0,len(lines)):
            re_obj = norm_re.search(lines[m]) 
            if re_obj != None:
                gr = re_obj.groups()
                print(lines[m])
                if not decrease:
                    newline = gr[0] + str(int(gr[1])+1) + gr[2] + gr[3]
                else:
                    newline = gr[0] + str(int(gr[1])-1) + gr[2] + gr[3]
                lines[m] = newline+'\n'
                print(lines[m])

            re_obj = tol_re.search(lines[m])
            if re_obj != None:
                gr = re_obj.groups()
                print(lines[m])
                if not decrease:
                    newline = gr[0] + str(int(gr[1])+1) + gr[2] + gr[3]
                else:
                    newline = gr[0] + str(int(gr[1])-1) + gr[2] + gr[3]
                lines[m] = newline+"\n"
                print(lines[m])
        outfile = open(filename,'w')
        outfile.writelines(lines)
        return lines


def execute_test(testname, elmer_build_dir=DEFAULT_BUILD_DIR, remake=True):
    build_dir = elmer_build_dir + "/fem/tests"
    if remake:
        makeval = subprocess.call(["make","-j4"], cwd=build_dir+"/../")
    process = subprocess.call("ctest", cwd=build_dir+"/"+testname)
    return process == 0


def convert_sif_file(filename):
    NRM_re = re.compile("NRM=[ ]*([+|-]?\d+\.?\d*[e|E]?[+|-]?\d*)")
    EPS_re = re.compile("EPS=[ ]*([+|-]?\d+\.?\d*[e|E]?[+|-]?\d*)")
    with open(filename,'r') as infile:
        lines = infile.readlines()
        TESTCASE_lines = find_elems(["*TEST CASE*","*END*"], lines, yesnomatch)
        NRM_lines = find_elems("*NRM=*", lines, fnmatch.fnmatchcase)
        if len(TESTCASE_lines) != len(NRM_lines):
            raise Exception("nonequal amount of TESTCASE and NRM keywords in %s file" % (filename))
        infile.close() 
        #infile = open(filename,'w')
        target_norms = [[(NRM_re.search(lines[n]).groups()[0]),None] for n in NRM_lines]
        for m in range(0,len(NRM_lines)):
            re_obj = EPS_re.search(lines[NRM_lines[m]]) 
            if re_obj != None:
                target_norms[m][1] = re_obj.groups()[0]
        offset = 0
        for r in range(0,len(TESTCASE_lines)):
            lines.insert(TESTCASE_lines[r]+offset, "Solver 1 :: Reference Norm = Real %s\n" % target_norms[r][0])
            offset = offset + 1
            if target_norms[r][1] != None:
                lines.insert(TESTCASE_lines[r]+offset, "Solver 1 :: Reference Norm Tolerance = Real %s\n" % target_norms[r][1])
                offset = offset + 1
        outfile = open(filename,'w')
        outfile.writelines(lines)
        return (NRM_lines, TESTCASE_lines, lines)


    raise Exception("convert_sif_file is not ready")

def get_sif_name(testname):
    with open(testname + "/" + "ELMERSOLVER_STARTINFO") as startinfofile:
        siffilename = startinfofile.readline().rstrip('\n')
        startinfofile.close()
        return(siffilename)

def convert_test(testname):
    """ Converts old test to follow new testing paradigm """
    print("Converting " + testname)
    cmakelists_file = open(testname+"/CMakeLists.txt", 'w')
    runtests_file = open(testname+"/runtest.cmake", 'w')
    makefile_dict = parse_makefile(testname+"/Makefile")

    add_module_pattern = """ADD_ELMERTEST_MODULE(%s %s %s)"""
    module_string = ""

    test_files = os.listdir(testname)
    modules = []
    copyfiles = ""

    # Find files to copy and transform the .sif file
    for test_file in test_files:
        if not(contains(test_file, ignore_file_glob, fnmatch.fnmatchcase)):
                copyfiles = copyfiles + " " + test_file
        if contains(test_file, source_file_glob, fnmatch.fnmatchcase):
            modules.append([testname,test_file.split('.')[0], test_file])


    # find sif file from ELMERSOLVER_STARTINFO
    siffilename = get_sif_name(testname)

    #    copyfiles = copyfiles + " " + siffilename.rstrip('\n')
    print("> " + testname+"/"+siffilename)
    convert_sif_file(testname+"/"+siffilename)
    module_string = module_string + "\nCONFIGURE_FILE( "+siffilename+" "+siffilename+" COPYONLY)"

    # find modules to compile
    for module in modules:
        module_string = module_string + "\n" + \
                (add_module_pattern % (module[0],module[1],module[2])) 


    cmakelists = CMakeLists_Template % (module_string, copyfiles, testname)

    # parse elmergrid command
    run_test_command = ""
    for run_rule in makefile_dict['run']:
        if(fnmatch.fnmatchcase(run_rule,"$(ELMER_GRID)*")):
            run_test_command = run_test_command + \
                    "execute_process(COMMAND %s)" % \
                    run_rule.replace("$(ELMER_GRID)", "${ELMERGRID_BIN}") + "\n"


    cmakelists_file.write(cmakelists)

    runtests = runtests_Template % (run_test_command)
    runtests_file.write(runtests)

    cmakelists_file.close()
    runtests_file.close()

    return True


def fixtests(testlist):
    max_solver_n = 5
    n = 1
    helpless_tests = []
    for test in testlist:
        while(execute_test(test, remake= n!=1) != True and n != max_solver_n):
            increase_solver_n(test)
            n = n+1
        if n == max_solver_n:
            helpless_tests.append(test)
        n = 1
    return helpless_tests

if __name__ == '__main__':
    if len(sys.argv) > 1:
        try:
            convert_test(sys.argv[1])
        except Exception as e:
            print(e)
            print("Failed to convert %s"% (sys.argv[1]))
    else:
        print("Usage: convert_tests.py <test-directory>")
        print("Convert traditional Elmer test to ctest (EXPERIMENTAL).")
