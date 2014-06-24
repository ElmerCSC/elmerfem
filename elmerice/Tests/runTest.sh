#!/bin/bash
#

###------------------------------------------------------------------------------###
# Def extensions used
###------------------------------------------------------------------------------###
if test "$SHL_EXT" = ""; then
    SHL_EXT=".so"
fi
if test "$OBJ_EXT" = ""; then
    OBJ_EXT=".o"
fi

if test "$EXE_EXT" = ""; then
    EXE_EXT=""
fi

export SHL_EXT
export OBJ_EXT
export EXE_EXT

###------------------------------------------------------------------------------###
# PATH ?
###------------------------------------------------------------------------------###
if test "$ELMER_HOME" = ""; then
    # assume that we are testing local version
    printf "\$ELMER_HOME undefined, setting it to ../src\n"
    export ELMER_HOME="../../src"
    export ELMER_LIB="../../src"
    export ELMER_INCLUDE="../../src"
    # elmergrid must be in path
    export ELMER_GRID="ElmerGrid"
    export ELMER_SOLVER="../../src/ElmerSolver"
    # elmergrid must be in path
    export ELMER_MESH2D="Mesh2D"
    export LD_LIBRARY_PATH=".:$ELMER_HOME:$ELMER_HOME/modules:$LD_LIBRARY_PATH"
    # assume that stuff has been installed here
    export PATH=/home/ltavard/Elmer/bin:$PATH
else 
    # ELMER_HOME is defined, so we'll just use that
    export ELMER_HOME=`echo $ELMER_HOME | sed 's+.:+/&+' | sed 's.:..' | sed 's.\\\./.g'`
    export ELMER_LIB="$ELMER_HOME/share/elmersolver/lib"
    export ELMER_INCLUDE="$ELMER_HOME/share/elmersolver/include"
    export ELMER_GRID="$ELMER_HOME/bin/ElmerGrid"
    export ELMER_SOLVER="$ELMER_HOME/bin/ElmerSolver"
    export ELMER_MESH2D="$ELMER_HOME/bin/Mesh2D"
    export LD_LIBRARY_PATH=".:$ELMER_HOME/lib:$LD_LIBRARY_PATH"
    export PATH=$ELMER_HOME/bin:$PATH
fi

export LD=elmerld
export FORTnosh=elmerf90-nosh
export FORT=elmerf90
export LIBS="-L$ELMER_LIB "
echo ""
echo "============================================================="
echo "                    Test ElmerIce Solvers                    "


###------------------------------------------------------------------------------###
# Options: help, all tests, selection repertory, selection solver, default case
###------------------------------------------------------------------------------###

#Delete old directory
if (test -e File_log); then
    rm -r File_log
fi

while true; do
    case "$1" in
	"")
	    echo ""
	    echo "Option [$1]: All tests without email option"
	    echo "============================================================="
	    echo ""
	    dirs=`find . -type d |sed -e 's/\// /g' |awk '{print $2}' |uniq| grep -v .svn | sort`
	    write_mail=NO
	    break
	    ;;
	#Check for optional argument.
	--all | -a) 
	    echo ""
	    echo "Option [$1]: All tests"
	    echo ""
	    echo "Give your mail to send result [N or email]"
	    read courriel
	    echo "============================================================="
	    echo ""
	    if (test $courriel != "N"); then
		write_mail=YES
	    fi
	    dirs=`find . -type d |sed -e 's/\// /g' |awk '{print $2}' |uniq| grep -v .svn | sort`
	    break
	    ;;
	--selection | -s) 
	    echo ""
	    echo "Option [$1]: Selection repertory"
	    echo ""
	    echo "Enter selection test-name with Space-separator"
	    read dirs
	    for dir in $dirs; do
		if (test ! -e $dir); then
		    echo ""
		    printf "Test case [%s] does not found \n" $dir
		    echo "============================================================="
		    echo ""
		    exit
		fi
	    done
	    echo ""
	    echo "Give your mail to send result [N or email]"
	    read courriel
	    echo "============================================================="
	    echo ""
	    if (test $courriel != "N"); then
		write_mail="YES"
	    fi
	    break
	    ;;
	--solvers | -so)
	    echo ""
	    echo "Option [$1]: Find repertory with particular solver"
	    echo ""
	    echo "Enter solver you want to test"
	    read solver
	    echo ""
	    echo "Give your mail to send result [N or email]"
	    read courriel
	    echo "============================================================="
	    echo ""
	    if (test $courriel != "N"); then
		write_mail="YES"
	    fi
	    dirs=`find . -name "*.sif" -print | xargs grep "$solver" |awk -F / '{print$2}'|uniq`
	    rm ListingSolver.txt
	    break
	    ;;
	--help | -h)
	    echo ""
	    echo "4 options: --all       | -a  -> Make all the tests"
	    echo "                             -> Option send-mail include"
	    echo "           --selection | -s  -> Make a selection of tests"
	    echo "                             -> Option send-mail include"
	    echo "           --solver    | -so -> Find repertories with this solver"
	    echo "                             -> Option send-mail include"
	    echo "           ''                -> Default option: all tests"
	    echo "           --help      | -h"
	    echo "============================================================="
	    echo ""
	    exit
	    ;;
	*)
	    echo "runTest.sh: usage error: unrecognize option"
	    echo "usage: runtest.sh --help"
	    exit
	    ;;
    esac
    shift
done

###------------------------------------------------------------------------------###
# Elmer has been installed
# Test if ElmerSolver_mpi has been intalled
# If "yes", parrallel tests will run else parallel runs are ignored
###------------------------------------------------------------------------------###
if ( ! test -e $ELMER_HOME/bin/ElmerSolver_mpi); then
    List_parallel_run=`find . -name "Makefile" -print | xargs grep "ElmerSolver_mpi" . | awk -F / '{print$2}'`
    for var in $List_parallel_run ; do
	dirs=`echo $dirs | sed -e "s/$var//"`
    done
    List_Ignored=$List_parallel_run
    echo "ElmerSolver_mpi not found, parallel test-cases are switched off:"
    echo $List_parallel_run
    echo ""
fi

###------------------------------------------------------------------------------###
# Test if gmsh has been intalled
# If "No", tests with gmsh are ignored
###------------------------------------------------------------------------------###
if ( test $(which gmsh | wc -l) -eq 0) ; then
    List_gmsh_run=`find . -name "Makefile" | xargs grep "GMSH" . |   awk -F / '{print$2}'`
    for var in $List_gmsh_run ; do
	dirs=`echo $dirs | sed -e "s/$var//"`
    done
    echo "gmsh not found, tests with gmsh are switched off:"
    echo $List_gmsh_run
    echo ""
else
    export GMSH=$(which gmsh)
fi


###------------------------------------------------------------------------------###
# Loop on repertory to realise each test
###------------------------------------------------------------------------------###

passed=0
nbtest=0
visited=0
iter_List=0
List_Failed=( )

# Compilation f90 -
if (test -e Compare.f90); then
    elmerf90-nosh Compare.f90 -o Compare >  /dev/null
    chmod 775 Compare$EXE_EXT
else
    echo "File Compare.f90 not found\n"
    exit
fi

# Test case
for dir in $dirs; do

    #initialisation
    Error=0
    n_files=0
    success=0
    Result=0
    Result_True=0

    #counter test
    nbtest=`expr $nbtest + 1`
    cwd=`pwd`

    # Move in repertory
    cd $dir

    # Delete old files
    if (test -e difference.txt); then
	rm difference.txt
    fi
    if (test -e OutputSIF_$dir.log); then
	rm OutputSIF_$dir.log
    fi
    if (test -e Output_$dir.log); then
	rm Output_$dir.log
    fi

    printf "test $nbtest : %25s " $dir
    
    # Run .sif and verification "all done"
    make run  > OutputSIF_$dir.log 2>&1

    if (test -e OutputSIF_$dir.log); then
	success=`grep "ALL DONE" OutputSIF_$dir.log | wc -l`
	n_files=1
	
	if (test $success -ge 1) ; then
	    
	    if (test -f results.dat); then
    	    # Transform file with only one " " between colomns
		sed -e 's/^ *//' results.dat > $dir.txt
		
		if (test -f $dir.txt); then
		    
		    if (test -e results.dat.names ); then
			n_arg=`tail -1 results.dat.names |  awk '{print $1}'| sed 's/://'`
			n_line_cpu=`grep -n 'cpu time' results.dat.names | awk '{print $1}' | sed 's/://'`
		    else
			Error=5
		    fi
		else
		    Error=4
		fi
	    else
		Error=3
	    fi
	else
	    Error=2
	fi
    else
	Error=1
    fi
    
    if (test $Error -eq 0) ; then #All conditions checked for comparison
	if (test ! $n_arg -eq 0); then 
            # Test how many arguments for comparison	
	    if  (test -n "$n_line_cpu" ); then
		nb_line_tot=`wc -l results.dat.names |awk '{print $1}'`
		n_col_cpu=`head -$nb_line_tot results.dat.names |tail -1| awk '{print $1}' | sed 's/://'`
	    else 
		n_col_cpu=0
	    fi
	    # Tests if exists specific TARGET in *.sif else TARGET's value is 1E-6
	    if (test `grep "TARGET" OutputSIF_$dir.log | uniq | awk -F = '{print$2}'`=""); then
		TARGET=1e-6
	    else
		TARGET=`grep "TARGET" OutputSIF_$dir.log | uniq | awk -F = '{print$2}'`
	    fi
            # Results of comparison in file difference.txt and errors print in Output_$dir.log
	    `../Compare$EXE_EXT $dir valid_$dir $n_arg $n_col_cpu $TARGET > Output_$dir.log`
	    FindError=`grep "ERROR" Output_$dir.log | wc -l`
	    FindDiff=`grep "DIFF-TIME" Output_$dir.log | wc -l`
	    # Values for screen output
	    Result=`tail -1 difference.txt | awk '{print $1}'`
	    Result_Valid=`tail -1 difference.txt | awk '{print $2}'`
	    Number_Argument=`tail -1 difference.txt | awk '{print $3}'`
	    Number_concat="$Number_Argument"": "
	    grep "$Number_concat" results.dat.names > temp_name.txt
	    Name_Var=`cut -d : -f 2,3 temp_name.txt`
	    rm temp_name.txt

	    if (test ! -e difference.txt); then
		List_Failed[$iter_List]=$dir
		((iter_List++))
		printf "		        [FAILED]\n"
		printf "File [difference.txt] not found: look at [File_log/Output_%s.log] for details\n" $dir
		if (test $visited -eq 0); then
		    `mkdir ../File_log`
		visited=1
		fi
		
	       #Copy files wich resume SIF, and problem
		cp ./Output_$dir.log ../File_log
		cp ./OutputSIF_$dir.log ../File_log
		rm ./Output_$dir.log  ./OutputSIF_$dir.log

            #If bad results
	    elif (test !  $FindError -eq 0) ; then
		cat Output_$dir.log > tmp01.log
   		cat results.dat.names tmp01.log > Output_$dir.log
		rm tmp01.log 
		
		List_Failed[$iter_List]=$dir
		((iter_List++))
		printf "		        [FAILED]\n"
		printf "Differences in results: look at [File_log/Output_%s.log] for details\n" $dir
		printf "For variable [%s]:  Found result [%s] - Expected result [%s]\n" "$Name_Var" "$Result" "$Result_Valid"
	
		if (test $visited -eq 0); then
		    mkdir ../File_log
		    visited=1
		fi
                #Copy files wich resume difference
		cp ./Output_$dir.log ../File_log
		cp ./OutputSIF_$dir.log ../File_log
		rm ./Output_$dir.log  ./OutputSIF_$dir.log
		
            # Difference between cpu-time is observed, test PASSED
	    elif (test ! $FindDiff -eq 0 ); then
		passed=`expr $passed + 1`
		printf "		[PASSED - Diff Time-]\n"
	    	
            # Nothing appened, test PASSED
	    else	
		passed=`expr $passed + 1`
		printf "		[PASSED]\n"
	    fi
	    make -i clean > /dev/null 2>&1
	else
	    printf "		        [FAILED]\n"
	    printf "Problem in number of arguments for comparison\n"
	    make -i clean > /dev/null 2>&1
	fi
    else
	case $Error in
	    1)
		printf "		        [FAILED]\n"
		printf "Simulation aborted: OutputSIF_$dir.log not found\n";;
	    2)
		printf "		        [FAILED]\n"
		printf "Simulation aborted: look at [File_log/OutputSIF_%s.log] for details\n" $dir;;
	    3)
		printf "		        [FAILED]\n"
		printf "File [results.dat] not found: look at [File_log/OutputSIF_%s.log] for details\n" $dir;;
	    4)
		printf "		        [FAILED]\n"
		printf "Problems with transforming results.dat in $dir.txt\n";;
	    5)
		printf "		        [FAILED]\n"
		printf "File [results.dat.names] not found: look at [File_log/OutputSIF_%s.log] for details\n" $dir;;
	    
	esac
	List_Failed[$iter_List]=$dir
	((iter_List++))
	if (test $visited -eq 0); then
	    mkdir ../File_log
	    visited=1
	fi
        #Copy files wich resume difference
	if (test -e OutputSIF_$dir.log); then
	    cp ./OutputSIF_$dir.log ../File_log
	    rm ./OutputSIF_$dir.log
	fi
	if (test -e Output_$dir.log); then
	    cp ./Output_$dir.log ../File_log
	    rm ./Output_$dir.log
	fi
	make -i clean > /dev/null 2>&1
    fi
    
    # Then come-back to global repertory
    cd $cwd
done

# Mail-option
if [ '$write_mail' = 'YES' ] && [ $visited -eq 1 ]; then
    `tar -zcvf File_log.tar.gz File_log | uuencode File_log.tar.gz | mail -s "Resume" $courriel`
    rm File_log.tar.gz
elif [ '$write_mail' = 'YES' ]; then
    echo "All tests passed" | mail -s "Resume Test" $courriel
fi

#Print Screen
echo ""
printf "Tests completed, passed: $passed out of total $nbtest tests \n"
echo ""
