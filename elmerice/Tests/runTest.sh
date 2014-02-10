#!/bin/sh
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
    #printf "ELMER_HOME=%s\n" $ELMER_HOME
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
	    dirs=`find . -type d |sed -e 's/\// /g' |awk '{print $2}' |uniq| sort`
	    break
	    ;;
	"")
	    echo ""
	    echo "Option [$1]: All tests without email option"
	    echo "============================================================="
	    echo ""
	    dirs=`find . -type d |sed -e 's/\// /g' |awk '{print $2}' |uniq| sort`
	    write_mail=NO
	    break
	    ;;
	--selection | -s) #Check for optional argument.
	    echo ""
	    echo "Option [$1]: Selection repertory"
	    echo ""
	    echo "Enter selection test-name with Space-separator"
	    read dirs
	    for dir in $dirs; do
		if (test ! -e $dir); then
		    printf "Test case %s does not found \n" $dir
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
	    find . -maxdepth 1 -mindepth 1 | grep -r --include="*.sif" $solver > ListingSolver.txt
	    dirs=`cut -d / -f 1 ListingSolver.txt |uniq`
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
# Loop on repertory to realise each test
###------------------------------------------------------------------------------###

passed=0
nbtest=0
visited=0
Error1=0
iter_List=0

# Compilation f90 -
if (test -e Compare.f90); then
    elmerf90-nosh Compare.f90 -o Compare >  /dev/null
    chmod 775 Compare$EXE_EXT
else
    echo "File Compare.f90 not found"
    exit
fi

# Test case
for dir in $dirs; do
    nbtest=`expr $nbtest + 1`
    cwd=`pwd`
    # Move in repertory
    cd $dir

    # Delete old files
    
    if (test -e difference.txt); then #clean old file
	rm difference.txt
    fi
    if (test -e Result_$dir.log); then
	rm Output_$dir.log
    fi

    printf "test $nbtest : %25s " $dir
    
    # Run .sif and verification "all done"
    make run  > OutputSIF_$dir.log 2>&1 
    success=`grep "ALL DONE" OutputSIF_$dir.log | wc -l`
    sleep 10

    if (test -e results.dat); then
	n_files=1
    else
	n_files=0
    fi
    
    if (test $success -ge 1) && (test ! $n_files -eq 0); then
       	# Transform file with only one " " between colomns
	if (test -f results.dat); then
	    sed -e 's/^ *//' results.dat > $dir.txt
	    #`mv results.dat $dir.txt`
	else #Pbm creation file results
	    Error1=1
	fi
	
	# Comparison and print ERREUR if bad results
	if (test -e results.dat.names ); then
	    n_arg=`tail -1 results.dat.names |  awk '{print $1}'| sed 's/://'`
	    n_line_cpu=`grep -n 'cpu time' results.dat.names | awk '{print $1}' | sed 's/://'`
	else
	    n_arg=0
	fi
	
        
	if (test ! $n_arg -eq 0); then 
            # Test how many arguments for comparison
            #      anf if particular files for comparison exists
	    if  (test -n "$n_line_cpu" ); then
		nb_line_tot=`wc -l results.dat.names |awk '{print $1}'`
		n_col_cpu=`head -$nb_line_tot results.dat.names |tail -1| awk '{print $1}' | sed 's/://'`
	    else 
		n_col_cpu=0
	    fi
	    # Results of comparison in file difference.txt and errors print in OutputTest_$dir.log
	    `../Compare$EXE_EXT $dir valid_$dir $n_arg $n_col_cpu > OutputTest_$dir.log`
	    FindError=`grep "ERROR" OutputTest_$dir.log | wc -l`
	    FindDiff=`grep "DIFF-TIME" OutputTest_$dir.log | wc -l`
	fi
    fi
    
    #If comparison failed and difference.txt do not exist (created by Compare.f90)
    if (test ! -e difference.txt); then
	List_Failed[$iter_List]=$dir
	((iter_List++))
	printf "		[FAILED]\n"
	printf "Comparison failed: look at [File_log/%s.log] for details\n" $dir
	echo 'File difference.txt not found' >> Output_$dir.log
	if (test $visited -eq 0); then
	    `mkdir ../File_log`
	    visited=1
	fi
	 #Copy files wich resume SIF, and problem
	cp ./Output_$dir.log ../File_log
	cp ./OutputSIF_$dir.log ../File_log
	
	make -i clean > /dev/null 2>&1

    #If bad results
    elif (test ! $FindError -eq 0) || (test $Error1 -eq 1); then
	if (test -e results.dat.names ); then
	    cat results.dat.names OutputTest_$dir.log > Output_$dir.log
	    rm OutputTest_$dir.log
	fi
	List_Failed[$iter_List]=$dir
	((iter_List++))
	printf "		[FAILED]\n"
	printf "look at [File_log/%s.log] for details\n" $dir
	
	if (test $visited -eq 0); then
	    mkdir ../File_log
	    visited=1
	fi
        #Copy files wich resume difference
	cp ./Output_$dir.log ../File_log

	make -i clean > /dev/null 2>&1
	
    # Difference between cpu-time is observed, test PASSED
    elif (test ! $FindDiff -eq 0 ); then
	passed=`expr $passed + 1`
	printf "		[PASSED - Diff Time-]\n"
	make -i clean > /dev/null 2>&1

    # Nothing appened, test PASSED
    else	
	passed=`expr $passed + 1`
        printf "		[PASSED]\n"
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
if [ ${#List_Failed[@]} != 0 ]; then
    echo "Tests FAILED:"
    for elt in "${List_Failed[*]}"; do echo $elt ;done
fi
echo ""
