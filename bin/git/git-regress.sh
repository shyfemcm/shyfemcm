#!/bin/sh
#
# regression test to check different commits
#
# compiles different commits and runs shyfem
# optioanlly compares results with external files
#
# specify the commits to be used the variable commits
# alternatively specify first and last commit and optionally step
# ...and call MakeCommitsBetween to start from first with step to last
#
#---------------------------------------------------------

shyfemdir=/home/georg/shyfemcm
simul=venlag
strfile=veniceLagoon.str
strfile=test03_tvd.str
run_command="mpirun -np 4 shyfem $strfile"
make_clean="NO"
do_compare="NO"
default_branch=develop

#---------------------------------------------------------

actdir=$( pwd )
overall=$actdir/overall_diff.txt

#---------------------------------------------------------

MakeCommitsBetween()
{
  local c1=$1
  local c2=$2
  local step=$3

  commits=""

  [ $# -lt 2 ] && echo "need two commits in function" && return 1

  [ -z "$step" ] && step=1
  [ $step -eq 0 ] && step=1

  echo "look for commits between: $c1 $c2" with step $step

  in_commits="NO"
  i=0
  iout=0

  cd $shyfemdir
  commits_all=$( $shyfemdir/bin/git/git-tags -all | sed -e 's/ .*//' )

  for commit in $commits_all
  do
    #echo "commit: $commit"
    if [ $in_commits = "NO" ]; then
      if [ $commit = $c1 ]; then
        in_commits="YES"
        commits="$commit"
        iout=$((iout+1))
      fi
    else
      i=$((i+1))
      n=$((i % step))
      if [ $n -eq 0 ]; then
        iout=$((iout+1))
        commits="$commits $commit"
      fi
      if [ $commit = $c2 ]; then
        if [ $n -ne 0 ]; then                   #always put last commit
          iout=$((iout+1))
          commits="$commits $commit"
        fi
        break
      fi
    fi
  done

  echo "commits gathered: $iout"
}

ShowCommits()
{
  local i=0
  for commit in $commits
  do
    i=$((i+1))
    echo "chosen commits: $i $commit"
  done
}

GetCommit()
{
  cd $shyfemdir
  logfile=$actdir/log_commit_$commit.log
  git co Rules.make
  git co $commit > $logfile 2>&1
  version=$( head -1 VERSION )
  lasttag=$( bin/git/git-tags | tail -1 )
  cd $actdir
}

ResetCommit()
{
  cd $shyfemdir
  git co Rules.make > /dev/null 2>&1
  git co $default_branch    > /dev/null 2>&1
  cd $actdir
}

#---------------------------------------------------------

Info()
{
  echo "version = $version"
  echo "lasttag = $lasttag"
  echo "commit  = $commit "
}

Compile()
{
  cd $shyfemdir
  logfile=$actdir/log_compile_$commit.log
  echo "compiling shyfem"
  make rules_mpi 	 > $logfile 2>&1
  [ $make_clean = "YES" ] && make cleanall		>> $logfile 2>&1
  make			>> $logfile 2>&1
  if [ $? -eq 0 ]; then
    echo "successful compilation"
  else
    echo "$commit   error in compilation"
    echo "$commit   error in compilation" >> $overall
    return 3
  fi
  cd $actdir
  return 0
}

Run()
{
  cd $actdir
  logfile=log_run_$commit.log
  echo "running shyfem with log to $logfile"
  $run_command > $logfile 2>&1
  status=$?
  if [ $status -eq 99 ]; then
    echo "successfull run finished with status = $status"
  else
    echo "$commit   error in run"
    echo "$commit   error in run" >> $overall
    return 5
  fi
  return 0
}

Compare()
{
  if [ $do_compare = "NO" ]; then
    echo "$commit   successfull run"
    echo "$commit   successfull run" >> $overall
    return 0
  fi
  echo "comparing results..."
  cd $actdir
  inffile=$simul.$commit.inf
  mv $simul.inf $inffile
  diff $inffile check/venlag.inf.8_0_8 > diff.txt
  errors8=$( cat diff.txt | wc -l )
  diff $inffile check/venlag.inf.8_0_10 > diff.txt
  errors10=$( cat diff.txt | wc -l )
  echo "---------------------------------"
  echo "differences in lines (8_0_8): $errors8"
  echo "differences in lines (8_0_10): $errors10"
  echo "---------------------------------"
  echo "$commit   $errors8  $errors10" >> $overall
}

DoRegress()
{
  echo "================================="
  GetCommit $commit
  Info
  Compile
  [ $? -eq 0 ] || return
  Run
  [ $? -eq 0 ] || return
  Info
  Compare
  ResetCommit
  echo "================================="
}

#---------------------------------------------------------

commits="d9a5cd98ef d1c9eb47b9"
commits="2934ed564e 4c36c05105"
commits="81d1e91388"
commits="f447e1096d"		#error
commits="4c9115fbff"		#error
commits="b06620e4f1"
commits="00da14193e"
commits="0a3c58cfa7"
commits="26999920db b859e7b1ef 88a2acc0d9"
commits="26999920db 103eef11f0 74e30fed50 39e7f9fab0 d0328bac3d cc93d161c0"
commits="74e30fed50 39e7f9fab0"
commits="e468e1eaa2 fa19003ecc e9111a719e f040e4b4bf 39e7f9fab0"
commits="db3f59f827 7af8214d2f 9b59c5db36 e3f224f34b b8030366d7 39e7f9fab0"
commits="8e812355f9 d9a5cd98ef 4c36c05105 e468e1eaa2 5178686a35"
commits="4c36c05105 4c9115fbff 0a3c58cfa7 077c803231 e468e1eaa2"
commits="4c36c05105 0a3c58cfa7"
commits="81d1e91388 0a3c58cfa7"
commits="88a2acc0d9 1f9e96c46f"

#---------------------------------------------------------

ResetCommit
MakeCommitsBetween $commits
ShowCommits

echo "---------------------------" >> $overall
for commit in $commits
do
  DoRegress $commit
done
echo "---------------------------" >> $overall

#---------------------------------------------------------

