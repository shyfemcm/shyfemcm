#!/bin/bash
#
# sets colors to be used in bash scripts
#
# example how to use:
#
# echo "... $name is ${green}installed${normal}"
#
#--------------------------------------------------------------

SetColors()
{
  local ncolors=$(tput colors)

  if [ -n "$ncolors" -a $ncolors -ge 8 ]; then
        bold="$(tput bold)"
        underline="$(tput smul)"
        standout="$(tput smso)"
        normal="$(tput sgr0)"
        black="$(tput setaf 0)"
        red="$(tput setaf 1)"
        green="$(tput setaf 2)"
        yellow="$(tput setaf 3)"
        blue="$(tput setaf 4)"
        magenta="$(tput setaf 5)"
        cyan="$(tput setaf 6)"
        white="$(tput setaf 7)"
  fi
}

#--------------------------------------------------------------

SetColors

#--------------------------------------------------------------

