#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# this script handles server specific settings
#
#--------------------------------------------------------------

debug="YES"
debug="NO"

#--------------------------------------------------------------

FindServer()
{
  host=$( hostname )

  if [ "$host" = galileo -o "$HPC_SYSTEM" = galileo ]; then
    server=galileo
  elif [ "$host" = fluxus -o "$SGE_CLUSTER_NAME" = fluxus ]; then
    server=fluxus
  else
    server=$host
  fi
}

#--------------------------------------------------------------

Usage()
{
  echo "Usage: check_server.sh [-h|-help] host"
}

FullUsage()
{
  Usage
}

ParseOptions()
{
  while [ -n "$1" ]
  do
    case $1 in
        -server)        what="server";;
        -check)         what="check";;
        -show)          what="show";;
        -reset)         what="reset";;
        -load)          what="load"; setting=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             echo "no such option: $1"; exit 1;;
        *)              break;;
    esac
    shift
  done

  compiler=$1
}

Sourced()
{
  if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "script is not sourced"
    echo "please run this script as \". ${0}\""
  fi
}

SetVars()
{
  if [ $server = "Caesium" ]; then
    export NETCDF_C_HOME=
    export NETCDF_FORTRAN_HOME=
    export NETCDFDIR=
    export METIS_HOME=$HOME/lib/metis
    export PETSC_HOME=
  else
    echo "unknown server: $server"
    exit 1
  fi
}

#--------------------------------------------------------------

Sourced

ParseOptions $*
FindServer

SetVars

#--------------------------------------------------------------

