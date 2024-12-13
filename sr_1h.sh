#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --partition=norm2

  ntasks=32
  mpi=mpi
  lx=6
  ly=6
  hole=1h
  holeexist=2h
  ifphi0=0
  dim=500
  if [ ${mpi} == mpi ]
  then
    file=mpi_square_${lx}_${ly}_${hole}_${dim}dim_${holeexist}holeexist_$(date '+%y%m%d')_$$
  else
    file=square_${lx}_${ly}_${hole}_${dim}dim_${holeexist}holeexist_$(date '+%y%m%d')_$$
  fi


  mkdir ${file}
  cp sr_${hole} ${file}
  cp 1h_Parameter.txt ${file}
  if [ ${ifphi0} == 1 ]
  then
    cp ./wavef/lx${lx}ly${ly}_${dim}dim_${holeexist}/psi_file ${file}
    cp ./wavef/lx${lx}ly${ly}_${dim}dim_${holeexist}/sites_file ${file}
  fi
  #cp ./basis/lx${lx}ly${ly}_${dim}dim_1h_ps1_${holeexist}exist/basis.h5 ${file}
  cd ${file}
  
  #cat ${inputfile} > energy.m \
  #&& echo "ntasks=${ntasks}" >> energy.m \
  #&& sed -i 's/^/%/' energy.m \
  #&& csplit -s ${inputfile} 9 \
  #&& cat xx00 ./../src/minimize/RVB.inp xx01 >> tempt.txt \
  #&& cut -d '#' -f1 tempt.txt >liang_test_no_comment_$$.txt \
  #&& rm tempt.txt \
  #&&
  #if [ ${mpi} == mpi ]
  #then
  #  mpirun ./sr_${hole} liang_test_no_comment_$$.txt>>energy.m
  #else
  ./sr_${hole} 1h_Parameter.txt >>energy.txt
  #fi 


  rm sr_${hole}
  if [ ${ifphi0} == 1 ]
  then
    rm psi_file 
    rm sites_file
  fi
  cd ..
  #rm slurm*.out
  mv ${file} ./data/${hole}
