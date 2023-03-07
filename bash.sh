

#PBS -l nodes=1:ppn=36
#PBS -l walltime=23:50:00
#PBS -l pmem=3gb
#PBS -A lp_elect_gen_modeling group
#PBS -m abe
#PBS -M alexander.hoogsteyn@kuleuven.be

cd $VSC_DATA/EU-ETS-hydrogen

julia MAIN.jl

