
#PBS -l nodes=1:ppn=36
#PBS -l walltime=23:59:00
#PBS -l pmem=3gb
#PBS -A lp_elect_gen_modeling group
#PBS -m abe
#PBS -M alexander.hoogsteyn@kuleuven.be

cd $VSC_DATA/EU-ETS-hydrogen

for i in {1..32}; do
julia --threads 9 MAIN.jl --start_scen $i --stop_scen $i --start_sens 1 --stop_sens 3 &
done
wait

