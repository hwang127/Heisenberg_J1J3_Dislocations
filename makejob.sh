#!/bin/csh


foreach nsamples (20000)
foreach seed ( ` seq 1 3 ` )
foreach nwait(0)
foreach L( 300)
foreach num(1000 2000)
foreach T (0.01)
#(  0.416 0.417 0.418 0.419  0.42 0.421 0.422 0.423 0.424)
foreach alpha(4 )
foreach J(-1.0)
foreach Jthird(3.0)

# Input file
cat >/home-2/hwang127@jhu.edu/scratch/Dislo/GS/jobs/in_seed_${seed}_T_${T}_alpha_${alpha}_num_${num}_L_${L}.txt << EOFm
nsamples=${nsamples}
nwait=0
L=${L}
num=${num}
J=${J}
Jthird=${Jthird}
T=${T}
alpha=${alpha}

seed=${seed}
outfilename=/home-2/hwang127@jhu.edu/scratch/Dislo/GS/runs/seed_${seed}_T_${T}_alpha_${alpha}_num_${num}_L_${L}_

EOFm

cat >/home-2/hwang127@jhu.edu/scratch/Dislo/GS/jobs/job_seed_${seed}_T_${T}_alpha_${alpha}_num_${num}_L_${L}<<EOFm
#!/bin/bash -l
# Batch Queue Script
#SBATCH --time=72:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --mail-user=hwang127@jhu.edu
#SBATCH --cpus-per-task=1
#SBATCH --account=olegt
/home-2/hwang127@jhu.edu/scratch/Dislo/GS/mcfile /home-2/hwang127@jhu.edu/scratch/Dislo/GS/jobs/in_seed_${seed}_T_${T}_alpha_${alpha}_num_${num}_L_${L}.txt 

EOFm

chmod 755 /home-2/hwang127@jhu.edu/scratch/Dislo/GS/jobs/job_seed_${seed}_T_${T}_alpha_${alpha}_num_${num}_L_${L}
sbatch /home-2/hwang127@jhu.edu/scratch/Dislo/GS/jobs/job_seed_${seed}_T_${T}_alpha_${alpha}_num_${num}_L_${L}

end
end
end
end
end
end
end
end
end
