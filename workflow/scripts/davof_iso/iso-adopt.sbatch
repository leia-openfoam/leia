#!/bin/bash
#SBATCH -J davof-iso-adopt
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=3600
#SBATCH -t 00:40:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The case stamps (workflow/scripts/case_stamp.py, the port agent's mechanism, merged
# 1478b631) adopted for the prepared studies of the isoAdvector verification: a snakemake
# dry run per study parses the study and stores the current stamp as the stamp of every
# finished case without one (adopted: true), then `case_stamp.py report`. Claude, 2026-11-07.
cd /work/scratch/tm83tomy/davof/agents/swept/iso/leia || exit 1
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. etc/leia-env.sh || exit 1
export PATH="$HOME/.local/bin:$PATH"
echo "leia $(git log --oneline -1 | cut -c1-80), $(date)"
for s in $(ls config/davof/iso*3D_*N*.yaml config/davof/isoShear2D_hexN*.yaml | sed 's|config/||; s|\.yaml$||'); do
  [ -d studies/$s ] || { echo "-- $s: no study directory"; continue; }
  echo "-- $s"; snakemake -n --quiet --workflow-profile profiles/slurm --nolock --rerun-triggers mtime --configfile config/$s.yaml > /dev/null 2> studies/$s/.adopt.err; echo "   dry run rc=$? ($(grep -c . studies/$s/.adopt.err) stderr lines)"
  python3 workflow/scripts/case_stamp.py report studies/$s 2>&1 | head -1 | cut -c1-200
done
c=studies/davof/isoRotation3D_hexN32/isoRotation3D_00000; echo "== a stored stamp"; ls -a $c | grep leia_stamp; python3 -c "import json; d=json.load(open('$c/.leia_stamp.json')); print({k: (v if not isinstance(v,(dict,list)) else type(v).__name__) for k, v in d.items()})" 2>&1 | cut -c1-600
echo "ISO-ADOPT-DONE $(date)"
