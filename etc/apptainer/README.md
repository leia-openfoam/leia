# leia in an Apptainer image (a machine without a scheduler)

For a Linux machine that has Apptainer but no OpenFOAM build and no SLURM. The study
workflow is unchanged: only the profile differs (`profiles/container48`,
`profiles/container16`, the local executor with the `mpirun` of the image).

## What is where

- The image: `leia-openfoam-v2606.def` (Ubuntu 24.04, OpenFOAM-v2606 of openfoam.com,
  Open MPI, Snakemake 9.23.1 and the packages of `requirements.txt` in `/opt/leia-venv`).
  `/opt/leia-image.txt` in the image records the package versions.
- The code: the leia clone and the blastAMR clone on the host disk. Both are built INSIDE
  the image into their own `platforms/` folders, so a code change needs no new image.
- blastAMR: `$HOME/OpenFOAM/blastAMR-v2606` as seen inside the image (the default of
  `etc/leia-env.sh`), the fork at `9c4f667`, branch `v2606-build`.

## Procedure

1. Build the image on the machine (no root rights; `apptainer build` fails under WSL):

   ```bash
   apptainer build --fakeroot leia-openfoam-v2606.sif etc/apptainer/leia-openfoam-v2606.def
   ```

2. Make one work folder `W` that holds `leia/` (this clone), `home/OpenFOAM/blastAMR-v2606/`
   and the image. Start every command with the work folder bound to `/leia` and a core set:

   ```bash
   taskset -c <cores> apptainer exec --cleanenv --no-home -B W:/leia --pwd /leia \
       leia-openfoam-v2606.sif bash -c 'export HOME=/leia/home; <command>'
   ```

3. Build, inside the image:

   ```bash
   source /usr/lib/openfoam/openfoam2606/etc/bashrc
   (cd $HOME/OpenFOAM/blastAMR-v2606 && WM_PROJECT_USER_DIR=$PWD ./Allwmake)
   cd /leia/leia && . ./etc/leia-env.sh && ./Allwmake
   ```

4. Run a study, inside the image, with the profile whose core count equals the size of
   the core set of step 2:

   ```bash
   cd /leia/leia && make studies-one STUDY=<study> PROFILE=profiles/container48
   ```

## Rules

- The core count of the profile must equal the core set of `taskset`. A second study at
  the same time gets its own core set and its own profile.
- Start long runs in `tmux` (`exec` in the pane command, not `setsid`), and record the pid.
  Stop only your own pids.
- A result from this image and a result from Lichtenberg are compared with
  `workflow/scripts/compare_metrics_csv.py` and a tolerance: the compilers and the CPUs
  differ, so equality at tolerance 0 is not expected.
- No host name, user name or address of a machine goes into this repository.
