#!/bin/bash
#SBATCH -J leia-amr-u1
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=3000
#SBATCH -t 00:30:00
#SBATCH -o int.u1.%j.out
#SBATCH -e int.u1.%j.err
# Unit gate U1: leiaTestMeshTransfer (the mesh transfer of psi on a blastAMR mesh), in
# serial and on 4 ranks, for the four test fields. Submit from the clone root:
#   sbatch --export=ALL,INT_SCRIPTS=<dir>,SERIAL=<prepared serial case>,PAR=<prepared case with 4 subdomains> transfer_gate.sbatch
# The cases are COPIED (time 0 and constant only); the study directories are not changed.
# The numbers are in int.u1.log; the last action is int.u1.done.
cd "$SLURM_SUBMIT_DIR" || exit 1
. "${INT_SCRIPTS:?}/env.sh"
rm -f int.u1.done
{
( cd applications/test/leiaTestMeshTransfer && wmake ) > int.u1.build.log 2>&1; echo "wmake rc=$?"; tail -3 int.u1.build.log
ls -la "$FOAM_USER_APPBIN/leiaTestMeshTransfer"
W=studies/slAmr/u1TransferGate; rm -rf $W; mkdir -p $W/serial $W/np4
cp -r $SERIAL/0 $SERIAL/constant $SERIAL/system $W/serial/
cp -r $PAR/constant $PAR/system $W/np4/
for p in $PAR/processor*; do mkdir -p $W/np4/$(basename $p); cp -r $p/0 $p/constant $W/np4/$(basename $p)/; done
for f in linear harmonic quadratic sphere; do
  ( cd $W/serial && leiaTestMeshTransfer -field $f > log.u1.$f 2>&1; echo "== serial $f rc=$?"; grep -E "^pass|leiaTestMeshTransfer field|PASS|FAIL|REPORT|FOAM FATAL" -A0 log.u1.$f | tail -18 )
  ( cd $W/np4 && mpirun -np 4 leiaTestMeshTransfer -parallel -field $f > log.u1.$f 2>&1; echo "== np4 $f rc=$?"; grep -c "Successfully distributed" log.u1.$f; grep -E "^pass|leiaTestMeshTransfer field|PASS|FAIL|REPORT|FOAM FATAL" log.u1.$f | tail -18 )
done
} > int.u1.log 2>&1
echo U1_FINISHED > int.u1.done
