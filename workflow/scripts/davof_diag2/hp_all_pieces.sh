#!/bin/bash
# analyse every HPCToolkit database under prof/runs of the pieces profile (diag2, 2026-11-16): CPUTIME and, where present, REALTIME
W=/work/scratch/tm83tomy/davof/agents/diag2
P=$W/prof
R=$W/leia/workflow/scripts
mkdir -p $P/an
for mode in ${1:-serial par8 par64 s2d}; do
  for c in $P/runs/$mode/*/; do
    tag=$(basename $c); db=$c/hpctoolkit-database
    [ -f $db/profile.db ] || continue
    n=$(grep -c '^Time = ' $c/log.solver)
    for met in CPUTIME REALTIME; do
      out=$P/an/${mode}_${tag}.$met
      echo "== $mode/$tag steps $n $met"
      /usr/bin/python3 -I $R/davof_prof/hpan.py --db $db --out $out --metric $met --steps $n --rules $R/davof_diag2/rules_pieces.py --src 'leia|OpenFOAM' --tree-min 0.4 > $out.log 2>&1 || { echo "FAILED"; tail -3 $out.log; }
      head -4 $out.summary.txt 2>/dev/null | tail -2 | cut -c1-400
    done
  done
done
