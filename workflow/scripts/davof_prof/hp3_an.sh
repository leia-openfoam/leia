#!/bin/bash
# the analysis of the round-3b profiles (prof/runs3, keys B.1-B.3) on the login node, then the
# side-by-side tables against the reference setting's profile (prof/an2) via hp_side.py
P=/work/scratch/tm83tomy/davof/agents/swept/prof
mkdir -p $P/an4
for mode in ${1:-serial par8 par64}; do
  for c in $P/runs3/$mode/*/; do
    tag=$(basename $c); db=$c/hpctoolkit-database
    [ -f $db/profile.db ] || continue
    n=$(grep -c '^Time = ' $c/log.advectorVoF)
    out=$P/an4/${mode}_$tag
    [ -f $out.summary.txt ] && [ $out.summary.txt -nt $db/profile.db ] && [ -z "$FORCE" ] && continue
    echo "== $mode/$tag steps $n"
    /usr/bin/python3 -I $P/hpan.py --db $db --out $out --steps $n --rules $P/rules_swept4.py --src 'TwoPhaseFlow|OpenFOAM' --tree-min 0.4 > $out.log 2>&1 || { echo "FAILED"; tail -5 $out.log; }
    head -4 $out.summary.txt | tail -2
  done
done
/usr/bin/python3 -I $P/hp_side.py $P /work/scratch/tm83tomy/davof/agents/swept/csv an2 an4 "407f969" "a11ebc1 + B.1-B.3" hp3
