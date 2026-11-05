#!/bin/bash
# the analysis of the round-4 profiles (prof/runs4, the seven keys) on the login node, then the
# side-by-side tables against 407f969 (prof/an2) and a11ebc1 + B (prof/an4) via hp_side.py
P=/work/scratch/tm83tomy/davof/agents/swept/prof
mkdir -p $P/an5
for mode in ${1:-serial par8 par64}; do
  for c in $P/runs4/$mode/*/; do
    tag=$(basename $c); db=$c/hpctoolkit-database
    [ -f $db/profile.db ] || continue
    n=$(grep -c '^Time = ' $c/log.advectorVoF)
    out=$P/an5/${mode}_$tag
    [ -f $out.summary.txt ] && [ $out.summary.txt -nt $db/profile.db ] && [ -z "$FORCE" ] && continue
    echo "== $mode/$tag steps $n"
    /usr/bin/python3 -I $P/hpan.py --db $db --out $out --steps $n --rules $P/rules_swept5.py --src 'TwoPhaseFlow|OpenFOAM' --tree-min 0.4 > $out.log 2>&1 || { echo "FAILED"; tail -5 $out.log; }
    head -4 $out.summary.txt | tail -2
  done
done
C=/work/scratch/tm83tomy/davof/agents/swept/csv
echo "######## against 407f969"; /usr/bin/python3 -I $P/hp_side.py $P $C an2 an5 "407f969" "fbbfaf9 + B + gains 1-5" hp4
echo "######## against a11ebc1 + B"; /usr/bin/python3 -I $P/hp_side.py $P $C an4 an5 "a11ebc1 + B" "fbbfaf9 + B + gains 1-5" hp4b
