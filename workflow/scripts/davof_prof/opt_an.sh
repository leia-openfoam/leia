#!/bin/bash
# the HPCToolkit analysis of the cycle TAG=$1 on the login node (the compute nodes' python lacks numpy)
P=/work/scratch/tm83tomy/davof/agents/swept/prof
T=$P/opt/$1/profile
[ -f $T/hpctoolkit-database/profile.db ] || { echo "$1: no database"; exit 0; }
n=$(grep -c '^Time = ' $T/log.advectorVoF)
RULES=${RULES:-$P/rules_swept2.py}
/usr/bin/python3 -I $P/hpan.py --db $T/hpctoolkit-database --out $T/an --steps $n --rules $RULES --src 'TwoPhaseFlow|OpenFOAM' --tree-min 0.4 > $T/an.log 2>&1 || { echo "$1: hpan failed"; tail -3 $T/an.log; exit 0; }
/usr/bin/python3 - $T/an.cat.csv $1 <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
tot = [r for r in rows if r["category"] == "(total)"][0]
print("%s: profile total %s s/step" % (sys.argv[2], tot["s_per_step_per_rank"]))
for r in sorted(rows, key=lambda r: -float(r["cputime_s"])):
    if r["category"].startswith("(") or float(r["pct"]) < 0.3: continue
    print("  %-78s %8s s/step %6s %%" % (r["category"][:78], r["s_per_step_per_rank"], r["pct"]))
EOF
grep -o 'faces by path: [^;]*\|pre-test: [^;]*\|trace cache: [^;]*\|cells: pieces [0-9]*' $T/log.advectorVoF | tail -4
