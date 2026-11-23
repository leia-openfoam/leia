import glob, collections, sys
R = "/work/scratch/tm83tomy/davof/agents/pointiso/leia/studies/davof/"
for tag in ("C02N1024", "C05N128", "B02N256"):
    for rowtag in ("r1", "r3"):
        f = glob.glob(R + "cl_dump_%s_%s/*_0????/postProcessing/contactLineError/dump/candidates.dat" % (rowtag, tag))[0]
        cand = collections.defaultdict(list)
        for l in open(f):
            w = l.split()
            if w[0] == "C":
                cand[(float(w[1]), int(w[2]))].append(dict(ch=int(w[3]), alpha=float(w[5]), mk=float(w[6]), x=float(w[10]), af=float(w[12])))
        for side in (0, 1):
            n = nm = nm6 = naf = 0
            for (t, s), lst in cand.items():
                if s != side:
                    continue
                ch = [c for c in lst if c["ch"]]
                if not ch:
                    continue
                c = ch[0]
                n += 1
                if c["mk"] < 1e-3:
                    nm += 1
                if c["mk"] == 0:
                    nm6 += 1
                if c["alpha"] < 1e-3 and c["af"] > 1e-3:
                    naf += 1
            print("%s %s side %s: steps %d; chosen with |m|/V^(2/3) < 1e-3: %d (exactly 0: %d); chosen with alpha_k < 1e-3 and wall alphaf > 1e-3: %d" % (tag, rowtag, "LR"[side], n, nm, nm6, naf))
