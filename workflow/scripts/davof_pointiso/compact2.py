import re, sys
t = open(sys.argv[1]).read().split("\n")
cur = None
for l in t:
    if l.startswith("### Case"):
        cur = l.replace("OpenFOAM Co ", "Co ").split(":")[0]; print(cur)
        continue
    m = re.match(r"\| \((\d)\) \| (\d+) \| (.*)\|$", l)
    if m and "of" in l:
        cells = [c.strip() for c in m.group(3).split("|")]
        vals = [re.match(r"(\S+) \(", c).group(1) for c in cells[:4]]
        diff = cells[4]
        print("  R%s N%-4s thL %s thR %s clL %s clR %s | differ %s" % (m.group(1), m.group(2), vals[0], vals[1], vals[2], vals[3], diff))
        continue
    m = re.match(r"\| \((\d)\) \| ([-0-9.]+) \(([-0-9.]+)\) \| ([-0-9.]+) \(([-0-9.]+)\) \| ([-0-9.]+) \(([-0-9.]+)\) \| ([-0-9.]+) \(([-0-9.]+)\) \|$", l)
    if m:
        print("  orders R%s: thL %s (std %s) thR %s (%s) clL %s (%s) clR %s (%s)" % m.groups())
        continue
    m = re.match(r"\| (\d+) \| \((1)\)/\((\d)\) \| (.*)\|$", l)
    if m:
        print("  ratio N%s (1)/(%s): thL, thR, clL, clR = %s" % (m.group(1), m.group(3), m.group(4).strip()))
