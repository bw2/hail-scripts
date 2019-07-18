import sys
import hashlib

def hash(s):
    return hashlib.sha224(s).hexdigest()[:12]

with open(sys.argv[1], "rt") as f, open(sys.argv[1].replace(".tsv", "") + ".obfuscated.tsv", "wt") as f2:
    is_header = True
    for line in f:
        if is_header:
            is_header=False
            f2.write(line)
            continue

        fields = line.split("\t")
        
        for i in [2,3,4]:
            fields[i] = hash(fields[i])

        f2.write("\t".join(fields))
