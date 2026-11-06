#!/usr/bin/env python3
import random

random.seed(7)
BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def rand(n):
    return "".join(random.choice(BASES) for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s, rate):
    out = []
    for b in s:
        r = random.random()
        if r < rate / 3:
            out.append(random.choice([c for c in BASES if c != b]))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(b)
            out.append(random.choice(BASES))
        else:
            out.append(b)
    return "".join(out)


segdup = rand(4000)
truth = {
    "chr1": rand(6000) + "CA" * 40 + rand(5000) + segdup + rand(7000) + "AAT" * 30 + rand(6000) + "GATA" * 25 + rand(5000),
    "chr2": rand(5000) + "AGG" * 35 + rand(6000) + mutate(segdup, 0.003) + rand(5000) + "TTAGGG" * 20 + rand(6000),
}

errors = {"chr1": [("CA" * 40, "CA" * 37), ("AAT" * 30, "AAT" * 32)], "chr2": [("AGG" * 35, "AGG" * 33)]}
ref = dict(truth)
for name, fixes in errors.items():
    for good, bad in fixes:
        ref[name] = ref[name].replace(good, bad, 1)


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


write_fasta("truth.fa", truth.items())
write_fasta("ref.fa", ref.items())

reads = []
for name, seq in truth.items():
    for i in range(len(seq) * 30 // 8000):
        n = random.randint(6000, 10000)
        s = random.randint(0, len(seq) - n)
        r = mutate(seq[s:s + n], 0.003)
        strand = "+"
        if random.random() < 0.5:
            r, strand = revcomp(r), "-"
        reads.append((f"read{len(reads) + 1} {name}:{s}-{s + n}{strand}", r))
write_fasta("read.fa", reads)
