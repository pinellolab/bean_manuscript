import subprocess

variants_to_phase = [
    # "rs3767844_Maj_ABE_192",
    "rs62084210_Maj_ABE_250",
    # "rs116734477_Maj_ABE_51",
    "rs35081008_Min_ABE_465",
    # "rs4390169_Maj_ABE_208",
]

ps = []

for variant in variants_to_phase:
    p = subprocess.Popen(["python", "scripts/atac_seq/phase_alleles.py", variant])
    ps.append(p)

for p in ps:
    p.wait()
