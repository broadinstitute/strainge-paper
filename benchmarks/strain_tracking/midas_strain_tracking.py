from pathlib import Path


MIDAS_DIR = Path("all_vs_all/midas/")

COVS = ["0.1x", "0.5x", "1x", "10x"]
TRAINING_SAMPLES = []
for i in range(1, 11):
    for j in [1, 3, 7, 9]:
        TRAINING_SAMPLES.append(f"equal{i}.t{j}")

print(",".join(TRAINING_SAMPLES))
