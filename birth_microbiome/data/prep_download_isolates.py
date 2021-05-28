import csv
from pathlib import Path

wget_cmd_tpl = "wget -N -c -q -P isolates/ ftp://{fname} 2> log/{log}.log"


with open("isolates_ena.txt") as f:
    reader = csv.DictReader(f, delimiter='\t')

    for row in reader:
        if "Enterococcus" not in row['scientific_name']:
            continue

        fastq = row['fastq_ftp']
        f1, f2 = fastq.split(';')

        print(wget_cmd_tpl.format(fname=f1, log=Path(f1).stem))
        print(wget_cmd_tpl.format(fname=f2, log=Path(f2).stem))
