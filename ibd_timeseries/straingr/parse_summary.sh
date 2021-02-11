head -n 1 `/bin/ls -1 *.summary.tsv | head -n 1` | sed "s/^/sample\t/" > summary.txt
while read SAMPLE; do
    FILE=${SAMPLE}.summary.tsv
    tail -n +2 $FILE | sed "s/^/${SAMPLE}\t/" >> summary.txt
done < ../samples.txt

head -n 1 summary.txt > chrom.summary.txt
grep -P "\t[1-5][0-9]{6}\t" summary.txt | grep -vP "\tTOTAL\t" >> chrom.summary.txt
