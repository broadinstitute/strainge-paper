head -n 1 `/bin/ls -1 compare/*.summary.tsv | head -n 1` > compare.summary.txt
cat compare.summary.txt > compare.summary.chrom.txt
for FILE in `awk '{ print $4 }' compare.sh`; do tail -n +2 $FILE >> compare.summary.txt ; done
grep "	[45][0-9]\{6\}	" compare.summary.txt >> compare.summary.chrom.txt
