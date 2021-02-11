qacct -j 2328095 > qacct.txt; echo "runtime	memory" > time_mem.txt; grep -e ru_utime -e maxvmem qacct.txt | awk '{ print $2 }' | paste - - >> time_mem.txt
