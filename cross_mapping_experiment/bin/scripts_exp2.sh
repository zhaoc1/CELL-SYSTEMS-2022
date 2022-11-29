cd /mnt/cz/simulation/GCF_000162155.1/5_exp2/paired
ls */*/temp/snps/*/*.bam | grep 92 |  xargs -Ixx bash -c 'samtools view xx | cut -f1 | cut -d'-' -f1 | sort | uniq -c | awk "$0"' '{print "xx", $0}' | grep cov_20X



cd /mnt/cz/simulation/GCF_000162155.1/2_reads/mngs
#grep "^@" cov_20_1.fastq | grep '/1' | cut -d'-' -f1 | sort | uniq -c
