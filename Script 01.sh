for LIB in SRR25744243 SRR25744244 SRR25744245 SRR25744246 SRR25744247 SRR25744248 SRR25744249;
do

/home/imt/Desktop/Programas/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files $LIB

done
