#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002875.2/
#09 de janeiro de 2024
#bwa index GCF_000002875.2_ASM287v2_genomic.fna #indexar o genoma de referência

for LIB in SRR25744247 SRR25744248 SRR25744249;
do

fastp -i fastq/${LIB}_1.fastq -I fastq/${LIB}_2.fastq -o fastq/${LIB}_1.trim.fastq -O fastq/${LIB}_2.trim.fastq #pré-processamento, avalia qualidade e realiza limpeza



bwa mem -t 4 ref/l_infantum/GCF_000002875.2_ASM287v2_genomic.fna fastq/${LIB}_1.trim.fastq fastq/${LIB}_2.trim.fastq | samtools sort --threads 4 -o map/${LIB}.bam #alinhamento contra o genoma de referência e posterior ordenamento

samtools depth map/${LIB}.bam > depth/${LIB}.txt #cálculo de profundidade


rm map/${LIB}.bam 

done

#samtools view -bS --threads 4 map/SRR25744243.sam > map/SRR25744243.bam

#samtools sort --threads 4 map/SRR25744243.bam -o map/SRR25744243_sorted.bam
#fastq/${LIB}_2.fastq
#fastq/$LIB_2.fastq
