

for LIB in SRR16972309;
do

/home/imt/Desktop/Programas/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files $LIB

fastp -i fastq/${LIB}_1.fastq -I fastq/${LIB}_2.fastq -o fastq/${LIB}_1.trim.fastq -O fastq/${LIB}_2.trim.fastq #pré-processamento, avalia qualidade e realiza limpeza

bwa mem -t 4 ref/l_infantum/GCF_000002875.2_ASM287v2_genomic.fna fastq/${LIB}_1.trim.fastq fastq/${LIB}_2.trim.fastq | samtools sort --threads 4 -o map/${LIB}.bam #alinhamento contra o genoma de referência e posterior ordenamento

samtools depth map/${LIB}.bam > depth/${LIB}.txt #cálculo de profundidade


rm map/${LIB}.bam 

Rscript /home/imt/Desktop/leishmania/calculo_somia.R /home/imt/Desktop/leishmania/depth/${LIB}.txt >/home/imt/Desktop/leishmania/somia/${LIB}_somia.txt

echo "finalizado somia"

Rscript /home/imt/Desktop/leishmania/calculo_janela.R /home/imt/Desktop/leishmania/depth/${LIB}.txt >/home/imt/Desktop/leishmania/janela/${LIB}_janela.txt

echo "finalizado janela"

done

echo "todos os processos foram finalizados"
