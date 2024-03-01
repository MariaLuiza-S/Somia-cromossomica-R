

for LIB in SRR16972309;
do

	#/home/imt/Desktop/Programas/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump --split-files $LIB

	#pré-processamento, avalia qualidade e realiza limpeza
	fastp -i /home/imt/Desktop/leishmania/fastq/${LIB}_1.fastq \
	-I /home/imt/Desktop/leishmania/fastq/${LIB}_2.fastq \
	-o /home/imt/Desktop/leishmania/fastq/${LIB}_1.trim.fastq \
	-O /home/imt/Desktop/leishmania/fastq/${LIB}_2.trim.fastq 

	#alinhamento contra o genoma de referência e posterior ordenamento
	bwa mem -t 4 \
	/home/imt/Desktop/leishmania/ref/l_infantum/GCF_000002875.2_ASM287v2_genomic.fna \
	/home/imt/Desktop/leishmania/fastq/${LIB}_1.trim.fastq \
	/home/imt/Desktop/leishmania/fastq/${LIB}_2.trim.fastq | \
	samtools sort --threads 4 \
	-o /home/imt/Desktop/leishmania/map/${LIB}.bam 
	
	#cálculo de profundidade
	samtools depth /home/imt/Desktop/leishmania/map/${LIB}.bam \
	> /home/imt/Desktop/leishmania/depth/${LIB}.txt


	rm /home/imt/Desktop/leishmania/map/${LIB}.bam 

	Rscript /home/imt/Desktop/leishmania/calculo_somia.R /home/imt/Desktop/leishmania/depth/${LIB}.txt \
	>/home/imt/Desktop/leishmania/somia/${LIB}_somia.txt

echo "finalizado somia ${LIB}"

	Rscript /home/imt/Desktop/leishmania/calculo_janela.R /home/imt/Desktop/leishmania/depth/${LIB}.txt \
	>/home/imt/Desktop/leishmania/janela/${LIB}_janela.txt

echo "finalizado janela ${LIB}"

done

echo "todos os processos foram finalizados"
