# configurando o base dir
BASEDIR=${PWD}

# Criando as pastas necessárias
#mkdir -p ${BASEDIR}/fastq
#mkdir -p ${BASEDIR}/depth
#mkdir -p ${BASEDIR}/map
#mkdir -p ${BASEDIR}/somia
#mkdir -p ${BASEDIR}/janela

# Definindo a LIB que será processada
for LIB in SRR20748434 ERR2300714 ERR966775 SRR12506675; do

    # Executando o comando 
    docker run -v ${BASEDIR}:/mnt r_ploidy:latest Rscript /mnt/calculo_janela2.R /mnt/depth/${LIB}.txt \
    > ${BASEDIR}/janela/${LIB}_janela2.txt

    echo "finalizado janela ${LIB}"

done

