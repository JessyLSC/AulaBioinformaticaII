#!/bin/bash

#Deve-se passar como primeiro argumento na linha de comando o diretório contendo os arquivos de entrada no formato .fastq 
#processados

input=$1

# Se naõ for passado diretório input, print então na tela a mensagem dizendo que está faltando o argumento, e se foi 
# passado mais é o diretório errado ou não contem o arquivo necessário, então print na tela dizendo que o argumento 
# passado é errado.

if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi

# output - É necessário passar como segundo argumento o diretório para armazenar o resultado do processo de montagem

output=$2

# Similarmente, se não foi passado nenhum argumento, print na tela que estão faltando, ou se foi passado e não é um
# diretório,print na tela que o argumento está incorreto

if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi

# Foi criado uma variável para representar o número de processadores a serem utilizados
num_threads="8"

# Foi criado outra variável para representar a quantidade de memória que se pode utilizar para o processo de montagem
mem_gb="10G"


# Arquivos e diretórios de saída (output) 
# Foram criadas varáveis para o diretório passada como outpur e os diretórios que serão criados durante o processo

basedir_out="${output}"

trinity_out="${basedir_out}/trinity_GG_assembled"

# Foram criados diretórios para as saídas dos programas que serão utilizados a seguir

mkdir -p ${trinity_out}

# Foi criado uma variável para os arquivos de alinhamento no formato bam

bamfiles=()

# Para a variável bamfiles, encontre dentro do diretório passado como input todos os arquivos com nome "Aligned.out.sorted.bam"

bamfiles=( $( find ${input} -name 'Aligned.out.sorted.bam' ) )

# Fazer uma fusão de todos os arquivos Aligned.out.sorted.bam encontrados

samtools merge -f ${basedir_out}/All.sorted.bam ${bamfiles[*]}

# Montagem dos transcritos utilizando o montador Trinity e o alinhamento (utilizando genoma de referência) como guia.
if [ ! -d ${trinity_out}/Trinity.timing ]; then
	
	echo -e "Assembling step (Trinity) ..."

	Trinity --KMER_SIZE 27 \
		--output ${trinity_out} \
		--seqType fq \
		--max_memory ${mem_gb} \
		--CPU ${num_threads} \
		--min_per_id_same_path 95 \
		--max_diffs_same_path  5 \
		--path_reinforcement_distance 5 \
		--group_pairs_distance 500 \
		--min_glue 5 \
		--min_contig_length 600 \
               	--min_kmer_cov 3 \
		--genome_guided_bam ${basedir_out}/All.sorted.bam \
		--genome_guided_max_intron 10000 \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
fi

# --KMER_SIZE: Tamanho dos kmers a serem originados (serão decompostos a partir das reads para serem então montados)
# --output: Diretório que conterá os dados de saída gerados.
# --seqType: Formato das reads que vão ser utilizadas no processo, no caso fastq.
# --max_memory: Quanto de espaço será destinado a esse processo.
# --CPU ${num_threads}: Número de processadores a serem utilizados.
# --min_per_id_same_path: minimo de indentidade em percentual para duas sequencias para  serem unidas em um só.
# --max_diffs_same_path: maximo de diferença permitida entre duas sequencias para combina-las
# --path_reinforcement_distance: minimo de sobreposição entre as reads.
# --group_pairs_distance: reads que não se encontrarem dentro dessa distancia, serão tratads como single-end
# --min_glue: número minimo de reads para que  cole ou junte os contigs gerados.
# --min_contig_length: tamanho minimo do contig
# --min_kmer_cov: minimo de kmers para serem montados na etapa Inchworm.
# --genome_guided_bam: arquivo de alinhamento a servir como guia (realizado na etapa acima)
# --genome_guided_max_intron: Tamanho maximo do intron do arquivo de alinhamento

