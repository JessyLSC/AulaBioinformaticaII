#!/bin/bash

# É necessário passar um argumento que servirá como input,um diretório contendo os arquivos no formato .fastq
input=$1

# Para validação desse argumento, se não foi passado, printe na tela a mensagem que está faltando um diretório de entrada,
# Se não existe esse diretório, printe na tela um mensagem dizendo que está errado o argumento passado.

if [ ! ${input} ]
then   
        echo "Missing input (renamed for Trinity) directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input (renamed for Trinity) directory ${input}"
                exit
        fi
fi

# É necessário passar como segundo argumento o diretório contendo o arquivo de montagem

trinity_output=$2

# É realizada a validação desse argumento, "trinity_output", similarmente ao argumento anterior

if [ ! ${trinity_output} ]
then   
        echo "Missing Trinity output directory"
        exit
else   
        if [ ! -d ${trinity_output} ]
        then   
                echo "Wrong Trinity output directory ${trinity_output}"
                exit
        fi
fi

#É necessário passar como terceiro argumento o diretório para armazenar o resultado da avaliação de abundância
output=$3

# Novamente é feita a validação desse argumento
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

# Criando uma variável do número de processadores a serem utilizados
num_threads="8"

# Arquivos e diretórios de saída (output) 
# Criando a variável para o diretório de saída:

abundance_out="${output}/abundance"

# Criando o diretório de saída utilizando a função mkdir
mkdir -p ${abundance_out}

# Criando as variáveis para as reads1 (left) e reads2(right):

left=()
right=()

# Printa a mensagem entre aspas na tela:
echo "Collecting reads step ..."

# Para a variável left, encontrar todos os arquivos que finalizam o nome com .prinseq_1.fastq:

left=($(find ${input} -type f -name '*.prinseq_1.fastq'))

# Remoção dos arquivos gerados que são desnecessaríos:

rm -f ${abundance_out}/samples.txt
rm -f ${abundance_out}/quant_files.txt
rm -f ${abundance_out}/groups.txt

echo -e "id\tname\tgroup" > ${abundance_out}/groups.txt

# A estrutura for foi utilizada para fazer um loop da função abaixo para todas reads2 de uma só vez: Declarar a variável
# right: reads 2.
for l in ${left[@]}; do
	repname=`basename ${l} | sed 's/\..*$//'`
	condname=`echo ${repname} | sed 's/[0-9]\+//'`
	r=`echo ${l} | sed 's/_1.fastq/_2.fastq/'`
	right=(${right[@]} ${r})

	echo -e "${condname}\t${abundance_out}/${repname}\t${l}\t${r}" >> ${abundance_out}/samples.txt
	echo -e "${abundance_out}/${repname}/quant.sf" >> ${abundance_out}/quant_files.txt

	echo -e "${repname}\t${repname}\t${condname}" >> ${abundance_out}/groups.txt
done


#echo ${left[*]}
#echo ${right[*]}

# Criação de duas variáveis:trinity_fasta e trinity_trans_map: Buscar dentro do diretório de saída do programa Trinity,
# o arquivo de montagem (FASTA) para a primeira variável e o arquivo gene_trans_map para a segunda:

trinity_fasta=`find ${trinity_output} -type f -name 'Trinity*.fasta'`
trinity_trans_map=`find ${trinity_output} -type f -name Trinity*.gene_trans_map`

echo "Estimating abundances ..."

# Script pertencente ao programa Trinity para estimação da abundancia:

${TRINITY_HOME}/util/align_and_estimate_abundance.pl 	--transcripts	${trinity_fasta} \
							--est_method	salmon \
							--salmon_add_opts "--validateMappings" \
							--samples_file	${abundance_out}/samples.txt \
							--gene_trans_map ${trinity_trans_map} \
							--prep_reference \
							--thread_count ${num_threads} \
							--seqType fq \
							--output_dir ${abundance_out} \
							 > ${abundance_out}/align_and_estimate_abundance.log.out.txt \
							2> ${abundance_out}/align_and_estimate_abundance.log.err.txt

## --transcripts: Arquivo FASTA de montagem: variável trinity_fasta
## --est_method: Método de estimação da abundancia, o método salmon é rápido e realiza o processo em 2 etapas: indexação da
# referência e quantificação.
## --salmon_add_opts: Opções a serem adicionadas para o método salmon
## --samples_file: Arquivo contendo as reads, separado por condição
## --gene_trans_map: Arquivo gene_trans_map: variável trinity_trans_map
## --prep_reference: Constrói um index
## --thread_count: Número de processadores a serem utilizados
## --seqType: Formato do arquivo de sequencias: fastq
## --output_dir: Diretório que ira conter os resultados gerados

echo "Constructing abundance matrix ..."

#Construção da matriz de abundancia, utilizando um script também presente no programa Trinity:

${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl	--est_method salmon \
							--gene_trans_map ${trinity_trans_map} \
							--name_sample_by_basedir \
							--cross_sample_norm none \
							--quant_files ${abundance_out}/quant_files.txt \
							--out_prefix ${abundance_out}/abundance \
							 > ${abundance_out}/abundance_estimates_to_matrix.log.out.txt \
							2> ${abundance_out}/abundance_estimates_to_matrix.log.err.txt 

# --name_sample_by_basedir: Nomeia a coluna da amostra pelo nome do diretório e não pelo nome do arquivo
# --cross_sample_norm: Método de normalização
# --quant_files: Diretório contendo uma lista de todos os arquivos alvos.
# --out_prefix: Valores a serem utilizados para o método de estimação


echo "Calculating Differentially Expressed Genes ..."

# Criação do diretório para análise de genes diferencialmente expressos:

mkdir -p ${abundance_out}/DEG

# Análise de genes diferencialmente expressos utilizando o DESeq, utilizando como entrada a matriz gerada no passo anterior:

run-DESeq2.R 	--in="${abundance_out}/abundance.gene.counts.matrix"  \
		--groups="${abundance_out}/groups.txt" \
		--out="${abundance_out}/DEG" \
		 > ${abundance_out}/DEG/run-DESeq2.log.out.txt \
		2> ${abundance_out}/DEG/run-DESeq2.log.err.txt


