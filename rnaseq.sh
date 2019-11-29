#!/bin/bash

# Número de processadores a ser utilizado durante o processo. Criando a variável, basta substituir nos comandos.

num_threads=8

# É necessário passar um diretório input de onde vai estar contido os dados de entrada para o alinhamento, na linha de
# comando ao chamar o script para ser executado.

indir=$1

# SE(if) não (!) for passado um diretório de entrada, print na tela a mensagem "Missing input directory."
if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi
# A função echo, mostra na tela a mensagem que etá entre aspas.

# SE for passado um argumento na linha de comando, mais esse não for o diretório, print na tela a mensagem:"Wrong input directory (${indir})."
if [ ! -d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

# É necessário passar um diretório de saída, ou seja, que vai conter os dados gerados.

outdir=$2

# SE  não for passado o argumento 2, o diretório de saída na linha de comando, print a mensagem:"Missing output directory."
if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE  o argumento passado não é um diretório, print na tela a mensagem dizendo:"Wrong output directory (${outdir})."

if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

# É necessário passar um terceiro argumento na linha de comando, o arquivo de anotação do genoma no formato gtf.

refgtf=$3

# SE ${refgtf} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO
if [ ! ${refgtf} ]; then
	echo "Missing GTF file."
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "Not found GTF file (${refgtf})."
	exit
fi

# É necessário passar um quarto argumento na linha de comando: Arquivo FASTA referente ao genoma a ser utilizado como referência.

refseq=$4

# SE ${refseq} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO
if [ ! ${refseq} ]; then
	echo "Missing GENOME fasta file."
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})."
	exit
fi

# Primeiramente é chamado para a execução o script preprocess.sh, que realizará o processamento dos dados e controle de
# qualidade. O mesmo necessita de 2 argumentos na linha de comando: O diretório de entrada de dados (Arquivos fastq) e
# um diretório para saída dos dados processados, que por sua vez serão o input do presente script.

./preprocess.sh "${indir}" "${outdir}"

# O comando mkdir cria diretórios. Assim os diretórios abaixo serão criados dentro do diretório passado como saída de dados
# na linha de comando ao chamar o script.

mkdir -p ${outdir}/star_index
mkdir -p ${outdir}/star_out_pe
mkdir -p ${outdir}/star_out_se
mkdir -p ${outdir}/star_out_final
mkdir -p ${outdir}/cufflinks
mkdir -p ${outdir}/cuffmerge
mkdir -p ${outdir}/stringtie
mkdir -p ${outdir}/stringmerge
mkdir -p ${outdir}/cuffcompare
mkdir -p ${outdir}/cuffquant

# A estrutura "for" serve como um loop, facilitando o script, pois assim não é necessário fazer para cada amostra
# separadamente.  A estrutura abaixo diz: para todos os dados contidos no diretório prinseq que finalizam com
# .atropos_final.prinseq_1.fastq, faça o alinhamento.
# "ls" é pra procurar dentro do diretório todos os arquivos que terminam com .atropos_final.prinseq_1.fastq.
# a função "sed" é utilizada para substituição.
# a função "touch" cria arquivos vazios.

for r1 in `ls ${outdir}/processed/prinseq/*.atropos_final.prinseq_1.fastq`; do
	r1_singletons=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_1_singletons.fastq/'`
	if [ ! -e "${r1_singletons}" ]; then
		touch ${r1_singletons}
	fi

	r2=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_2.fastq/'`

	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi

	r2_singletons=`echo ${r2} | sed 's/prinseq_2.fastq/prinseq_2_singletons.fastq/'`
	if [ ! -e "${r2_singletons}" ]; then
		touch ${r2_singletons}
	fi

	name=`basename ${r1} | sed 's/.atropos_final.prinseq_1.fastq//'`

# Se não exister o index do genoma, então faça o index.
# Para a construção do index a ferramenta STAR será utilizada.
# O index facilitará o mapeamento das reads contra o genoma de referência.

	if [ ! -e "${outdir}/star_index/SAindex" ]; then
		echo "Indexing genome (${refseq}) ..."


		STAR 	--runThreadN        ${num_threads} \
			--runMode           genomeGenerate \
			--genomeFastaFiles  ${refseq} \
			--genomeDir         ${outdir}/star_index \
			--sjdbGTFfile       ${refgtf} \
#			--genomeSAindexNbases 12 \ Apenas para genomas muito pequenos (escalonamento)
			--sjdbOverhang      149 \
		 > ${outdir}/star_index/STAR.index.log.out.txt \
		2> ${outdir}/star_index/STAR.index.log.err.txt

	fi

# --runThreadN: Número de processadores a ser utilizado, foi criado uma variável no inicio do script.
# --runMode: Direciona o STAR a construir o index do genoma.
# --genomeFastaFiles: Arquivo FASTA referente ao genoma de referência.
# --genomeDir: Diretório de saída que irá conter o index gerado.
# --sjdbGTFfile: Arquivo gtf - anotação do genoma de referência.
# --sjdbOverhang: Especifica o tamanho da sequencia genomica levando em consideração a anotação para ser usado na construção
# do banco de dados de junções de splicing. Esse tamanho geralmente é igual ao tamanho das reads -1. Logo, 149.
#  > log.out.txt: saída que contém as etapas da criação do index.
# 2> log.err.txt: segunda saída que contem erros que venham a ocorrer durante o processo.

	echo "STAR alignment PE with sample ${name}: ${r1} & ${r2} ..."

	mkdir -p ${outdir}/star_out_pe/${name}

# Para o alinhamento das reads contra o genoma de referência também será utlizado a ferramenta STAR.

	STAR	--runThreadN        ${num_threads} \
		--genomeDir         ${outdir}/star_index \
		--readFilesIn       ${r1} ${r2} \
		--outSAMstrandField intronMotif \
		--outFilterIntronMotifs RemoveNoncanonical \
		--sjdbGTFfile       ${refgtf} \
		--outFilterMultimapNmax 20 \
		--outFileNamePrefix ${outdir}/star_out_pe/${name}/ \
		--outSAMtype        BAM Unsorted \
		 > ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.out.txt \
		2> ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.err.txt

# --readFilesIn: Reads (no caso já processadas) a serem utilizadas no mapeamento.
# --outSAMstrandField: Para dados unstranded,as ferramentas Cufflinks/Cuffdiff requerem alinhamento de splicing com o
# o atributo de fita XS, que será gerado com esse parâmetro e a opção intronMotif.
# --outFilterIntronMotifs: Remove introns não canônicos.
# --outFilterMultimapNmax: Número máximo de alinhamentos multiplos permitido por read, se exceder esse valor, a read é
# considerada não mapeada; Foi utilizado o valor default.
# --outFileNamePrefix: prefixo a ser dado para o arquivo gerado.
# --outSAMtype: arquivo de saída do tipo BAM sem estar ordenado por nome ou por coordenada.

	echo "STAR alignment SE with sample ${name}: ${r1_singletons} & ${r2_singletons} ..."

	mkdir -p ${outdir}/star_out_se/${name}

# Para as reads que não foram encontradas seu par, denomindas singletons, realizar o alinhamento:

	STAR	--runThreadN        ${num_threads} \
		--genomeDir         ${outdir}/star_index \
		--readFilesIn       ${r1_singletons},${r2_singletons} \
		--sjdbGTFfile       ${refgtf} \
		--outSAMtype        BAM Unsorted \
		--outFilterMultimapNmax 20 \
		--outSAMstrandField intronMotif \
		--outFileNamePrefix ./$outdir/star_out_se/${name}/ \
		 > ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.out.txt \
		2> ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.err.txt

	echo "Merging STAR alignment PE & SE (${name}) ..."

# STAR_out é a junção do alinhamento das reads paired-end e das singletons

	mkdir -p ${outdir}/star_out_final/${name}

# Combinar resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons):
# Utilizar a função samtools merge para junção de ambos os alinhamentos.
# -f: caso tenha o arquivo, sobreescreva.
# -n: ordena o alinhamento por nome.
# @: número de processadores.

	 samtools merge -@ ${num_threads} -f -n  ${outdir}/star_out_final/${name}/Aligned.out.bam \
                                                ${outdir}/star_out_pe/${name}/Aligned.out.bam \
                                                ${outdir}/star_out_se/${name}/Aligned.out.bam \
	 > ${outdir}/star_out_final/${name}/samtools.merge.log.out.txt \
	2> ${outdir}/star_out_final/${name}/samtools.merge.log.err.txt

	echo "Sorting STAR alignment final (${name}) ..."

# Ordenar o resultado do alinhamento por coordenadas genômicas, pois é uma exigência para executar o cufflinks posteriormente.
# Para a ordenação utilzar a função samtools sort:

	 samtools sort -@ ${num_threads} -o      ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
                                                ${outdir}/star_out_final/${name}/Aligned.out.bam \
	 > ${outdir}/star_out_final/${name}/samtools.sort.log.out.txt \
	2> ${outdir}/star_out_final/${name}/samtools.sort.log.err.txt

	echo "Collecting alignment statistics (${name}) ..."

# Gerando as estatisticas do alinhamento utilizando o script SAM_nameSorted_to_uniq_count_stats.pl.
# É necessário passar  o arquivo de entrada (BAM final, contendo todos os alinhamentos) e o nome para o arquivo de saída.

	SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/star_out_final/${name}/Aligned.out.bam > ${outdir}/star_out_final/${name}/Aligned.stats.txt

# Após o mapeamento das reads contra o genoma de referência, o transcriptoma será então montado utilizando o cufflinks:

	echo "Running Cufflinks (${name}) ..."

	mkdir -p ${outdir}/cufflinks/${name}

	cufflinks --output-dir ${outdir}/cufflinks/${name} \
		  --num-threads ${num_threads} \
		  --GTF-guide ${refgtf} \
		  --frag-bias-correct ${refseq} \
		  --multi-read-correct \
		  --library-type fr-unstranded \
		  --frag-len-mean 250 \
		  --frag-len-std-dev 40 \
		  --total-hits-norm \
		  --compatible-hits-norm \
		  --max-frag-multihits 20 \
		  --min-isoform-fraction 0.25 \
		  --max-intron-length 9000 \
		  --min-intron-length 80 \
		  --overhang-tolerance 10 \
		  --max-bundle-frags 999999 \
		  --max-multiread-fraction 0.65 \
		  --overlap-radius 20 \
		  --3-overhang-tolerance 300 \
		  ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
		 > ${outdir}/star_out_final/${name}/cufflinks.log.out.txt \
		2> ${outdir}/star_out_final/${name}/cufflinks.log.err.txt

# --output-dir: Diretório para a saída de dados.
# --num-threads: Número de processadores a ser utilizado.
# --GTF-guide: Arquivo de anotação referente ao genoma de referência no formato gtf.
# --frag-bias-correct: Arquivo FASTA do genoma de referência.
# --multi-read-correct: Esse parametro solicita uma estimaçao inicial das reads mapeadas em multiplas posições.
# --library-type: Indica o tipo de biblioteca.
# --frag-len-mean: Média do tamanho dos fragmentos.
# --frag-len-std-dev: Desvio padrão referente ao tamnho dos fragmentos.
# --overhang-tolerance: Número de pb permitidos a entrar em um intron, caso seja mapeado corretamente.
# --total-hits-norm: Conta todos os fragmentos, incluindo aqueles não compativeis com nenhum transcrito de referência,
# levando ao número de mapeamentos, utilizando o FPKM.
# –compatible-hits-norm: Conta somente os fragmentos compativeis com algums transcrito de referência, levando ao número
# de mapeamento, utilizando o FPKM.
# --min-isoform-fraction: Após calcular a abundancia de isoformas para um gene, Cufflinks elimina transcritos que apresentem
# baixa abundancia, pois esses transcritos podem não ser montados de maneira confiável, podendo ser artefatos.Também elimina
# introns que apresentem poucos alinhamentos de splicing.
# --max-intron-length: Tamanho máximo dos introns.
# --min-intron-length: Tamanho minimo dos introns.
# --max-bundle-frags: Número maximo de fragmentos que um locus deve apresentar para que então passe para o proximo.
# --max-multiread-fraction: Porção de reads que podem ser mapeadas em multiplas posições.
# --overlap-radius: Transcrito que estão separados por menos que essa distancia são unidos e o gap é preenchido.
# --3-overhang-tolerance: Número de pares de bases que se excedem na porção 3', para ser permitido que se una ao transcrito ja montado.
# Será também realizada uma montagem utilizando o software Stringtie:

	echo "Running StringTie (${name}) ..."

	mkdir -p ${outdir}/stringtie/${name}

# Deve-se passar primeiramente o arquivo de alinhamento

	stringtie ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
		-G ${refgtf} \
		-o ${outdir}/stringtie/${name}/transcripts.gtf \
		-p ${num_threads} \
		-f 0.25 \
		-a 10 \
		-j 5 \
		-c 6 \
		-g 20 \
		-M 0.65 \
		-A ${outdir}/stringtie/${name}/gene_abundance.txt

# -G: Arquivo de anotação referente ao genoma de referência no formato gtf.
# -o: Arquivo de saída - formato gtf.
# -p: Número de processadores a serem utilizados no processo.
# -f: Fração minima de abundância de isoformas.
# -a: Junções que não apresentam reads que alinharam através delas em pelo menos essa quantidade de bases em ambos os lados são filtradas. 
# -c: Minimo de cobertura permitido para os transcritos preditos.
# -g: Transcrito que estão separados por menos que essa distancia são unidos e o gap é preenchido.
# -M: Porção de reads que podem ser mapeadas em multiplas posições.
# -A: Saida do arquivo de abudancia dos genes - formato txt.

done

# Junção de todos os arquivos de saída "transcripts.gtf" do cufflinks em um só - assembly_GTF_list.txt, utilizando a 
# ferramenta cuffmerge.

echo "Running cuffmerge ..."

# O comando find encontra todos os arquivos com esse nome no diretório dado.

find ${outdir}/cufflinks/ -name 'transcripts.gtf' > ${outdir}/cuffmerge/assembly_GTF_list.txt

cuffmerge -o ${outdir}/cuffmerge/ \
	--ref-gtf ${refgtf} \
	--ref-sequence ${refseq} \
	--min-isoform-fraction 0.20 \
	--num-threads ${num_threads} \
	   ${outdir}/cuffmerge/assembly_GTF_list.txt \
	 > ${outdir}/cuffmerge/cuffmerge.log.out.txt \
	2> ${outdir}/cuffmerge/cuffmerge.log.err.txt

# -o: Saída do arquivo gerado.
# -ref-gtf: Arquivo de anotação referente ao genoma no formato gtf.
# - ref-sequence: Arquivo no formato FASTA do genoma de referência.
# -min-isoform-fraction: minimo de abundancia de isoformas, abaixo disso são eliminadas.
# -num-threads: Número de processadores a serem utilizados.


# Similarmente ao processo anterior, se realizará a junção dos arquivos gerados do StringTie em um só arquivo, utilizando a
# ferramenta stringtie merge:

echo "Running stringtie merge ..."

find ${outdir}/stringtie/ -name 'transcripts.gtf' > ${outdir}/stringmerge/assembly_GTF_list.txt

stringtie --merge \
	-G ${refgtf} \
	-o ${outdir}/stringmerge/merged.gtf \
	-c 1 \
	-T 1 \
	-f 0.25 \
	-g 20 \
	-i \
	${outdir}/stringmerge/assembly_GTF_list.txt

# -T: minimo de transcritos TPM a serem utilizados na junção.
# -i: mantém os transcritos unidos com retenção dos introns.


# Em seguida será realizada uma comparação das montagens obtida com a referência utilizando a ferramenta cuffcompare:

cuffcompare	-r ${refgtf} \
		-s ${refseq} \
		-o ${outdir}/cuffcompare/stringmerge \
		${outdir}/stringmerge/merged.gtf \
		 > ${outdir}/stringmerge/cuffcompare.log.out.txt \
		2> ${outdir}/stringmerge/cuffcompare.log.err.txt

# Criação da variávek biogroup_label:

biogroup_label=()

for bamfile in `ls ${outdir}/star_out_final/*/Aligned.out.sorted.bam`; do
	name=`basename $(dirname ${bamfile})`

# Quantificação dos transcritos utilizando a ferramenta cuffquant:

	echo "Running cuffquant using sample ${name} with ${outdir}/stringmerge/merged.gtf as reference ..."

	mkdir -p ${outdir}/cuffquant/${name}

	cuffquant 	--output-dir ${outdir}/cuffquant/${name} \
			--frag-bias-correct ${refseq} \
			--multi-read-correct \
			--num-threads ${num_threads} \
			--library-type fr-unstranded \
			--frag-len-mean 250 \
			--frag-len-std-dev 40 \
			--max-bundle-frags 9999999 \
			--max-frag-multihits 20 \
			${outdir}/stringmerge/merged.gtf \
			${bamfile} \
		 > ${outdir}/cuffquant/${name}/cuffquant.log.out.txt \
		2> ${outdir}/cuffquant/${name}/cuffquant.log.err.txt

# Os parâmetros são similares ao cufflinks.

	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done

biogroup_files=()

# Para realizar a análise de genes diferencialmente expressos,utilizaremos os arquivo de abundância gerados na etapa anterior.

echo "Running Differential Expression Analysis ..."

# Será utilizado novamente a estrutura for, para fazer com todas as amostras de uma só vez.
# Preparo dos dados para cuffdiff

for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()
	for cxbfile in `ls ${outdir}/cuffquant/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done
	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}
echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

# Normalização dos resultados utilizando o cuffnorm:

mkdir -p ${outdir}/cuffnorm/

cuffnorm 	--output-dir ${outdir}/cuffnorm \
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		--num-threads ${num_threads} \
		--library-type fr-unstranded \
 		--library-norm-method geometric \
		--output-format simple-table \
 		${outdir}/stringmerge/merged.gtf \
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffnorm/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm/cuffdiff.log.err.txt

# -output-format: formato da saída: Tabela simples
# Restante dos parametros já foram citados no cufflinks.

echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

# Análise de genes diferencialmente expressos utilizando o cuffdiff

mkdir -p ${outdir}/cuffdiff/ 

cuffdiff 	--output-dir ${outdir}/cuffdiff \
 		--labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		--frag-bias-correct ${refseq} \
 		--multi-read-correct \
 		--num-threads ${num_threads} \
 		--library-type fr-unstranded \
 		--frag-len-mean 250 \
 		--frag-len-std-dev 40 \
 		--max-bundle-frags 9999999 \
 		--max-frag-multihits 20 \
 		--total-hits-norm \
 		--min-reps-for-js-test 2 \
 		--library-norm-method geometric \
 		--dispersion-method per-condition \
 		--min-alignment-count 12 \
 		${outdir}/stringmerge/merged.gtf \
 		${biogroup_files[*]} \
 	 	> ${outdir}/cuffdiff/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffdiff/cuffdiff.log.err.txt

# --min-reps-for-js-test: Não será testada a expressão genica diferencial se não houver esse número minimo de réplicas.
# --library-norm-method: Normalização do método. Geométrica: FPKMs e contagem dos fragmentos são escalados com base na
# mediana das médias geométricas da contagem de fragmentos em todas as bibliotecas.
# --dispersion-method: Método de dispersão. Por condição: Cada réplica recebe seu próprio modelo, não é realizado um pool.
# --min-alignment-count:Número minimo de alinhamentos em um loco para esse ser submetido ao teste.
# Restante dos parametros já foram citados anteriormente
