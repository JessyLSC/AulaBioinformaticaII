#!/bin/bash


# É necessário passar um diretório  de entrada de dados ao chamar o script para execução:

indir=$1


# ! Significa negação, no caso abaixo, se não foi passado diretório de entrada print a mensagem "Missing input directory"

if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi

# SE  o que foi passado não é um diretório, print a mensagem: "Wrong input directory"

if [ !-d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

# Também é necessário passar um diretório de saída, ou seja onde os arquivos gerados ficarão:

outdir=$2

# SE não foi passado o diretório de saída, print uma mensagem dizendo que o argumento está faltando:

if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE o que foi passado não é um diretório, print uma mensagem dizendo que argumento está errado:
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

# A função mkdir serve para criação de diretórios,no caso dentro do diretório que foi passado como output será criada uma
# pasta chamada processed, e dentro da mesma, serão criadas outras 3 pastas: atropos, prinseq e fastqc.
# A última conterá outras duas pastas: pre e pos. Ou seja, antes do processamento as reads passarão por um controle de
# qualide(pre) e após o processamento, serão novamente submetidas ao controle de qualidade (pos), pondendo realizar uma
# comparação de ambos posteriormente:

mkdir -p ${outdir}/processed/fastqc/pre
mkdir -p ${outdir}/processed/atropos
mkdir -p ${outdir}/processed/prinseq
mkdir -p ${outdir}/processed/fastqc/pos

# A estrutura for é utilizada para facilitar o processo, com ela pode-se realizar o comando para todas as reads,
# não sendo necessário a execução para cada uma delas. Assim, no comando abaixo, para cada read no formato .fastq será
# executado o mesmo comando: Pegar a base de cada nome (o que difere) e realizar o controle de qualidade (comando fastqc),
# gerando 2 arquivos um de erro (log.err.txt), caso venha a ocorrer, e outro de saída que mostra o andamento do processo
# (log.out.txt).

for r1 in `ls ${indir}/*_R1.fastq`; do

	r2=`echo ${r1} | sed 's/_R1.fastq/_R2.fastq/'`

	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi

	name=`basename ${r1} | sed 's/_R1.fastq//'`
	echo -e "FastQC pre-evaluation using sample ${name}: ${r1} & ${r2} ...\n"

	fastqc -t 2 \
   		${r1} \
   		-o ${outdir}/processed/fastqc/pre/ \
		 > ${outdir}/processed/fastqc/pre/${name}_R1.log.out.txt \
		2> ${outdir}/processed/fastqc/pre/${name}_R1.log.err.txt

	fastqc -t 2 \
   		${r2} \
		-o ${outdir}/processed/fastqc/pre/ \
		 > ${outdir}/processed/fastqc/pre/${name}_R2.log.out.txt \
		2> ${outdir}/processed/fastqc/pre/${name}_R2.log.err.txt

# A função echo é para printar na tela a mensagem que vem em seguida entre aspas. Isso contribui para saber em que etapa 
# o processo se encontra.

	echo -e "Running atropos (insert) for adapter trimming using sample ${name}: ${r1} & ${r2} ...\n"

# É necessário inicializar o ambiente pyenv dentro desta sessão, pois para a execução do atropos é necessário um 
# ambiente Python, dado o seu desenvolvimento

        eval "$(pyenv init -)"

# Para ativar o ambiente Python para execução do atropos:
	pyenv activate atropos

# Execução do Atropos - trimagem das reads, eliminando os adaptadores e sequencias de baixa qualidade. O atropos funciona
# utilizando 2 algortimos de trimagem o "insert" e o "adapter". O insert realiza a trimagem após alinhamento das reads
# por match das reads, o que sobrar é adaptadore e então ocorre a trimagem. e as reads que não foram trimadas nessa etapa
# passa para o segundo algortimo o adapter, que realiza a trimagem baseando-se na correspondendia do alinhamento dos adaptadores. 

# -e: Esse parâmetro fornece o valor máximo da taxa de erro permitida para matches de adaptadores para que o alinhamento 
# seja considerado que obteve sucesso. É o número de erros dividido pelo comprimento da região de matching.

# -n: Nessa opção deve-se optar por até quantos adaptadores deve ser removidos de cada read. 
# No caso do nosso comando, o parâmtro diz: remova até dois adaptadores de cada read.

# -m: Essa opção descarta reads trimadas que apresentem tamanho menor que X (valor que você opta).
#  No comando abaixo, reads menores que 15 pb serão descartadas.

# Op-order: Esse parâmetro descreve a ordem que as operações de trimagem são aplicadas. As letras significam:
# -A: Trimagem do adaptador,C: Corte,G: Trimagem NextSeq,Q: Qualidade da trimagem,W: Sobrescrever leituras de baixa qualidade.
# O default é WCGQA entretanto no comando acima a ordem foi alterada para GAWCQ devido a experiências prévias.

# Match-read-wildcards: Modo de interpretação da IUPAC está habilitado, assim ele identifica caracteres coringas quando
# não se sabe qual base será inserida por exemplo.

# -O: Nesse parâmetro você fornece um valor minimo que indica que se a sobreposição entre a read e o adaptador for menor
# que esse valor, a read não será modificada.

# -q: Nesse parâmetro deve ser fornecido um valor de qualidade de bases, o qual abaixo desse as bases das extremidades
# serão removidas. Se quer que ambas as extremidades sejam trimadas, deve-se passar 2 valore, o primeiro para a extremidade
# 5’ e o segundo para a extremidade 3’. No comando acima apenas um valor foi passado: 25, assim apenas bases com score de
# qualidade abaixo de 25 da extreimidade 3’ serão trimadas.

# -T: Fornece o número de threads a ser utilizado para trimagem das reads.

# Correct-mismatches: Esse parâmetro indica como lidar com os mismatches durante o alinhamento.
# As opções são “Liberal, conservativa, e N”. O N significa substituir a base que deu mismatch por N.
# As opções de correção ‘Liberal e conservativa’  envolvem a susbtituição de bases por aquela que apresenta a maior
# qualidade. Elas diferem apenas quando os scores de qualidade são iguais. A correção liberal substitui de acordo com a
# melhor mediana das qualidades das bases, enquanto que a conservativa significa manter inalterada.

# Pair-filter: Nesse parâmetro deve-se escolher qual das reads (pair-end) deve corresponder ao critério de filtragem para
# que ocorra a filtragem. Você pode escolher que ambas as reads entrem no critério ou qualquer uma. No caso foi escolhido qualquer uma. 

# -a: Nesse parametro deve-se fornecer a sequência do adaptador a ser removido na primeira read na região 3’.
# -A: Nesse parametro deve-se fornecer a sequência do adaptador a ser removido na segunda read na região 3’.

# -o: Local de saída da read 1 trimada.
# -p: Local de saída da read 2 trimada.
# -pe1: Primeiro arquivo de entrada
# -pe2: Segundo arquivo de entrada

# Untrimmed-output: Indica o local de saída das reads 1 que não foram trimadas.
# Untrimmed-paired-output: Indica o local de saída das reads 2 que não foram trimadas.

	atropos trim --aligner insert \
            	     -e 0.1 \
                     -n 2 \
                     -m 15 \
                     --op-order GAWCQ \
                     --match-read-wildcards \
                     -O 25 \
                     -q 28 \
                     -T 8 \
                     --correct-mismatches conservative \
                     --pair-filter any \
                     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT  \
                     -o ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
                     -p ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
                     -pe1 ${r1} \
                     -pe2 ${r2} \
                     --untrimmed-output        ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
                     --untrimmed-paired-output ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \
                     > ${outdir}/processed/atropos/${name}.atropos.log.out.txt \
                     2> ${outdir}/processed/atropos/${name}.atropos.log.err.txt

# As reads que não foram trimadas na etapa anterior são submetidas a esse outro algoritmo "adapter"

	echo -e "Running atropos (adapter) for adapter trimming using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq ...\n"

# -pe1: É o arquivo de entrada da read 1 para a trimagem, entretanto são o arquivo de saída do comando anterior,
# que não foram trimadas pelo algoritmo ‘insert’.
# -pe2: É o arquivo de entrada da read 2 para a trimagem, entretanto são o arquivo de saída do comando anterior,
#  que não foram trimadas pelo algoritmo ‘insert’.

# -g:  Nesse parametro deve-se fornecer a sequência do adaptador a ser removido na primeira read na região 5’.
# -G: Nesse parametro deve-se fornecer a sequência do adaptador a ser removido na segunda read na região 5’.

	atropos trim  --aligner adapter \
                      -e 0.1 \
                      -n 2 \
                      -m 1 \
                      --match-read-wildcards \
                      -O 3 \
                      -q 28 \
                      --pair-filter both \
                      -pe1 ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
                      -pe2 ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \
                      -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
                      -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  \
                      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                      -G CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
                      -T 8 \
                      -o  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq  \
                      -p  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
                      >  ${outdir}/processed/atropos/${name}.atropos_adapter.log.out.txt \
                      2>  ${outdir}/processed/atropos/${name}.atropos_adapter.log.err.txt

# Após esses comandos, ambas as saidas dos dois comandos são concatenadas em um só arquivo, utlizando o comando cat:
# -read1 trimada com algoritmo insert + read 1 trimada com algoritmo adapter.
# -read 2 trimada com algoritmo insert + read 2 trimada com algoritmo adapter.

	echo -e "Merging atropos adapter trimming results using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq + ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq ...\n"

	cat       ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
        	  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
   		> ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq

	cat       ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
        	  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
   		> ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq

# Em seguida, os arquivos gerados antes da concantenação foram removidos pois não serão utilizados.

	echo -e "Removing useless atropos results ...\n"

	rm -f ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
	      ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
	      ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
	      ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq

# Em seguida as reads são submetidas ao programa PrinSeq, que é uma ferramenta para garantir que os dados obtidos não estejam
# comprometidos por sequencias de baixa qualidade, artefatos ou contaminantes que podem levar a conclusões errôneas.
# Essa ferramenta, realiza tanto a trimagem, filtragem de sequencias de baixa qualidade, quanto gera um resumo estatistico
# das caracteristicas das sequencias, como um controle de qualidade, parecido com o fastqc.

	echo -e "PrinSeq processing: ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq ...\n"

# Após as reads serem processadas pelo Atropos, servirão de entrada para o prinseq (paramentros:fastq e fastq2).
# -min_len: Filtrar as sequências menores que 25pb.
# -noniupac: Filtrar sequências que apresentem caracteres que não seja A,T,C,G ou N.

# -lc-method: Método a ser utilizado para filtragem de sequências de baixa complexidade, pode ser "dust" ou "entropy". No caso,
# foi selecionado o dust: Os scores são computados de acordo com a frequencia de diferentes trinucleotideos (0-100).
# Scores elevados implicam em baixa complexidade.

# -lc_threshold: Valor a ser usado para filtragem de sequencias de baixa complexidade, no caso, como escolhido o método
#  "dust", é escolhido um valor máximo aceito. Acima desse score, as sequencias serão filtradas.

# -out_format: Utilizado para escolher o formato do output. No comando abaixo foi selecionado o formato fastq.

# -out_good: Por default a saida de dados é o diretorio em que se encontram as reads. Para alterar o diretório de
# saída e nome do arquivo de saída utiliza-se esse parametro. No comando abaixo o diretório prinseq como output e o nome
# dos arquivos de saída.

# -out_bad: Por default, a saída é o diretório input. Para alterar deve-se usar esse parametro. No caso foi utilizado o null
# para que esses arquivos bad não sejam gerados.

# -trim_qual_window: Tamnho da "janela" a ser usado para calcular o score de qualidade. No caso, está dizendo para parar na
# primeira base que falhar.

# -trim_qual_step: De quanto em quanto a janela desliza avaliando os scores de qualidade sem faltar nenhum.
# -trim_qual-right: Trima as sequências que apresentam score de qualidade menor que o threshold na extremidade 3'.
# -trim_qual_type: Tipo de calculo a ser usado para se obter o score de qualidade, foi selecionado a média.
# -trim_qual_rule: Regra a ser utilizada para comparação do score de qualidade gerado.Foi selecionado lt (less than).

	prinseq-lite.pl -fastq  ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq \
			-fastq2 ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq \
			-out_format 3 \
			-trim_qual_window 1 \
			-trim_qual_step 1 \
			-trim_qual_right 30 \
			-trim_qual_type mean \
			-trim_qual_rule lt \
			-out_good ${outdir}/processed/prinseq/${name}.atropos_final.prinseq \
			-out_bad  null \
			-lc_method dust \
			-lc_threshold 30 \
			-min_len 25 \
			-trim_tail_right 5 \
			-trim_tail_left 5 \
			-ns_max_p 80 \
			-noniupac \
			 > ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.out.log \
			2> ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.err.log

# Realização do controle de qualidade após ter sido realizado o pre-processamento dos dados

	echo -e "FastQC pos-evaluation using sample ${name}: ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq & ${outdir}/processed/prinseq/${name}_2.atropos_final.prinseq.fastq ...\n"

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.err.txt

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}_2.atropos_final.prinseq.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.err.txt

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_1_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq \
	   	   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.err.txt
	fi

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_2_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq \
		   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.err.txt
	fi


done
