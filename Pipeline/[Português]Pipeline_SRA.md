# Coorte HBOC — Configuração do Ambiente, Download de Amostras e Criação do Arquivo BED

Este guia apresenta o passo a passo utilizado para configurar um ambiente de trabalho em um servidor Linux com GPU, baixar arquivos FASTQ de amostras HBOC (Síndrome de Câncer de Mama e Ovário Hereditário) do NCBI, construir um arquivo BED das regiões de interesse (BRCA1, BRCA2, TP53) e realizar processamento, alinhamento e chamada de variantes com o NVIDIA Parabricks.

## Pré-requisitos

- Servidor Linux com `bash`, `wget` e acesso à internet para `repo.anaconda.com`, `eutils.ncbi.nlm.nih.gov` e `trace.ncbi.nlm.nih.gov`
- Acesso `sudo` (para instalação do Conda e permissões de grupo do Docker)
- GPU(s) NVIDIA com drivers instalados, Docker e NVIDIA Container Toolkit (para a etapa do Parabricks)
- Licença/acesso ao container do NVIDIA Parabricks (`nvcr.io/nvidia/clara/clara-parabricks`)

---

## 1. Instalar o Conda

Se o `conda` ainda não estiver disponível:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

Verifique a instalação:

```bash
conda info --envs
```

---

## 2. Criar o ambiente `sra-tools`

```bash
conda create -n sra-tools -c bioconda -c conda-forge sra-tools entrez-direct -y
conda activate sra-tools
```

Verifique a instalação:

```bash
prefetch --version
fasterq-dump --version
esearch -help
```

---

## 3. Resolver o accession BioSample → SRA Run

As amostras selecionadas fazem parte do estudo "Familial history and prevalence of BRCA1, BRCA2 and TP53 pathogenic variants in HBOC Brazilian patients from a public healthcare service" do Instituto Nacional de Câncer (INCA) [PMID: 36329109; PMCID: PMC9633799] e estão disponíveis no BioProject [PRJDB13663](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB13663).

Você pode baixar individualmente os arquivos a partir de um accession de BioSample ou de amostra SRA (ex: `DRS309369`), encontre o(s) accession(s) de run associado(s):

```bash
esearch -db sra -query "DRS309369" | efetch -format runinfo
```

Isso retorna uma tabela no formato CSV, incluindo a coluna `Run` (ex: `DRR378923`), layout da biblioteca, plataforma, BioProject e metadados da amostra.

> **Dica:** avisos ocasionais do tipo `curl: (56) OpenSSL SSL_read` vindos dos servidores do NCBI durante essa etapa geralmente são resets de conexão transitórios e não indicam falha — verifique se a tabela de saída esperada ainda foi impressa logo abaixo dos avisos.

Repita para cada amostra de interesse. Confirme os campos `BioSample` e `SRA` de cada resultado para evitar confundir pacientes/amostras diferentes que tenham números de accession parecidos.

Ou você pode baixar todos os runs do estudo com o acesso do BioProject (ex: `PRJDB13663`)

```bash
esearch -db sra -query "PRJDB13663" | efetch -format runinfo > all_runs.csv
```

O arquivo CSV gerado terá todos os accessions e os metadados de cada amostra. Confira em seguida quantas e quais amostras existem antes de baixar. Recomendável verificar se existe espaço em disco disponível para armazenar os dados. Pode ser que demore um pouco.

```bash
wc -l all_runs.csv
cut -d',' -f1 all_runs.csv | tail -n +2
```

---

## 4. Baixar e converter os arquivos FASTQ

## 4.1 Baixar arquivos em lote
```bash
mkdir -p fastq_output
tail -n +2 all_runs.csv | cut -d',' -f1 | while read -r ACC; do
    echo "Baixando ${ACC}..."
    prefetch "${ACC}"
    fasterq-dump "${ACC}" --split-files -O ./fastq_output
done
```

## 4.2 Baixar arquivos individualmente
```bash
prefetch <accession_SRR_ou_DRR>
fasterq-dump <accession_SRR_ou_DRR> --split-files -O ./fastq_output
```

- `--split-files` separa os reads pareados em `_1.fastq` e `_2.fastq`
- Adicione `--gzip` depois, ou use `gzip` via pipe, para comprimir a saída (o `fasterq-dump` não comprime nativamente)
- Se os dados estiverem sob acesso controlado do dbGaP, será necessária uma chave de repositório `.ngc` aprovada:
  ```bash
  prefetch <accession> --ngc /caminho/para/seu.ngc
  ```

---

## 5. Construir um arquivo BED dos genes de interesse (BRCA1, BRCA2, TP53)

### 5.1 Exportar coordenadas de genes/éxons do UCSC Table Browser

1. Acesse o [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
2. Selecione o genoma de referência: **GRCh38/hg38** (precisa ser igual à sua referência FASTA)
3. Group: *Genes and Gene Predictions* → Track: *GENCODE V49* (ou *RefSeq Curated*)
4. Em "Define regions of interest," clique em "Paste list" e informe os nomes dos genes: `BRCA1`, `BRCA2`, `TP53`
5. Em `Output format` selecione `BED` (ou a tabela estendida no estilo genePred, como usado abaixo)
6. Adicione um `Output filename` (Ex: `Genes_of_interest`)
7. Baixe o arquivo

> **Importante:** a exportação padrão do UCSC é uma tabela no estilo genePred (colunas como `name`, `chrom`, `exonStarts`, `exonEnds`), não um arquivo BED padrão de 3–12 colunas. Ela precisa ser convertida antes de ser usada com `samtools`/`bedtools`.

### 5.2 Converter para formato BED por éxon

```bash
awk 'BEGIN{OFS="\t"}
NR>1 {
  split($9, starts, ",");
  split($10, ends, ",");
  for(i=1; i<=$8; i++) {
    print $2, starts[i], ends[i], $1"_exon"i, ".", $3
  }
}' Genes_of_interest.bed | sort -k1,1 -k2,2n > BRCA_TP53_exons.bed
```

### 5.3 Mesclar éxons sobrepostos entre transcritos

```bash
bedtools merge -i BRCA_TP53_exons.bed -s -c 4,6 -o distinct,distinct > BRCA_TP53_merged.bed
```

### 5.4 Verificação rápida

```bash
wc -l BRCA_TP53_merged.bed
cut -f1 BRCA_TP53_merged.bed | sort -u
```

Esperado: apenas registros em `chr13` (BRCA2) e `chr17` (BRCA1, TP53). Qualquer outro cromossomo indica um erro de parsing em alguma etapa anterior.

---

## 6. Configurar acesso Docker + GPU para o Parabricks

### 6.1 Verificar visibilidade da GPU dentro do container

```bash
docker run --rm --gpus all nvcr.io/nvidia/clara/clara-parabricks:<tag> nvidia-smi
```

### 6.2 Verificar o `pbrun` e a montagem dos volumes

```bash
docker run --rm --gpus all nvcr.io/nvidia/clara/clara-parabricks:<tag> pbrun version

docker run --rm --gpus all \
    --volume /caminho/para/fastq_output:/workdir \
    --volume /caminho/para/output_dir:/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:<tag> \
    ls -lh /workdir /outputdir
```

> Substitua `<tag>` pela versão do Parabricks instalada (confirme com o teste do `nvidia-smi` acima, ou com `docker images`).

---

## 7. Rodar o pipeline germline do Parabricks

Um script em lote (`run_germline.sh`) automatiza isso por amostra:

- Descobre os pares de amostra (`<accession>_1.fastq`, `<accession>_2.fastq`)
- Pula amostras já concluídas (verifica se o VCF de saída já existe)
- Roda `pbrun germline` com a referência, os VCFs de known-sites e o arquivo BED de intervalos
- Registra em log o tempo de execução e o status de saída de cada amostra

> **Nota:** a saída do Docker é redirecionada para um arquivo de log (`>> "${LOG}" 2>&1`) em vez de ser impressa em tempo real. Para acompanhar o progresso ao vivo, use `tail -f` no arquivo de log a partir de uma sessão de terminal separada.

---

## Referências

Matta BP, Gomes R, Mattos D, Olicio R, Nascimento CM, Ferreira GM, Brant AC, Boroni M, Furtado C, Lima V, Moreira MÂM, Dos Santos ACE. Familial history and prevalence of BRCA1, BRCA2 and TP53 pathogenic variants in HBOC Brazilian patients from a public healthcare service. Sci Rep. 2022 Nov 3;12(1):18629. doi: 10.1038/s41598-022-23012-3. PMID: 36329109; PMCID: PMC9633799.

https://docs.nvidia.com/clara/parabricks/about-parabricks 

https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
