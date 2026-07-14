#!/bin/bash
set -euo pipefail

FASTQ_DIR="/home/mkolarik/Practice/fastq_output"
REF_DIR="/nfs/theseus/Pipelines_parabricks/gatk_ref"
REFERENCE="${REF_DIR}/Homo_sapiens_assembly38.fasta"
KNOWN_SITES_1="${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
KNOWN_SITES_2="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_3="${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz"
INTERVAL_FILE="/nfs/theseus/Local_Collaborators/mlkolarik/Training/BRCA_TP53_merged.bed"
BRCA_DIR="/nfs/theseus/Local_Collaborators/mlkolarik/Training/"
OUTPUT_DIR="${BRCA_DIR}/output/Germline"
IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

# Descoberta de amostras: arquivos no formato DRRxxxxxx_1.fastq / DRRxxxxxx_2.fastq
mapfile -t SAMPLES < <(ls "${FASTQ_DIR}"/*_1.fastq 2>/dev/null | xargs -I{} basename {} | sed 's/_1\.fastq$//' | sort -u)
#SAMPLES=("DRR378918")  # amostra menor, ideal para smoke test

TOTAL=${#SAMPLES[@]}
echo "Amostras: ${TOTAL} | Modo: sequencial | Inicio: $(date '+%Y-%m-%d %H:%M:%S')"
CONCLUIDAS=0
FALHAS=0

for SAMPLE in "${SAMPLES[@]}"; do
    LOG="${LOG_DIR}/${SAMPLE}.log"
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}"

    if [[ -f "${SAMPLE_OUT}/${SAMPLE}.vcf" ]]; then
        echo "SKIP: ${SAMPLE}"
        CONCLUIDAS=$((CONCLUIDAS+1))
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Iniciando $((CONCLUIDAS+FALHAS+1))/${TOTAL}: ${SAMPLE}" | tee -a "${LOG}"

    R1="${FASTQ_DIR}/${SAMPLE}_1.fastq"
    R2="${FASTQ_DIR}/${SAMPLE}_2.fastq"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "ERRO: sem FASTQ para ${SAMPLE}"
        FALHAS=$((FALHAS+1))
        continue
    fi

    START=$(date +%s)
    docker run --rm \
        --gpus all \
        --volume "${FASTQ_DIR}:${FASTQ_DIR}:ro" \
        --volume "${REF_DIR}:${REF_DIR}:ro" \
        --volume "${BRCA_DIR}:${BRCA_DIR}:ro" \
        --volume "${SAMPLE_OUT}:${SAMPLE_OUT}" \
        "${IMAGE}" pbrun germline \
        --ref "${REFERENCE}" \
        --in-fq "${R1}" "${R2}" \
        --knownSites "${KNOWN_SITES_1}" \
        --knownSites "${KNOWN_SITES_2}" \
        --knownSites "${KNOWN_SITES_3}" \
        --interval-file "${INTERVAL_FILE}" \
        --out-bam "${SAMPLE_OUT}/${SAMPLE}.bam" \
        --out-variants "${SAMPLE_OUT}/${SAMPLE}.vcf" \
        --out-recal-file "${SAMPLE_OUT}/${SAMPLE}.txt" \
        --num-gpus 8 \
        >> "${LOG}" 2>&1
    EC=$?
    EL=$(( $(date +%s) - START ))
    FMT="$((EL/3600))h $(((EL%3600)/60))m $((EL%60))s"

    if [[ ${EC} -eq 0 ]]; then
        CONCLUIDAS=$((CONCLUIDAS+1))
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK ${CONCLUIDAS}/${TOTAL} (${FMT}): ${SAMPLE}" | tee -a "${LOG}"
    else
        FALHAS=$((FALHAS+1))
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FALHOU exit=${EC} (${FMT}): ${SAMPLE}" | tee -a "${LOG}"
    fi
done

echo "Fim: $(date '+%Y-%m-%d %H:%M:%S') | OK: ${CONCLUIDAS} | Falhas: ${FALHAS} | Total: ${TOTAL}"
