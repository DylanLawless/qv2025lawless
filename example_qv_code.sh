#!/bin/bash

# Source master settings (including VQSR) and custom QV1 settings
source ./variables_master.sh
source ./variables_qv1.sh

# Run VQSR for SNPs

# 1. Calculate VQSLOD tranches for SNPs using VariantRecalibrator
gatk --java-options "${JAVA_OPTS}" VariantRecalibrator \
-R ${REF} \
-V ${vcf_file} \
--resource:hapmap,known=${vqsr_snp_hapmap_known},training=${vqsr_snp_hapmap_training},truth=${vqsr_snp_hapmap_truth},prior=${vqsr_snp_hapmap_prior} ${hapmap} \
--resource:omni,known=${vqsr_snp_omni_known},training=${vqsr_snp_omni_training},truth=${vqsr_snp_omni_truth},prior=${vqsr_snp_omni_prior} ${omni} \
--resource:1000G,known=${vqsr_snp_1000g_known},training=${vqsr_snp_1000g_training},truth=${vqsr_snp_1000g_truth},prior=${vqsr_snp_1000g_prior} ${thousandG} \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--mode SNP \
-O ${OUTPUT_DIR}/chr${INDEX}_snp1.recal \
--tranches-file ${OUTPUT_DIR}/chr${INDEX}_output_snp1.tranches

# 2. Filter SNPs on VQSLOD using ApplyVQSR

gatk --java-options "${JAVA_OPTS}" ApplyVQSR \
...continued