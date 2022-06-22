#!/bin/bash
# eggd_tso500_rescue 1.0.0

set -exo pipefail

main() {

    echo "Value of gvcf: '$gvcf'"
    echo "Value of hotspot_vcf: '$hotspot_vcf'"

    # Download input files
    dx download "$gvcf" & \
    dx download "$hotspot_vcf" & \
    dx download "$fasta_tar"

    # Unpack fasta tar
    tar xzf $fasta_tar_name

    # Get sample prefix
    sample_prefix=${gvcf_name/_MergedSmallVariants.genome.vcf/}

    # Remove reference calls from gvcf
    bcftools view -m2 $gvcf_name -o ${sample_prefix}.vcf

    # Remove chr prefix from patient vcf to match the reference genome
    # This app was created for output of TSO500 which only contains chr1-22,X & Y
    # It should be reaccessed for different input vcfs.
    awk '{gsub(/chr/,""); print}' ${sample_prefix}.vcf > ${sample_prefix}_noChr.vcf

    # Normalise and left align filtered vcf
    bcftools norm -m -any -f genome.fa ${sample_prefix}_noChr.vcf \
        -o ${sample_prefix}_norm.vcf

    # Create a vcf of all NON-PASS which match the OPA hotspots
    bcftools filter -i 'FILTER!="PASS"' ${sample_prefix}_norm.vcf.gz   \
        | bcftools filter -m + -s 'OPA' --mask-file ${hotspot_vcf_name} - \
        | bcftools filter -i 'FILTER~"OPA"' - -Oz -o ${sample_prefix}.rescued.vcf.gz

    # Create a vcf with only PASS variants
    bcftools view -f .,PASS ${sample_prefix}_norm.vcf  -Oz \
    -o ${sample_prefix}_pass.vcf.gz

    # Zip and index vcf files to use with bcftools isec command
    bcftools index ${sample_prefix}_pass.vcf.gz
    bcftools index ${sample_prefix}.rescued.vcf.gz


    # Concatenate OPA flagged non-pass variant vcf with pass vcf
    bcftools concat -a  ${sample_prefix}_pass.vcf.gz ${sample_prefix}.rescue.vcf.gz \
        -Oz --threads $(nproc --all) -o ${sample_prefix}_withLowSupportHotspots.vcf.gz

    # Upload output vcf
    filtered_vcf=$(dx upload ${sample_prefix}_withLowSupportHotspots.vcf.gz --brief)
    dx-jobutil-add-output filtered_vcf "$filtered_vcf" --class=file

    echo "Done"
}
