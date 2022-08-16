#!/bin/bash

set -exo pipefail

_validate_inputs() {
    : '''
    Check given input files to ensure correct ones passed for non pass
    recuing or rescuing filtered variants

    Arguments:
        None

    Returns:
        0 if inputs all valid, 1 if:
            - non_pass_rescue and rescue_filtered booleans both specified
            - required files for non_pass_rescue mode not provided
            - required files for resuce_filtered mode not provided
    '''
    mark-section "Validating inputs"

    if [[ ("$non_pass_rescue" != "True" && "$rescue_filtered" != "True") \
        || ("$non_pass_rescue" == "True" && "$rescue_filtered" == "True") ]]; then
        dx-jobutil-report-error "Error: one of --non-pass-rescue or --rescue-filtered
            are required and are not mutually exclusive. Please check inputs and re-run
            with only one specified."
        exit 1
    fi

    if [[ "$non_pass_rescue" == "True" \
       && ( -z $fasta_tar_name || -z $gvcf || -z $rescue_vcf ) ]]; then
        dx-jobutil-report-error "Error: non_pass_rescue specified but not all required files
            of fasta_tar, gvcf and rescue_vcf passed"
        exit 1
    fi

    if [[ "$rescue_filtered" \
       && ( -z $filtered_vcf || -z $unfiltered_vcf || -z $rescue_vcf ) ]]; then
        dx-jobutil-report-error "Error: rescue_filtered specified but not all required files
            of filtered_vcf, unfiltered_vcf and rescue_vcf passed"
        exit 1
    fi
}


_rescue_non_pass() {
    : '''
    Rescues non-pass variants from gvcf against given vcf of variants to retain, uploads
    VCF of PASS + rescued variants and exits app.


    Globals:
        gvcf: given gVCF to rescue non-PASS variants from
        filter_tag: tag to apply to FILTER field of rescued variants

    Arguments:
        None

    Outputs:
        None
    '''
    local sample_prefix
    sample_prefix=$(cut -d'_' -f1 <<< "$gvcf_name")

    # ensure all VCFs are bgzipped and indexed
    _compress_and_index "$gvcf_name"
    _compress_and_index "$rescue_vcf_name"

    # set names of compressed VCFs incase they now differ to app input
    gvcf_name="${gvcf_name/.gz/}.gz"
    rescue_vcf_name="${rescue_vcf_name/.gz/}.gz"

    # Remove reference calls from gvcf
    bcftools view -m2 "$gvcf_name" -o "${sample_prefix}.vcf"

    # Remove chr prefix from sample vcf
    awk '{gsub(/chr/,""); print}' "${sample_prefix}.vcf" > "${sample_prefix}_noChr.vcf"

    # Normalise and left align filtered vcf
    bcftools norm -m -any -f genome.fa "${sample_prefix}_noChr.vcf" \
        -o "${sample_prefix}_norm.vcf"

    # Create a vcf of all NON-PASS which match the rescue sites
    bcftools filter -i 'FILTER!="PASS"' "${sample_prefix}_norm.vcf"   \
        | bcftools filter -m + -s "$filter_tag" --mask-file "${rescue_vcf_name}" - \
        | bcftools filter -i "FILTER~\"${filter_tag}\"" - -Oz -o "${sample_prefix}.rescued.vcf.gz"

    # Create a vcf with only PASS variants
    bcftools view -f PASS "${sample_prefix}_norm.vcf"  -Oz -o "${sample_prefix}_pass.vcf.gz"

    # index vcf files
    bcftools index "${sample_prefix}_pass.vcf.gz" "${sample_prefix}.rescued.vcf.gz"

    # Concatenate flagged non-pass variant vcf with pass vcf
    bcftools concat -a "${sample_prefix}_pass.vcf.gz ${sample_prefix}.rescued.vcf.gz" \
        -Oz -o "${sample_prefix}_withLowSupportHotspots.vcf.gz"

    # Upload output vcf
    output_vcf=$(dx upload "${sample_prefix}_withLowSupportHotspots.vcf.gz" --brief)
    dx-jobutil-add-output output_vcf "$output_vcf" --class=file
}


_rescue_filtered() {
    : '''
    Rescues filtered out variants from a given unfiltered - filtered vcf pair
    against given VCF of variants to retain, uploads rescued VCF and exits app.

    Globals:
        filtered_vcf : vcf of filtered variants to append to for output vcf
        unfiltered_vcf : vcf unfiltered variants used to rescue from
        rescue_vcf : vcf of known variants to rescue against
        filter_tag : tag to apply to FILTER field of rescued variants

    Arguments:
        None

    Outputs:
        None
    '''
    local sample_prefix
    sample_prefix=$(cut -d'_' -f1 <<< "$unfiltered_vcf")

    # ensure all VCFs are bgzipped and indexed
    _compress_and_index "$filtered_vcf_name"
    _compress_and_index "$unfiltered_vcf_name"
    _compress_and_index "$rescue_vcf_name"

    # set names of compressed VCFs incase they now differ to app input
    filtered_vcf_name="${filtered_vcf_name/.gz/}.gz"
    unfiltered_vcf_name="${filtered_vcf_name/.gz/}.gz"
    rescue_vcf_name="${rescue_vcf_name/.gz/}.gz"

    # remove chr prefix from sample vcf
    awk '{gsub(/chr/,""); print}' "$unfiltered_vcf_name" > "${sample_prefix}_noChr.vcf"

    # normalise and left align filtered vcf
    bcftools norm -m -any -f genome.fa "${sample_prefix}_noChr.vcf" \
        -o "${sample_prefix}_norm.vcf"

    # rescue variants against given rescue vcf
    bcftools filter -m + -s "$filter_tag" --mask-file "$rescue_vcf_name" "${sample_prefix}_norm.vcf" \
        | bcftools filter -i "FILTER~${filter_tag}" - -Oz -o "${sample_prefix}.rescued.vcf.gz"

    # combine unfiltered vcf with rescued variants to output
    outname=$(sed 's/.vcf\|.gz//g' <<< "$filtered_vcf_name").vcf.gz
    bcftools concat -a -d all "$filtered_vcf_name" "${sample_prefix}.rescued.vcf.gz" -Oz -o "$outname"

    # Upload output vcf
    output_vcf=$(dx upload "$outname" --brief)
    dx-jobutil-add-output output_vcf "$output_vcf" --class=file

}


_compress_and_index() {
    : '''
    Decompresses given vcf if already compressed, then compresses with
    bgzip and indexes with bcftools index

    Arguments:
        - vcf to compress and index

    Outputs:
        - compressed vcf + index
    '''
    input_vcf="$1"

    if [[ "$input_vcf" == *.gz ]]; then
        gunzip "$input_vcf"
        input_vcf=${input_vcf/.gz/}
    fi

    bgzip "$input_vcf"
    bcftools index "${input_vcf}.gz"
}


main() {
    echo "Inputs specified:"
    echo "Value of gvcf: '$gvcf'"
    echo "Value of filtered_vcf: '$filtered_vcf'"
    echo "Value of unfiltered_vcf: '$unfiltered_vcf'"
    echo "Value of rescue_vcf: '$rescue_vcf'"
    echo "Value of fasta tar: '$fasta_tar_name'"
    echo "Value of boolean non_pass_rescue: '$non_pass_rescue'"
    echo "Value of boolean rescue_filtered: '$rescue_filtered"
    echo "Value of FILTER tag: '$filter_tag'"

    _validate_inputs

    mark-section "Downloading inputs"

    # build string of given files and download
    input_files=$(tr -s '[:space:]' <<< "${gvcf} ${filtered_vcf} ${unfiltered_vcf} ${rescue_vcf} ${fasta_tar_name}")
    xargs -P 4 -n1 dx download <<< "$input_files"

    # Unpack fasta tar if given
    if [[ -z $fasta_tar_name ]]; then tar xzf "$fasta_tar_name"; fi


    mark-section "Rescuing variants"

    if [[ "$non_pass_rescue" == "True" ]]; then
        # passed gvcf to process
        _rescue_non_pass
    fi

    if [[ "$rescue_filtered" == "True" ]]; then
        # rescuing variants from a filtered VCF and unfiltered vcf
        _rescue_filtered
    fi

    echo "Done!"
    mark-success
}