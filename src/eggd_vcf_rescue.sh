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

    if [[ ("$non_pass_rescue" != "true" && "$rescue_filtered" != "true") \
        || ("$non_pass_rescue" == "true" && "$rescue_filtered" == "true") ]]; then
        dx-jobutil-report-error "Error: one of --non-pass-rescue or --rescue-filtered
            are required and are not mutually exclusive. Please check inputs and re-run
            with only one specified as true."
        exit 1
    fi

    if [[ "$non_pass_rescue" == "true" \
       && ( -z $fasta_tar_name || -z $gvcf || -z $rescue_vcf ) ]]; then
        dx-jobutil-report-error "Error: non_pass_rescue specified but not all required files
            of fasta_tar, gvcf and rescue_vcf passed"
        exit 1
    fi

    if [[ "$rescue_filtered" == "true" \
       && ( -z $filtered_vcf || -z $unfiltered_vcf || -z $rescue_vcf ) ]]; then
        dx-jobutil-report-error "Error: rescue_filtered specified but not all required files
            of filtered_vcf, unfiltered_vcf and rescue_vcf passed"
        exit 1
    fi
}


_decompress() {
    : '''
    Decompresses given vcf if compressed, else returns 0

    Arguments:
        vcf to decompress

    Outputs
        decompressed vcf file
    '''
    if [[ "$input_vcf" == *.gz ]]; then
        gunzip "$input_vcf"
    fi
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
    local input_vcf="$1"

    if [[ "$input_vcf" == *.gz ]]; then
        gunzip "$input_vcf"
        input_vcf=${input_vcf/.gz/}
    fi

    bgzip "$input_vcf"
    bcftools index "${input_vcf}.gz"
}


_get_sample_prefix() {
    : '''
    Gets prefix from samplename, splits on underscores and takes first field and
    also strips .vcf.gz if no underscore present

    Arguments
        file : filename to get prefix from

    Returns
        string : file prefix to stdout
    '''
    sed 's/.vcf\|.gz//g' <<< "$(cut -d'_' -f1 <<< "$1")"
}


_strip_chr_prefix() {
    : '''
    Strips chr prefixes from given vcf file, will first decompress file if compressed

    Arguments
        file : vcf to strip prefixes from
        string : name for output file

    Outputs
        decompressed vcf with chr prefixes removed
    '''
    local input_vcf="$1"
    local outname="$2"

    _decompress "$input_vcf"
    awk '{gsub(/chr/,""); print}' "${input_vcf/.gz/}" > "$outname"
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

    # # ensure all VCFs are bgzipped and indexed
    # _compress_and_index "$gvcf_name"
    # _compress_and_index "$rescue_vcf_name"

    # set names of compressed VCFs incase they now differ to app input
    # gvcf_name="${gvcf_name/.gz/}.gz"
    # rescue_vcf_name="${rescue_vcf_name/.gz/}.gz"

    # Remove reference calls from gvcf
    bcftools view -m2 "$gvcf_name" -o "${sample_prefix}.vcf"

    # Remove chr prefix from sample vcf
    _strip_chr_prefix "${sample_prefix}.vcf" "${sample_prefix}_noChr.vcf"

    # Normalise and left align filtered vcf
    bcftools norm -m -any -f genome.fa "${sample_prefix}_noChr.vcf" \
        -o "${sample_prefix}_norm.vcf"

    # Create a vcf of all NON-PASS which match the rescue sites
    bcftools filter -i 'FILTER!="PASS"' "${sample_prefix}_norm.vcf"   \
        | bcftools filter -m + -s "$filter_tag" --mask-file "${rescue_vcf_name}" - \
        | bcftools filter -i "FILTER~\"${filter_tag}\"" - -Oz -o "${sample_prefix}.rescued.vcf.gz"

    # Create a vcf with only PASS variants
    bcftools view -f PASS "${sample_prefix}_norm.vcf" -Oz -o "${sample_prefix}_pass.vcf.gz"

    # index vcf files
    bcftools index "${sample_prefix}_pass.vcf.gz"
    bcftools index "${sample_prefix}.rescued.vcf.gz"

    # Concatenate flagged non-pass variant vcf with PASS vcf
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
    # get prefix from splitting on underscore, strip .vcf.gz in case of no underscore
    local sample_prefix
    sample_prefix=$(_get_sample_prefix "$unfiltered_vcf_name")


    if [[ "$strip_chr" == 'true' ]]; then
        # remove chr prefix from sample vcf
        filtered_outname=$(_get_sample_prefix $filtered_vcf_name).filtered.noChr.vcf
        unfiltered_outname=$(_get_sample_prefix $unfiltered_vcf_name).unfiltered.noChr.vcf
        rescue_outname=$(_get_sample_prefix $rescue_vcf_name).noChr.vcf

        _strip_chr_prefix "$filtered_vcf_name" "$filtered_outname"
        _strip_chr_prefix "$unfiltered_vcf_name" "$unfiltered_outname"
        _strip_chr_prefix "$rescue_vcf_name" "$rescue_outname"

        filtered_vcf_name="$filtered_outname"
        unfiltered_vcf_name="$unfiltered_outname"
        rescue_vcf_name="$rescue_outname"
    fi

    # normalise and left align vcfs
    bcftools norm -m -any -f genome.fa "${filtered_vcf_name}" \
        -o "${sample_prefix}.filtered.norm.vcf"
    bcftools norm -m -any -f genome.fa "${unfiltered_vcf_name}" \
        -o "${sample_prefix}.unfiltered.norm.vcf"

    filtered_vcf_name="${sample_prefix}.filtered.norm.vcf"
    unfiltered_vcf_name="${sample_prefix}.unfiltered.norm.vcf"


    # rescue variants against given rescue vcf
    rescue_vcf="${sample_prefix}.rescued.vcf.gz"
    bcftools filter -m + -s "$filter_tag" --mask-file "$rescue_vcf_name" "$unfiltered_vcf_name" \
        | bcftools filter -i "FILTER~\"${filter_tag}\"" - -Oz -o "$rescue_vcf"

    _compress_and_index "$filtered_vcf_name"
    _compress_and_index "$rescue_vcf"

    # combine filtered vcf with rescued variants to output
    outname="${sample_prefix}.vcf.gz"
    bcftools concat -a -d all "${filtered_vcf_name}.gz" "$rescue_vcf" -Oz -o "$outname"

    # Upload output vcf
    output_vcf=$(dx upload "$outname" --brief)
    dx-jobutil-add-output output_vcf "$output_vcf" --class=file
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
    dx-download-all-inputs --parallel
    find ~/in/ -type f -name "*" -print0 | xargs -0 -I {} mv {} /home/dnanexus

    # Unpack fasta tar if given
    if [[ $fasta_tar_name ]]; then time tar -I pigz -xf "$fasta_tar_name"; fi

    mark-section "Rescuing variants"

    if [[ "$non_pass_rescue" == "true" ]]; then
        # passed gvcf to process
        _rescue_non_pass
    fi

    if [[ "$rescue_filtered" == "true" ]]; then
        # rescuing variants from a filtered VCF and unfiltered vcf
        _rescue_filtered
    fi

    echo "Done!"
    mark-success
}
