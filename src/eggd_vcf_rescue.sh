#!/bin/bash

set -exo pipefail

_validate_inputs() {
    : '''
    Check given input files to ensure correct ones passed for non pass
    recuing or rescuing filtered variants. Also ensures no invalid characters
    are present in filter_tag and filter_tag_description to be written to
    vcf header.

    Globals
        rescue_non_pass : boolean to run in non_pass rescue mode
        rescue_filtered : boolean to run in rescue filter mode
        gvcf : given gVCF to rescue non-PASS variants from
        filtered_vcf : vcf of filtered variants to append to for output vcf
        unfiltered_vcf : vcf of unfiltered variants used to rescue from

    Arguments:
        None

    Returns:
        0 if inputs all valid, 1 if:
            - rescue_non_pass and rescue_filtered booleans both specified
            - required files for rescue_non_pass mode not provided
            - required files for resuce_filtered mode not provided
    '''
    mark-section "Validating inputs"

    if [[ $filter_tag ]]; then
        filter_tag=$(sed  -e 's/[ ;:=\]/_/g' -e "s/[\'\"]//g"  <<< "$filter_tag")
    fi

    if [[ $filter_tag_description ]]; then
        filter_tag_description=${filter_tag_description//\"/\'}
    fi

    if [[ ("$rescue_non_pass" != "true" && "$rescue_filtered" != "true") \
        || ("$rescue_non_pass" == "true" && "$rescue_filtered" == "true") ]]; then
        dx-jobutil-report-error "Error: one of --non-pass-rescue or --rescue-filtered
            are required and are not mutually exclusive. Please check inputs and re-run
            with only one specified as true."
        exit 1
    fi

    if [[ "$rescue_non_pass" == "true" && -z $gvcf ]]; then
        dx-jobutil-report-error "Error: rescue_non_pass specified but not all required files
            of fasta_tar, gvcf and rescue_vcf passed"
        exit 1
    fi

    if [[ "$rescue_filtered" == "true" && ( -z $filtered_vcf || -z $unfiltered_vcf) ]]; then
        dx-jobutil-report-error "Error: rescue_filtered specified but not all required files
            of filtered_vcf, unfiltered_vcf and rescue_vcf passed"
        exit 1
    fi
}


_decompress() {
    : '''
    Decompresses given vcf if compressed, else returns 0

    Arguments:
        file : vcf to decompress

    Outputs
        file : decompressed vcf file
    '''
    if [[ "$input_vcf" == *.gz ]]; then
        pigz -d "$input_vcf"
    fi
}


_compress_and_index() {
    : '''
    Decompresses given vcf if already compressed, then compresses with
    bgzip and indexes with bcftools index

    Arguments:
        file : vcf to compress and index

    Outputs:
        files : compressed vcf + index
    '''
    local input_vcf="$1"

    _decompress "$input_vcf"
    input_vcf=${input_vcf/.gz/}

    bgzip "$input_vcf"
    bcftools index -f "${input_vcf}.gz"
}


_modify_header() {
    : '''
    Modifies vcf header after calling bcftools filter to update FILTER description
    added for the given filter_tag to link to current DNAnexus job

    Arguments:
        file : vcf file with applied filter_tag

    Outputs:
        file : vcf file with modified header
    '''
    local vcf="$1"

    # get header from vcf
    bcftools view -h "$vcf" > header.txt

    # find ##FILTER line in header added by bcftools filter, add to description
    header_line=$(grep "^##FILTER=<ID=${filter_tag}" header.txt)
    description_addition="variants rescued against given variant positions "
    description_addition+="in eggd_vcf_rescue (DNAnexus job: $DX_JOB_ID)"
    if [[ "$filter_tag_description" ]]; then description_addition+=". ${filter_tag_description}"; fi

    modified_header_line=${header_line/\">/; $description_addition\">}
    sed -i "s/$header_line/$modified_header_line/" header.txt

    # write updated header back to the file
    bcftools reheader -h header.txt -o "tmp.vcf.gz" "$vcf"
    mv "tmp.vcf.gz" "$vcf"
}


_get_sample_prefix() {
    : '''
    Gets prefix from samplename, splits on underscores and takes first field.
    Will fall back to just stripping .vcf.gz if no underscore present.

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
    bcftools annotate --rename-chrs chr_prefix_map.txt "${input_vcf/.gz/}" > "$outname"
}

_filter_variants() {
    : '''
    Function: filters VCF based on filter_string inputted by the user

    Globals:
        filtered_vcf_name : vcf of filtered variants to further filter
        filter_string: string that will filter variant
    '''

    # check number of enteries before variant quality filtering
    num_var=$(grep -v ^"#"  "${filtered_vcf_name}" | wc -l)
    echo "VCF has $num_var variants before filtering variants"

    _compress_and_index $filtered_vcf_name
    eval ${filter_string} "${filtered_vcf_name}.gz" -o "$filtered_vcf_name"

    # check number of enteries after variant quality filtering
    rm "${filtered_vcf_name}.gz"
    num_var=$(grep -v ^"#"  "${filtered_vcf_name}" | wc -l)
    echo "VCF has $num_var variants after filtering variants"
}


_rescue_non_pass() {
    : '''
    Rescues non-pass variants from gvcf against given vcf of variants to retain, uploads
    VCF of PASS + rescued variants and exits app.


    Globals:
        gvcf: given gVCF to rescue non-PASS variants from
        rescue_vcf : vcf of known variants to rescue against
        filter_tag: tag to apply to FILTER field of rescued variants

    Arguments:
        None

    Outputs:
        None
    '''
    local sample_prefix
    sample_prefix=$(_get_sample_prefix "$gvcf_name")

    # Remove reference calls from gvcf
    bcftools view -m2 "$gvcf_name" -o "${sample_prefix}.vcf"

    if [[ "$strip_chr" == 'true' ]]; then
        # remove chr prefix from sample vcf and rescue vcf
        _strip_chr_prefix "${sample_prefix}.vcf" "${sample_prefix}.tmp.vcf"
        mv "${sample_prefix}.tmp.vcf" "${sample_prefix}.vcf"

        _strip_chr_prefix "$rescue_vcf_name" "tmp.vcf"
        mv tmp.vcf "${rescue_vcf_name/.gz/}"
        rescue_vcf_name="${rescue_vcf_name/.gz/}"  # overwrite global variable to new file
    fi

    # Normalise and left align filtered vcf
    bcftools norm -m -any -f genome.fa "${sample_prefix}.vcf" -o "${sample_prefix}_norm.vcf"

    # Create a vcf of all NON-PASS which match the rescue sites
    bcftools filter -i 'FILTER!="PASS"' "${sample_prefix}_norm.vcf"   \
        | bcftools filter -m + -s "$filter_tag" --mask-file "${rescue_vcf_name}" - \
        | bcftools filter -i "FILTER~\"${filter_tag}\"" - -Oz -o "${sample_prefix}.rescued.vcf.gz"

    # sense check for logs how many variants were rescued
    echo "Total variants rescued: $(zgrep -v '^#' ${sample_prefix}.rescued.vcf.gz | wc -l)"

    # modified description for FILTER field added by bcftools filter to better explain
    # provenance of filter_tag
    _modify_header "${sample_prefix}.rescued.vcf.gz"

    # if variant filtering is set using BCFtools filtering command,
    # then run the filter_variants command
    if [[ $filter_string ]]; then
        filtered_vcf_name="${sample_prefix}.rescued"
        _filter_variants
    fi


    # Create a vcf with only PASS variants
    bcftools view -f PASS "${sample_prefix}_norm.vcf" -Oz -o "${sample_prefix}_pass.vcf.gz"

    # index vcf files
    bcftools index "${sample_prefix}_pass.vcf.gz"
    bcftools index "${sample_prefix}.rescued.vcf.gz"

    # Concatenate flagged non-pass variant vcf with PASS vcf
    bcftools concat -a "${sample_prefix}_pass.vcf.gz" "${sample_prefix}.rescued.vcf.gz" \
        -Oz -o "${sample_prefix}_withLowSupportHotspots.vcf.gz"

    # Upload output vcf
    output_vcf=$(dx upload "${sample_prefix}_withLowSupportHotspots.vcf.gz" --brief)
    dx-jobutil-add-output output_vcf "$output_vcf" --class=file
}


_rescue_filtered() {
    : '''
    Rescues filtered out variants from a given unfiltered - filtered vcf pair
    against given VCF of variants to rescued, uploads concatenated VCF and exits app.

    Globals:
        filtered_vcf : vcf of filtered variants to append to for output vcf
        unfiltered_vcf : vcf of unfiltered variants used to rescue from
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
        filtered_outname=$(_get_sample_prefix $filtered_vcf_name).filtered.noChr.vcf
        unfiltered_outname=$(_get_sample_prefix $unfiltered_vcf_name).unfiltered.noChr.vcf
        rescue_outname=$(_get_sample_prefix $rescue_vcf_name).rescue.noChr.vcf

        # remove chr prefix from sample vcf
        _strip_chr_prefix "$filtered_vcf_name" "$filtered_outname"
        _strip_chr_prefix "$unfiltered_vcf_name" "$unfiltered_outname"
        _strip_chr_prefix "$rescue_vcf_name" "$rescue_outname"

        # overwrite global variables for convenience to point to new files
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

    # if variant filtering is set using BCFtools filtering command,
    # then run the filter_variants command
    if [[ $filter_string ]]; then
        _filter_variants
    fi

    # rescue variants against given rescue vcf
    rescue_vcf="${sample_prefix}.rescued.vcf.gz"
    bcftools filter -m + -s "$filter_tag" --mask-file "$rescue_vcf_name" "$unfiltered_vcf_name" \
        | bcftools filter -i "FILTER~\"${filter_tag}\"" - -Oz -o "$rescue_vcf"

    # modified description for FILTER field added by bcftools filter to better explain
    # provenance of filter_tag
    _modify_header "$rescue_vcf"

    # sense check for logs how many variants were rescued
    echo "Total variants rescued: $(zgrep -v '^#' $rescue_vcf | wc -l)"

    # compress and index for bcftools concat
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
    echo "Value of boolean rescue_non_pass: '$rescue_non_pass'"
    echo "Value of boolean rescue_filtered: '$rescue_filtered"
    echo "Value of FILTER tag: '$filter_tag'"
    echo "Value of FILTER string: '$filter_string'"

    _validate_inputs

    mark-section "Downloading inputs"
    dx-download-all-inputs --parallel
    find ~/in/ -type f -name "*" -print0 | xargs -0 -I {} mv {} /home/dnanexus
    tar -I pigz -xf "$fasta_tar_name"

    mark-section "Rescuing variants"

    if [[ "$rescue_non_pass" == "true" ]]; then
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
