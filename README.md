<!-- dx-header -->
# eggd_vcf_rescue (DNAnexus Platform App)
-----------------------------------------
## What does this app do?

Rescues variants against a given vcf of positions.

This app may be run in 2 modes:

- `rescue_filtered`: When run in this mode an unfiltered and filtered vcf must be provided, variants will be rescued from the unfiltered vcf against the given rescue vcf and concatenated with the filtered vcf.
- `rescue_non_pass`: When run in this mode a gvcf must be provided, PASS and non-PASS variants will be extracted from the gvcf (with reference calls removed), and variant positions present in the given rescue_vcf rescued against the non-PASS variants.


## What are typical use cases for this app
To rescue variants filtered out, or of low quality, against positions to never exclude in a rescue vcf (i.e. mutation hotspots, known pathogenic variants etc.)


## What data are required for this app to run?
**Files**:

- `rescue_vcf` (required): vcf of known sites to resuce against
- `fasta_tar` (required): tar of reference fasta and index
- `gvcf` (`rescue_non_pass` mode): gvcf to extract PASS variants and rescue non-PASS variants from
- `filtered_vcf` (`rescue_filtered` mode): vcf of filtered variants to concatenate rescued variants with
- `unfiltered_vcf` (`rescue_filtered` mode): vcf of unfiltered sites to rescue variants from


**Modes**:

- `rescue_filtered` (`bool`): when run in this mode an unfiltered and filtered vcf must be provided, variants will be rescued from the unfiltered vcf and concatenated with the filtered vcf. Mutually exclusive with `rescue_non_pass`.
- `rescue_non_pass` (`bool`): when run in this mode a gvcf must be provided, PASS and non-PASS variants will be extracted from the gvcf, and rescued against the non-PASS variants. Mutually exclusive with `rescue_filtered`.

**Optional**:

- `strip_chr` (`bool`): if true, will strip chr prefixes from input vcfs. Should be specified if given reference fasta does not contain chr prefixes.
- `filter_tag` (`string`): tag to add to `FILTER` field of rescued variants (default: `rescued`), this will be appended to any existing FILTER fields (`bcftools filter -m +`).


## What does this app output?

- **rescue_non_pass mode**: a vcf file with all PASS variants and rescued low quality variants from the provided rescue vcf.
- **rescue_filtered mode**: a vcf of variants from the given `filtered_vcf` combined with variants rescued from `unfiltered_vcf`.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.
