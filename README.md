<!-- dx-header -->
# eggd_tso500_rescue (DNAnexus Platform App)
-----------------------------------------
## What does this app do?

It takes the genomic vcf which is outputed by the TSO500 local app v2.2 and filters out the reference calls and low support variants but keeps low support variants that are present in a given whitelist.

## What are typical use cases for this app
During the validation of the TSO500 local app v2.2 analysis pipeline, some known variants were excluded during filtering. This app attempts to rescue in filtered out variants if they are a known variant of interest given a whitelist of variants.

## What data are required for this app to run?
Required inputs for this app:
* a genomic gvcf (.genome.vcf)
* a compressed whitelist vcf containing variants needed to always be filtered in (*.vcf.gz)

## What does this app output?
This app outputs:
* a vcf file with all PASS variants in addition to filtered in low quality variants given a whitelist of variants provided.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.