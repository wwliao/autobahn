# Autobahn
Variant Calling Suite for Cancer Genomics

Given a cancer and matched normal pair, will call somatic variants. 
SNVs
* Varscan
* GATK

Small Indels
* Pindel

Structural Variants
* Delly2
* Lumpy
* Manta

Once variants have been called, will aggregate all results into a single file for easy downstream analysis

`autobahn_config` provides necessary path to various softwares used. 
