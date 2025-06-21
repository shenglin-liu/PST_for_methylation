# PST_for_methylation
Filter away SNPs that are out of Hardy-Weinberg equilibrium (HWE) from a VCF file.

## Usage
On a linux system (R must be installed), download the script to any folder you desire.

Change the execution mode by running: 
`chmod a+x [path_to]/VCF_HWF.r`

Run the script by typing: 
`[path_to]/VCF_HWF.r $vcf_file $popmap_file $min_pop_HWE $out_vcf`

* vcf_file: the input file in VCF format.
* popmap_file: a population map file in the same format as in [Stacks](https://catchenlab.life.illinois.edu/stacks/manual/).
* min_pop_HWE: minimum number of populations in HWE for a SNP to be kept.
* out_vcf: name of the new VCF file to output.

## Requirements for the input
* The genotype fields of the VCF file must start with GT tag.
* The SNPs must be biallelic, with 0 representing REF and 1 for ALT.
* There is no limitation for the size of the VCF file.

## To cite
