#MSMC
#Multiple Sequential Markovian Coalescent model
Developed by Stephan Schiffels and Richard Durbin.

Extends on PSMC on several ways:
* Base model (SMC') is better - it accurately estimates the recombination rate, which PSMC does not
* Can model more than two sequences - will have more coalescence events in the recent past (<20ky), so can more 
accurately estimate recent population history.
* Exploits phase information

MSMC overcomes the computational complexity of determining the underlying tree by simplifying the relationships between samples at a given locus to:

1. The TMRCA of any two sequences i.e. most recent coalescence, and the identity of the sequences
2. The total length of all singleton branches in the tree - branches that give rise to variants of minor allele count 1

![MSMC image 1](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/ng.3015-F1.jpg)

As more sequences are added, the TMRCA tends to decrease -> more recent population inferences.

![MSMC image 2](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/ng.3015-F2.jpg)

**NOTE:** As the number of sequences increases, there is an increasing bias towards smaller ancestral population size.

**NOTE:** We also see the "smoothing" of instantaineous events over thousands of years with MSMC.

MSMC uses the "relative cross coalescence rate", the ratio between the cross-population and within-population coalescence rates, to measure genetic seperation between two populations. It should be close to 1 when mixed and close to 0 when fully seperated.

MSMC *can* be run with unphased sequences; however it shows some biases at either end. The relative cross coalescence rate in particular is affected when using out-of-phased data.

![MSMC image 3](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/ng.3015-SF5.jpg)


##Running MSMC
1. You must generate a mappability mask for your reference if you have not done so before (one for humans should already be available). This is both time and space consuming. See makeMappabilityMask.py from msmc-tools and SNPable - http://lh3lh3.users.sourceforge.net/snpable.shtml

2. Starting with a bam file (one for each sample), run bamCaller.py (from msmc-tools), producing a vcf file and a mask file. You must provide a mean coverage estimate:
`samtools mpileup -q 20 -Q 20 -C 50 -u -r <chr> -f <ref.fa> <bam> | bcftools call -c -V indels |\`  
`./bamCaller.py <mean_cov> <out_mask.bed.gz> | gzip -c > <out.vcf.gz>`
3. The bed.gz sample mask file may need to be converted to a bed.txt.gz file in the correct format for msmc - this involves some reshaping of files.
4. For each chromosome, you must run generate_multihetsep.py from msmc-tools to generate the msmc input files. You must have the sample vcf and mask files, as well as a reference mapability mask. If running with multiple samples, add additional vcf and sample mask files. For more a single individual, the samples must all be phased.
`./generate_multihetsep.py --mask=covered_sites_sample1_chr1.bed.txt.gz \`                  `--mask=mappability_mask_chr1.bed.txt.gz sample1_chr1.vcf.gz >sample1_chr1.msmc`
5. Then run msmc (single sample):
`./msmc_linux_64bit -o  sample1 <all_msmc_files_for_autosomes>`
For multiple samples, add the --fixedRecombination flag.

You can also investigate population seperation at this stage, by defining subpopulations. In msmc1, these must be of equal size (and only two subpops). Refer to the msmc guide for this.

R code to visualize output:
`mu<-1e-8`

`gen<-5`

`sample1<-read.table(“<sample1.final.txt>”, header=T)`

`sample2<-read.table(“<sample1.final.txt>”, header=T)`

`plot(sample1$left_time_boundary/mu*gen, (1/sample1$lambda)/mu, log="x", ylim=c(0,100000),type="n", xlab="Years ago", ylab="effective population size")`

`lines(sample1$left_time_boundary/mu*gen, (1/sample1$lambda)/mu, type="s", col="blue")`

`lines(sample2$left_time_boundary/mu*gen, (1/sample2$lambda)/mu, type="s", col="green")`

`legend("topleft", legend=c("SOM","HAR", "HL", "Bed3", "Hxh2"), col=c("red", "blue", "green", "black", "orange"), lty=c(1,1,1,1,1))`

**NOTE:** You must scale the output by the mutation rate at this point.

![example msmc image](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/all_samples_2016-6-9.jpeg)

##Resources
* MSMC paper: http://www.nature.com/ng/journal/v46/n8/full/ng.3015.html
* MSMC github: https://github.com/stschiff/msmc
* MSMC guide/tutorial: https://github.com/stschiff/msmc/blob/master/guide.md
* MSMC tools: https://github.com/stschiff/msmc-tools
