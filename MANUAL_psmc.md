#PSMC
PSMC utilizes the distribution of Times since the Most Recent Common Ancestor (TMRCA) between the two alleles at loci across the genome (autosomes and X) to estimate the demographic history of the population that individual is from. 
* Recombination breaks up large chunks that have a common TMRCA
* Whole genomes -> multiple loci -> more detailed TMRCA distribtuion -> insight into older demographic events. 

![PSMC image 1](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/nature10231-f1.2.jpg)

###Parameters
* Mutation rate (eg 10^-8)
* Generation time

###Limitations of PSMC
* Poor resolution for events <20yka (humans). Few recombinations events in that time period -> high variance and inaccurate
* Poor resolution of rapid events. May spread bottlenecks over tens of thousands of years preceding event, etc.
![PSMC image 2](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/nature10231-f2.2.jpg)
* False Positives (false heterozygotes - sequencing errors or ancient damage) pushes the estimate further along the X axis, and greatly affects recent population estimates (overestimates).
* False Negatives (heterozygotes called as homozygous reference - due to low coverage) pull the estimate towards the origin. However, it is roughly equivalent to a lower neutral mutation rate, so can be corrected for.
* Hypermutation regions (eg segmental duplications) can overestimated ancient population sizes, but rest of curve is unaffected.

###Other features
* Robust to recombination hotspots
* Robust to large stretches of ambiguous bases ("N"s)
* Masking coding regions (lower mutation rate) gives a nearly-identical history 
* Coalesence occures less frequently in a structured population. Fewer coalescences is equiavelent to a larger effective population size -> inflates population size estimate for structured populations.
* Can run on female X chromosomes - must scale plot by 4/3.
* Combine two X chromosomes from two male individuals to create pseudo-diploid X chromosomes -> inference on timing of population splits.
![PSMC image 3](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/psmc.png)

###Running PSMC

First we must generate a .fq file that represent bins of sequence and whether a heterozygote is found there.

`samtools mpileup -C50 -uf ref.fa aln.bam | bcftools view -c -O v -o - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz`

**Note:**Here we exclude all sites with less than a depth of 10 and any with more than 100. In general, I exclude sites with less than a third the mean coverage, and more than twice the mean coverage. So it varies from sample to sample.

**Note:**The version of samtools I'm using here is v1.2. If you use the incorrect version, it will cause problems in the pipeline.

We convert this fq file to a psmcfa file:

`utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa`

We then run psmc. From the psmc github page: "In particular, the `-p' option specifies that there are 64 atomic time intervals and 28 (=1+25+1+1) free interval parameters. The first parameter spans the first 4 atomic time intervals, each of the next 25 parameters spans 2 intervals, the 27th spans 4 intervals and the last parameter spans the last 6 time intervals. The '-p' and '-t' options are manually chosen such that after 20 rounds of iterations, at least ~10 recombinations are inferred to occur in the intervals each parameter spans. Impropriate settings may lead to overfitting. The command line in the example above has been shown to be suitable for modern humans."

`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa`

Finally we plot the trace. During this we should account for false negatives (missed heterozygotes) by scaling the mutation rate. To estimate the FNR, I calculate the proportion of the genome that falls within the the -d and -D options from the previous step eg the proportion of the genome falling within 10 and 100 depth (yes, this is imperfect). I do this using GATK's DepthOfCoverage.

`utils/psmc_plot.pl -M "sample=0.1" output_name diploid.psmc`

You can also plot multiple samples at the same time:

`utils/psmc_plot.pl -M "sample1=0.1, sample2=0.15" output_name sample1.psmc sample2.psmc`

You can also perform bootstrapping - sample subsequences with replacement, estimating demographic history for these subsequences (do not try to account for FNR!):

`utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa`
`utils/splitfa diploid.psmcfa > split.psmcfa`
`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa`
`seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc split.fa | sh`
`cat diploid.psmc round-*.psmc > combined.psmc`
`utils/psmc_plot.pl -pY50000 combined combined.psmc`

###Coverage
Orlando (2013) downsampled a 21.1X horse genome, and determined the FNR correction needed to return to the original demographic history.
* 15.8X - ~0.06
* 10.5X - ~0.26 - however older events become less accurate

Nadachowska-Brzyska (2016) developed guidelines for PSMC usage, baed on their analysis of ~200 Flycatcher genomes.

* Mean genome coverage of ≥18X
* Per-site minimum coverage (-d) of 10
* No more than 25% missing data (following filtering)


##Resources
* PSMC paper: http://www.nature.com/nature/journal/v475/n7357/full/nature10231.html
* PSMC github: https://github.com/lh3/psmc
* PSMC coverage discussion in Orlando 2013: http://www.nature.com/nature/journal/v499/n7456/extref/nature12323-s1.pdf
* PSMC coverage discussion in Nadachowska-Brzyska (2016): http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4793928/
* My PSMC script (may need modifying): https://github.com/Xevkin/scripts_for_goat_project/blob/master/run_psmc.sh
* 
