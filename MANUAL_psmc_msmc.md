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
* Masking coding regions (lower mutation rate) gives a *slightly* different history (shape is the same) -> not masking coding regions has minimal effect on estimations.
* Coalesence occures less frequently in a structured population. Fewer coalescences is equiavelent to a larger effective population size -> inflates population size estimate for structured populations.
* Can run on female X chromosomes - must scale plot by 4/3.
* Combine two X chromosomes from two male individuals to create pseudo-diploid X chromosomes -> inference on timing of population splits.
![PSMC image 3](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/psmc.png)

###Running PSMC
From the psmc github page:

  samtools mpileup -C50 -uf ref.fa aln.bam | bcftools view -c - \
  vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz

  utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa
  psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa
  utils/psmc2history.pl diploid.psmc | utils/history2ms.pl > ms-cmd.sh #you do not need this step (generates ms input)
  utils/psmc_plot.pl diploid diploid.psmc


##Papers
PSMC paper: http://www.nature.com/nature/journal/v475/n7357/full/nature10231.html
