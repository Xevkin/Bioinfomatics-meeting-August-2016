#BEAST
##Bayesian Evolutionary Analysis by Sampling Trees 
BEAST performs Bayesian phylogenetic analysis on molecular sequences. Current version is BEAST2. Data for BEAST must be prepared using Beauti.

BEAST can:

1. investigate phylogenies using different clock models (eg relaxed, strict) and site models
2. estimate mutation rates and divergence times
3. produce Bayesian phylogenetic trees with posterior probablilites and 95% HPDs 
4. take into account sampling time of sequences (e.g. ancient samples; viruses sampled over several decades)
5. take into account known priors (e.g. mutation rates; split time with outgroup)
6. partition data (e.g. into coding and non-coding, different genes etc)
7. produce Bayesian skyline plots

This is just the base prograe BEAST. Additional modules allow you to:

1. *BEAST: allows estimation of species trees using multilocus data
2. Perform phylogeographic analysis
3. Use SNP and AFLP data
4. Use fossil sample data
5. Compare different site models (by calculating the Bayes Factor)


##Before BEAST: Model Selection
A good idea before running Beauti/BEAST is to run a Model Selection program such as modelgenerator or jmodeltest. These will suggest a site model that best fits with your data, and other parameters that can be used as priors in BEAST (e.g. the proportion of invariant sites; the parameters of the gamma distribution; kappa, the transition/transversion rate ratio).

##Data Format
To prepare the BEAST input file (.xml), you must first use Beauti. The input to Beauti is the .nexus format. The ASSUMPTIONS block allows you to define partitions (and other details). This looks like:


	#NEXUS

	begin data;
		dimensions ntax=88 nchar=16647;
		format datatype=dna missing=? gap=- interleave;
	matrix
	A1_01            GTTGATGTAGCTTAAACTTAAAGCAAGGCACTGAAAATGCCTAGATGAGTGTACCAACTCCATAAACACA
	A1a_02           GTTGATGTAGCTTAAACTTAAAGCAAGGCACTGAAAATGCCTAGATGAGTGTACCAACTCCATAAACACA
	A1a_03           GTTGATGTAGCTTAAACTTAAAGCAAGGCACTGAAAATGCCTAGATGAGTGTACCAACTCCATAAACACA
	[...]
	goat_reference   ATCTTTACTCCAGCCAAGGTAAATATATAAGTGCCTGGGTCTTTTACATGGTAAGTG
	bezoar_reference ATCTTTACTCCAGCCAAGGTAAATATATAAATGCCTGGGTCTTTTAC----------

	;
	end;

	begin ASSUMPTIONS;

		charset coding = 2743-3697 3908-4949 5331-6875 7018-7701 7773-9394 9464-9809 9880-11547 11750-14081 14155-15292;

		charset RNA = 1-640 1026-2740 3699-3836 3839-3907 4950-5016 5018-5086 5088-5160 5194-5329 6876-6941 6949-7016 7705-7771 9395-9463 9811-9879 11548-11678 11680-11749 14082-14150 15298-15432;  

		charset D_loop = 15433-16647;

		charset remainder = 641-1025 2741-2742 3698 3837-3838 5017 5087 5161-5193 5330 6942-6948 7017 7702-7704 7772 9810 11679 14151-14154 15293-15297;

	END;


To generate a nexus file, first run Multiple Sequence Alignment (e.g. MUSCLE, CLUSTAL) on the samples of interest.Then convert the resulting fasta file to nexus using an online tool (e.g. [here] (http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_nexus.php)). From here you may need to add an assumptions block.

Other notes on input data:

1. BEAST ignores missing data.
2. BEAST core programme using single locus data - *BEAST for multilocus 

##Running Beauti
![Beauti](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/beauti.png)
To prepare the .xml input for BEAST, we run Beauti. Beauti is broken into several tabs:

Partitions 
------ 
Allows us to link or unlink the tree, clock or site models of any partitions we might have.

Tip Dates 
------ 
Allows us to define the age of our samples (set "Before Present"). Note that "Height" does not seem to work. Tip sampling allows you to estimate the age of samples but requires some manual set up not through Beauti.

Site Model 
------
The Site Model allows you to define how genetic variation the mutation rate is expected to change across the sequences in question. The default model is a gamma distribution. Gamma category count describes how much rate heterogenity you expect/allow in your data (1 mean none; 4 is sufficient). Shape also describes rate variation (see https://github.com/BEAST2-Dev/MGSM/wiki). You can also define/provide a prior for the proportion of invariant sites in the data. 

Finally you can define the specific nucleotide substitution model for the data, and any parameters it might have. 

We can define/allow different site models for different partitions.

Clock Model
------
Here we define how we expect the mutation rate to change across the phylogeny. Can define a prior mutation rate, and set different models and priors for each partition.

1. Strict Clock: constant molecular clock
2. Relaxed Clock Log Normal: clock rate can change along branches
3. Relaxed Clock Exponential: clock rate can change at nodes
4. Random local clock: different parts of the phylogeny can have their own clock rate

Priors
------
You can further modify priors here. For example you can set minimum and maximum values for particular priors, and set distributions for the priors. For the tree prior, we can set trees more suitable to within-species dataset (Coalescent Contstant Population) or between-species (Yule Model).

MCMC
------
Here we define the MCMC settings for the run. We want to perform a run with good mixing (no autocorelation of parameter estimates, varying around a line like a hairy caterpillar) that gives good Effective Sample Sizes (ESS - effectively independent draws from the posterior distribution) of each parameter (>100, but the higher the better). 

1. Chain Length is the number of times the MCMC algorithim estimates the parameters of the model. The number of chains needed varies from dataset to dataset. 50 million is what I find useful for ~85 whole mitochondria.
2. Pre Burn-In is the number of chains we run and then discard before the run starts properly. Random starting state may have very low posterior probability, so we run some number then discard them (e.g. 1 million). Burn-in can also be defined at a later stage (logannotator).
3. Trace Log: how often we record the parameter estimates in the .log output files. We want about 10,000 samples (so for 50M chains we record every 5,000). 
4. Tree log: similar, but records the tree estimate in a .tree output file.

Best practise is to run several identical runs (e.g. 2, 4) and combine the output (using logcombiner). We want different runs to equilibrate at the same end point.

##Running BEAST
BEAST itself can be run through the command line or via a sparse GUI. Command line allows you to give BEAST additional threads. This step is long - 50 million chains with 4 threads may take ~1 day.

##Analyzing output
The BEAST output is assessed using tracer. If multiple identical runs have been performed, merge .log and .tree using logcombiner, then load the log files in tracer. We want to see good mixing of the parameters ("hairy caterpillar" mcmc traces) and high ESS for parameters (preferably >200). This suggests the run has went to convergence - it is sampling parameters from a stationary distribution.

####What you want to see
![What you want to see](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/good_trace.png)


####What you don't want to see
![What you don't want to see](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/bad_trace.png)

We can view the mean estimate of each parameter, and other statistics such as the 95% Highest Posterior Distribution (shortest interval in the paramter space that contains 95% of the posterior probability ).

###Trees
Tree files can be combined using logcombiner, and should be summarized using treeannotator - allows you to define number of tree to be burned-in, and the tree heights to be reported (should use maximum clade credibilty tree and median heights). Visualize using figtree - can see 95% HPD for node heights or branch lengths, etc.

![node HPDs](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/tree_HPD.png)

###Resources
BEAST2: http://beast2.org/
BEAST usergroup: https://groups.google.com/forum/#!forum/beast-users
