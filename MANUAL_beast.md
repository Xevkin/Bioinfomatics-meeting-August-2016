#BEAST
##Bayesian Evolutionary Analysis by Sampling Trees
BEAST performs Bayesian phylogenetic analysis on molecular sequences. It can:

1. investigate phylogenies using different clock models (eg relaxed, strict) and site models
2. estimate mutation rates and divergence times
3. produce Bayesian phylogenetic trees with posterior probablilites and 95% HPDs 
4. take into account sampling time of sequences (i.e. ancient samples)
5. take into account known priors (e.g. mutation rates; split time with outgroup)
6. partition data (i.e. into coding and non-coding, different genes etc)
7. produce Bayesian skyline plots

This is just the base programme BEAST. Additional modules allow you to:

1. *BEAST: allows estimation of species trees using multilocus data
2. Perform phylogeographic analysis
3. Use SNP and AFLP data
4. Use fossil sample data
5. Compare different site models (by calculating the Bayes Factor)


###Data Format
BEAST accepts the .nexus format. The ASSUMPTIONS block allows you to define partitions (and other details). This looks like:


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


To generate a nexus file, first run Multiple Sequence Alignment (e.g. MUSCLE, CLUSTAL) on the samples of interest.Then convert the resulting fasta file to nexus using an online tool (e.g. [here] (http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_nexus.php)).



###Running BEAST

