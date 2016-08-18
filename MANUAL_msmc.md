#MSMC - Multiple Sequential Markovian Coalescent model
Developed by Stephan Schiffels and Richard Durbin.

Extends on PSMC on several ways:
* Base model (SMC') is better (transition and emission rate calculations)
* Can model more than two sequences - will have more coalescence events in the recent past (<20ky), so can more 
accurately estimate recent population history.
* Exploits phase information

MSMC overcomes the computational complexity of determining the underlying tree by simplifying the relationships between samples at a given locus to:
1. The TMRCA of any two sequences i.e. most recent coalescence, and the identity of the sequences
2. The total length of all singleton branches in the tree - branches that give rise to variants of minor allele count 1

![MSMC image 1](https://github.com/Xevkin/Bioinfomatics-meeting-August-2016/blob/master/ng.3015-F1.jpg)

##Resources
* MSMC paper: http://www.nature.com/ng/journal/v46/n8/full/ng.3015.html
* MSMC github: https://github.com/stschiff/msmc
* MSMC guide/tutorial: https://github.com/stschiff/msmc/blob/master/guide.md
