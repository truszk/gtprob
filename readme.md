# gtprob: Rapid computation of gene tree probabilities for monophyletically concordant gene trees under the multispecies coalescent #

This program can be used to compute the probability of a monophyletically concordant gene tree topology given a species tree. Alternatively, given a species tree and the number of sampled genes from each leaf species, it can be used to compute the probability that the gene tree is monophyletically concordant with the species tree. It implements an algorithm described in Truszkowski, Scornavacca and Pardi[1]. Monophyletic concordance was first defined by Rosenberg[2].

## Requirements: ##

- dendropy (>=4.4.0)
- numpy
- scipy

All of the above can be installed using the requirements file by typing

`pip install -r requirements.txt`

## Usage: ##

To compute the probability of a gene tree given a species tree, type

`python ./calc_prob_concordant.py gene_tree_file species_tree_file`

where `gene_tree_file` is the gene tree topology in Newick format (see **File Formats** below) and `species_tree_file` is the species tree with branch lengths in coalescent units. 

To compute the probability of monophyletic concordance given a species tree, type

`python ./calc_prob_concordant.py gene_samples_file species_tree_file`

where `gene_samples_file` specifies the number of samples from each species.

##File formats##

A **species tree file** is a standard Newick file with branch lengths in coalescent units.

A **gene tree file** should contain a gene tree topology in Newick format, without the branch lengths. Each sample from `species_x` should have `species_x` as a substring of its name.

A **gene samples file** specifies the number of samples from each species in the following format:

	species_1 num_samples_in_species_1
	species_2 num_samples_in_species_2
	...
	species_n num_samples_in_species_n

## Please cite: ##

[1] Truszkowski J., Scornavacca C.,Pardi F. Computing the probability of gene trees concordant with the species tree in the multispecies coalescent. Submitted.

## References: ##

[1] Truszkowski J., Scornavacca C.,Pardi F. Computing the probability of gene trees concordant with the species tree in the multispecies coalescent. Submitted.

[2] Rosenberg N. The probability of topological concordance of gene trees and species trees.Theoretical Population Biology, 61(2):225–247, 2002.


