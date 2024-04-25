# DNA SEQUENCE MODEL

### General info
An object-oriented implementation of DNA sequence evolution, written with generality in mind -- i.e., it should not be difficult to expand the capabilities of the code in the future as necessary. Nucleotide bases (A, C, G, T) are stored as integer indices (0, 1, 2, 3 respectively). The substitution models use a continuous-time Markov chain algorithm to generate a new sequence from an ancestral sequence after a certain amount of time has passed.

### Running the code
Originally written for Python 3.11.5. If using Anaconda, these should all be installed by default, but if not, the required packages are:
- NumPy
- SciPy
- Matplotlib

To run the code, enter `python main.py` into the console with an activated conda environment. Un/comment certain lines in the main method to run the various simulations (refer to the file for further information).

### Citations
* H. Honma et al., "Mutation tendency of mutator *Plasmodium berghei* with proofreading-deficient DNA polymerase δ" (2016), Nature. Retrieved from https://www.nature.com/articles/srep36971.
* N. Lanchier, *Stochastic Modeling* (2017), Springer, p. 101-124.
* P. Lemey et al., *The Phylogenetic Handbook: A Practical Approach to Phylogenetic Analysis and Hypothesis Testing* (2009), Cambridge University Press, p. 111-123.