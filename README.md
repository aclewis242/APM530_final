# DNA SEQUENCE MODEL

### General info
An object-oriented implementation of DNA sequence evolution, written with generality in mind -- i.e., it should not be difficult to expand the capabilities of the code in the future as necessary. Nucleotide bases (A, C, G, T) are stored as integer indices (0, 1, 2, 3 respectively). The substitution models use a continuous-time Markov chain algorithm to generate a new sequence from an ancestral sequence after a certain amount of time has passed.

### Model info
- *JC69: Jukes & Cantor 1969*

    The simplest nucleotide substitution model. It assumes that all nucleotide bases are equally frequent and equally likely to be replaced by a different base. Though these assumptions are biologically unrealistic, they render the model extremely easy to implement and use, which makes it useful for ensuring the code base (objects, methods, etc.) is working as intended.

- *K80: Kimura 1980*

    A slightly more complex model allowing for different rates of transitions (purine -> purine, pyrimidine -> pyrimidine) and transversions (purine <-> pyrimidine).

### Running the code
Originally written for Python 3.11.5. If using Anaconda, these should all be installed by default, but if not, the required packages are:
- NumPy
- SciPy
- Matplotlib

To run the code, enter `python main.py` into the console with an activated conda environment.

### Citations
* K. Strimmer, A. von Haeseler, "Genetic distances and nucleotide substitution models." *The Phylogenetic Handbook* (2009), ch. 4, Cambridge University Press.