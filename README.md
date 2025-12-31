# 3rd-Order Markov Gene Classifier

A probabilistic model for distinguishing protein-coding DNA sequences from shuffled
coding DNA using a third-order Markov approach.

## Overview
This project implements a third-order Markov chain model to classify whether a DNA
sequence originates from a protein-coding gene or from a randomly shuffled version
of a coding sequence. By modeling local nucleotide dependencies, the classifier
assigns log-likelihood scores that reflect gene-like sequence structure.

## Training Data
The model is trained using three known protein-coding genes:
- recA
- SRY
- ACE2

For each gene, a shuffled version of the sequence is generated and used as
negative training data. Coding sequences serve as positive samples, while shuffled
sequences serve as controls.

## Testing and Evaluation
Model performance is evaluated using ten independent protein-coding genes:
- BDNF
- CYP2C19
- HLA-DRB1
- IL6
- INS
- LRRK2
- MAPT
- MYC
- STAT3
- TNF

Each test gene is paired with its shuffled counterpart as a control group.
Log-likelihood ratios are computed for both groups.

Statistical significance is assessed using Welch’s *t*-test to compare the
log-likelihood distributions between gene and control sequences.

## Results
The results are summarized in `TestResult.html`. Log-likelihood ratios for both
gene and control groups follow approximately normal distributions. Welch’s
*t*-test indicates a statistically significant difference between protein-coding
sequences and shuffled controls, demonstrating the model’s ability to distinguish
gene-like sequence structure.

## Acknowledgement
The programs are partly wrote and debugged from Chat GPT.
The gene sequences are downloaded from NCBI.
