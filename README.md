# ErwiniaCRISPRAligner
### Written by: Michael Parcey

ErwiniaCRISPRAligner is a series of scripts which are used to extract and align CRISPR spacers from the genome of Erwinia amylovora. This pipeline was created for a specific study will not be supported or developed past submission. The scripts are not compiled so they can be freely used and modified.
If used or modified, please reference:
Parcey M, Gayder S, Castle AJ, Svircev AM. 2022. Function and Application of the CRISPR-Cas System in the Plant Pathogen Erwinia amylovora. Appl Envir Micro. e02513-21. doi: 10.1128/aem.02513-21

## REQUIREMENTS

### The scripts are written in Biopython and require to following to run:
* Python3
* Numpy (Biopython requirement)
* Biopython

### The scripts require the sequences under consideration to be in fasta format with the extension ".fasta". Each genome requires its own file, and there should only be one, concatinated sequence per file. All the sequence files being analyzed need to be in one folder.

## Usage

### Place the CRISPRScripts.tar.gz and the CRISPRAligner.sh in the folder containing the sequence files. From terminal, execute:
### "bash CRISPRAligner.sh Erwinia" 
### You can replace Erwinia with whatever you want the prefix of the data files to be.

## DATA
This code will produce a folder named Results. In this folder will be several files and folders:

* Erwinia.fasta - contains all the sequences used.
* Erwinia.Error.fasta - contains any spacers which were excluded because they were too large or contain null base calls.
* Erwinia.Results.csv - spreadsheet containing how many sequences were found in each cassette for each strain.

* CRR1, CRR2, CRR4 are all folders which contain data files of the spacers for a specific array.
* Erwinia.CRR1.fasta - contains all of the CRISPR spacers, the strain they are from, and their genomic positions.
* Erwinia.CRR1Aligned.csv - a visual representation of the spacer alignment in a spreadsheet. The top line is the consensus sequence and the each strain will show which CRISPR region and spacer it matches with (e.g. CRR1S2).
* Erwinia.CRR1Aligned.fasta - the alignment of the sequences for all of the spacers. Can be analyzed to SNPs (MEGA X).
* Erwinia.CRR1Binary.phy - a phylip file of the spacer alignment in binary format. Can be used with IQ Tree to produce a phylogeny of the CRISPR array. Note phylip format only allows 10 characters in the name of the sequence.
* Erwinia.CRR1Unique.fasta - all of the unique CRISPR spacers which were found in this array.

* The AllCRR folder contains the same data files as the CRR folders but the data represents all of the CRISPR cassettes ordered CRR1, CRR2, then CRR4.
