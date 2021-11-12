#!/bin/bash

# Written by: Michael Parcey (mp17ll@brocku.ca)

# Create Folders
if [ -d $PWD/Results ]; then
	rm -rf $PWD/Results
fi
mkdir -p $PWD/Results
mkdir -p $PWD/Results/AllCRR
mkdir -p $PWD/Results/CRR1
mkdir -p $PWD/Results/CRR2
mkdir -p $PWD/Results/CRR4

tar -zxvf CRISPRScripts.tar.gz

# Combine Fasta Files
echo "Combining Fasta Files..."
FILES="$PWD/*.fasta"
for f in $FILES
	do
	echo "Processing $f file..."
	python3 WriteFasta.py $f $1
	done
	
# Find CRISPR
echo "Finding CRISPR Spacers..."
python3 EaCRRFinder.py $1".fasta" $1"."

# Align Spacers
echo "Aligning CRISPR Spacers..."
python3 CRRAligner.py "CRR1/"$1".CRR1.fasta" "CRR1/"$1".CRR1Aligned.csv" "CRR1/"$1".CRR1Aligned.fasta" "CRR1/"$1".CRR1Unique.fasta" "CRR1/"$1".CRR1Binary.phy" "CRR1"
python3 CRRAligner.py "CRR2/"$1".CRR2.fasta" "CRR2/"$1".CRR2Aligned.csv" "CRR2/"$1".CRR2Aligned.fasta" "CRR2/"$1".CRR2Unique.fasta" "CRR2/"$1".CRR2Binary.phy" "CRR2"
python3 CRRAligner.py "CRR4/"$1".CRR4.fasta" "CRR4/"$1".CRR4Aligned.csv" "CRR4/"$1".CRR4Aligned.fasta" "CRR4/"$1".CRR4Unique.fasta" "CRR4/"$1".CRR4Binary.phy" "CRR4"


# Combine Alignments
echo "Combining Aligned CRISPR Cassettes"
python3 CombineCRRAligned.py "CRR1/"$1".CRR1Aligned.fasta" "CRR2/"$1".CRR2Aligned.fasta" "CRR4/"$1".CRR4Aligned.fasta" "AllCRR/"$1".AllCRRAligned.fasta"
python3 CombineCRRcsv.py "CRR1/"$1".CRR1Aligned.csv" "CRR2/"$1".CRR2Aligned.csv" "CRR4/"$1".CRR4Aligned.csv" "AllCRR/"$1".AllCRRAligned.csv"
python3 CombineCRRBinary.py "CRR1/"$1".CRR1Binary.phy" "CRR2/"$1".CRR2Binary.phy" "CRR4/"$1".CRR4Binary.phy" "AllCRR/"$1".AllCRRBinary.phy"

# Clean-up
rm -rf WriteFasta.py
rm -rf CombineCRRAligned.py
rm -rf CombineCRRcsv.py
rm -rf CombineCRRBinary.py
rm -rf EaCRRFinder.py
rm -rf CRRAligner.py