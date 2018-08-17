# RNA3DCNN: Local and Global Quality Assessments of RNA 3D Structures Using 3D Deep Convolutional Neural Networks

In our work, the evaluation of an RNA structure is divided into the assessment of each nucleotide. The nucleotide to be evaluated and its surrounding atoms are treated as a 3D image that is fed into the 3D CNNs, and the output is an RMSD-based nucleotide unfitness score characterizing how poorly the nucleotide fits into its surrounding environment.

## Prerequisites
1. Anaconda https://www.anaconda.com/download/
2. Keras https://keras.io/#installation
3. biopython https://biopython.org/wiki/Download

## Usage

1. To assess an RNA, use flag "-pn pdbname" 
2. To assess a list of RNAs, use flag "-pl pdblist"
3. To choose different 3DCNN models, use flag "-model RNA3DCNN_MD.hdf5" or "-model RNA3DCNN_MDMC.hdf5"
4. To print scores of each nucleotide and total scores, use flag "-local 1"
5. To print only total scores, use flag "-local 0"

For example:<br />
python Main.py -pl pdblist -model RNA3DCNN_MD.hdf5 -local 1<br />

The command above will assess RNAs in pdblist using RNA3DCNN_MD CNN model, and print scores of each nucleotide and the total scores of each RNA.<br />
It is very easy for you to read the Main.py and change the form of output.
