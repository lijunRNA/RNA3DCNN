import sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from PixelateResidue import *
from ModifyName import *
from keras.models import load_model

"""local and global assessment of RNA"""

def load_CNN_model(model_name):
    """ load CNN model for scoring each residue."""

    model = load_model(model_name)
    return model


def preprocess_input(input_pixels):
    """ preprocess input of CNN """

    input_pixels[:, 1] /= 31.0
    input_pixels[:, 2] += 1.0
    input_pixels[:, 2] /= 2.2
    return input_pixels


if __name__ == '__main__':
    if len(sys.argv) != 7:
        print ('Error! Usage: python main.py -pn pdbname -model cnn_model.hdf5 -local 0/1 or python main.py -pl pdblist -model cnn_model.hdf5 -local 0/1')
        exit()
    if sys.argv[1] == '-pn':
        pdblist = sys.argv[2].split()
    elif sys.argv[1] == '-pl':
        pdblist = open(sys.argv[2]).read().splitlines()
    else:
        print ('Error! Usage: python main.py -pn pdbname -model cnn_model.hdf5 or python main.py -pl pdblist -model cnn_model.hdf5 -local 0/1')
        exit()

    if sys.argv[5] == '-local':
        if sys.argv[6] == '0':
            local_assessment = False
        else:
            local_assessment = True
    else:
        print ('Error! Usage: python main.py -pn pdbname -model cnn_model.hdf5 or python main.py -pl pdblist -model cnn_model.hdf5 -local 0/1')

    cnn_model_name = sys.argv[4]
    cnn_model = load_CNN_model(cnn_model_name)
    cnn_model.summary()

    for RNA in pdblist:
        p = PDBParser(QUIET=True)
        s = p.get_structure(RNA, RNA)
        """ only consider the first model in PDB file """
        model = s[0]
        residues = list(model.get_residues())
        length = len(residues)
        modify_residue_atom_name(residues)
        pixels = np.zeros((length, 3, NBINS, NBINS, NBINS))
        pixels = pixelate_atoms_in_box(model, pixels)

        pixels = preprocess_input(pixels)
        score_residue = cnn_model.predict(pixels)

        if local_assessment:
            print ("Scores for each nucleotide in " + RNA + ":")
            print (score_residue)
            print ("Total score for " + RNA + " is ", np.sum(score_residue))
        else:
            print ("Total score for " + RNA + " is ", np.sum(score_residue))
