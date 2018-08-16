from Bio.PDB.Vector import Vector
import numpy as np

""" pixelate each residue and its neighbors """

""" atom name """
ATOMS_NAME = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'C3\'', 'O3\'',
              'O4\'', 'C1\'', 'C2\'', 'O2\'', 'N9', 'C8', 'N7', 'C5', 'C6',
              'O6', 'N1', 'C2', 'N2', 'N3', 'C4', 'N6', 'O4', 'O2', 'N4']

""" atom mass """
MP = 30.9738
MO = 15.9994
MC = 12.0107
MN = 14.0067

""" atom charge """
QP = 1.1662
QOP1 = -0.776
QOP2 = -0.776
QO5S_1 = -0.4989
QO5S_2 = -0.6223
QC5S = 0.0558
QC4S = 0.1065
QC3S = 0.2022
QO3S_1 = -0.5246
QO3S_2 = -0.6541
QO4S = -0.3548
QC1S_A = 0.0394
QC1S_G = 0.0191
QC1S_U = 0.0674
QC1S_C = 0.0066
QC2S = 0.067
QO2S = -0.6139
QN9_A = -0.0251
QN9_G = 0.0492
QC8_A = 0.2006
QC8_G = 0.1374
QN7_A = -0.6073
QN7_G = -0.5709
QC5_A = 0.0515
QC5_G = 0.1744
QC5_U = -0.3635
QC5_C = -0.5215
QC6_A = 0.7009
QC6_G = 0.477
QC6_U = -0.1126
QC6_C = 0.0053
QO6 = -0.5597
QN1_A = -0.7615
QN1_G = -0.4787
QN1_U = 0.0418
QN1_C = -0.0484
QC2_A = 0.5875
QC2_G = 0.7657
QC2_U = 0.4687
QC2_C = 0.7538
QN2 = -0.9672
QN3_A = -0.6997
QN3_G = -0.6323
QN3_U = -0.3549
QN3_C = -0.7584
QC4_A = 0.3053
QC4_G = 0.1222
QC4_U = 0.5952
QC4_C = 0.8185
QN6 = -0.9019
QO4 = -0.5761
QO2_U = -0.5477
QO2_C = -0.6252
QN4 = -0.953

DICT_MASS = {'P': MP, 'OP1': MO, 'OP2': MO, 'O5\'': MO, 'C5\'': MC, 'C4\'': MC,
             'C3\'': MC, 'O3\'': MO, 'O4\'': MO, 'C1\'': MC, 'C2\'': MC,
             'O2\'': MO, 'N9': MN, 'C8': MC, 'N7': MN, 'C5': MC, 'C6': MC,
             'O6': MO, 'N1': MN, 'C2': MC, 'N2': MN, 'N3': MN, 'C4': MC,
             'N6': MN, 'O4': MO, 'O2': MO, 'N4': MN}

DICT_CHARGE = {'P': QP, 'OP1': QOP1, 'OP2': QOP2, 'C5\'': QC5S, 'C4\'': QC4S,
               'C3\'': QC3S, 'O4\'': QO4S, 'C2\'': QC2S, 'O2\'': QO2S,
               'O6': QO6, 'N2': QN2, 'N6': QN6, 'O4': QO4, 'N4': QN4}

DICT_CHARGE_O5S = {'nonhead': QO5S_1, 'head': QO5S_2}

DICT_CHARGE_O3S = {'nontail': QO3S_1, 'tail': QO3S_2}

DICT_CHARGE_C1S = {'A': QC1S_A, 'G': QC1S_G, 'U': QC1S_U, 'C': QC1S_C}

DICT_CHARGE_N9 = {'A': QN9_A, 'G': QN9_G}

DICT_CHARGE_C8 = {'A': QC8_A, 'G': QC8_G}

DICT_CHARGE_N7 = {'A': QN7_A, 'G': QN7_G}

DICT_CHARGE_C5 = {'A': QC5_A, 'G': QC5_G, 'U': QC5_U, 'C': QC5_C}

DICT_CHARGE_C6 = {'A': QC6_A, 'G': QC6_G, 'U': QC6_U, 'C': QC6_C}

DICT_CHARGE_N1 = {'A': QN1_A, 'G': QN1_G, 'U': QN1_U, 'C': QN1_C}

DICT_CHARGE_C2 = {'A': QC2_A, 'G': QC2_G, 'U': QC2_U, 'C': QC2_C}

DICT_CHARGE_N3 = {'A': QN3_A, 'G': QN3_G, 'U': QN3_U, 'C': QN3_C}

DICT_CHARGE_C4 = {'A': QC4_A, 'G': QC4_G, 'U': QC4_U, 'C': QC4_C}

DICT_CHARGE_O2 = {'U': QO2_U, 'C': QO2_C}

""" 3D-image size """

BOX_SIDE_LENGTH = 32.0
HALF_BOX_SIDE_LENGTH = BOX_SIDE_LENGTH/2.0
TRIAL_HALF_BOX_SIDE_LENGTH = HALF_BOX_SIDE_LENGTH + 5.0
BIN_WIDTH = 1.0
HALF_BIN_WIDTH = BIN_WIDTH/2.0
NBINS = int(BOX_SIDE_LENGTH/BIN_WIDTH)


def check_if_lacking_atoms(res, atomname):
    """check if lacking certain atoms in residue"""
    if not res.has_id(atomname):
        print ('There is no atom ' + atomname +\
                ' in residue ' + str(res.get_id()[1]) +\
                res.get_resname()[2] +\
                ' in chain ', res.get_full_id()[2] +\
                ' in PDB ' + res.get_full_id()[0] + '.')
        return True
    else:
        return False


def calc_local_reference(res, local_ref):
    """calculate the local reference system based on the centered residue.
    local_ref is a list of ori, vx, vy, vz. """

    if check_if_lacking_atoms(res, 'C1\''):
        return False
    local_ref[0] = res['C1\''].get_vector()

    res_name = res.get_resname()[2]
    if res_name == 'A' or res_name == 'G':
        if check_if_lacking_atoms(res, 'N9'):
            return False
        local_ref[1] = res['N9'].get_vector() - local_ref[0]
    elif res_name == 'U' or res_name == 'C':
        if check_if_lacking_atoms(res, 'N1'):
            return False
        local_ref[1] = res['N1'].get_vector() - local_ref[0]
    else:
        print ("non-canonical residue name", res.get_resname())
        exit()

    if check_if_lacking_atoms(res, 'O5\''):
        return False
    if check_if_lacking_atoms(res, 'C5\''):
        return False
    r_o5s = res['O5\''].get_vector()
    r_c5s = res['C5\''].get_vector()
    #local_ref[2] = (r_o5s + r_c5s)/2 - local_ref[0]
    local_ref[2][0] = (r_o5s[0] + r_c5s[0])/2.0 - local_ref[0][0]
    local_ref[2][1] = (r_o5s[1] + r_c5s[1])/2.0 - local_ref[0][1]
    local_ref[2][2] = (r_o5s[2] + r_c5s[2])/2.0 - local_ref[0][2]
    local_ref[3] = local_ref[1] ** local_ref[2]
    local_ref[3].normalize()
    local_ref[1].normalize()
    local_ref[2] = local_ref[3] ** local_ref[1]
    local_ref[2].normalize()

    return True


def lattice_1d_point(dis, int_coord, prob):
    """ lattice a 1d point """

    int_coord[0] = int((dis + HALF_BIN_WIDTH)/BIN_WIDTH) - 1
    prob[0] = ((int_coord[0] + 1) * BIN_WIDTH + HALF_BIN_WIDTH - dis)\
        / BIN_WIDTH
    int_coord[1] = int_coord[0] + 1
    prob[1] = 1.0 - prob[0]

    if int_coord[0] < -1 or int_coord[1] > NBINS:
        print ('Atoms overflow the box.')
        exit()


def lattice_3d_point(vec_3d, int_coord_x, int_coord_y, int_coord_z, prob_x,
                     prob_y, prob_z):
    """ lattice a 3d point """

    lattice_1d_point(vec_3d[0], int_coord_x, prob_x)
    lattice_1d_point(vec_3d[1], int_coord_y, prob_y)
    lattice_1d_point(vec_3d[2], int_coord_z, prob_z)


def is_neighbor_residue(local_ref, residue):
    """ if residue is a neighbor:
            return True
        else:
            return False
    """
    if check_if_lacking_atoms(residue, 'C1\''):
        return True

    vec_r = residue['C1\''].get_vector() - local_ref[0]
    if abs(vec_r * local_ref[1]) > TRIAL_HALF_BOX_SIDE_LENGTH:
        return False
    if abs(vec_r * local_ref[2]) > TRIAL_HALF_BOX_SIDE_LENGTH:
        return False
    if abs(vec_r * local_ref[3]) > TRIAL_HALF_BOX_SIDE_LENGTH:
        return False

    return True


def is_neighbor_atom(local_ref, atom, local_coordinate):
    """ if atom is a neighbor:
            update local_coordinate
            return True
        else:
            return False
    """
    vec_r = atom.get_vector() - local_ref[0]
    dis_x = vec_r * local_ref[1]
    if abs(dis_x) >= HALF_BOX_SIDE_LENGTH:
        return False
    dis_y = vec_r * local_ref[2]
    if abs(dis_y) >= HALF_BOX_SIDE_LENGTH:
        return False
    dis_z = vec_r * local_ref[3]
    if abs(dis_z) >= HALF_BOX_SIDE_LENGTH:
        return False

    local_coordinate[0] = dis_x + HALF_BOX_SIDE_LENGTH
    local_coordinate[1] = dis_y + HALF_BOX_SIDE_LENGTH
    local_coordinate[2] = dis_z + HALF_BOX_SIDE_LENGTH

    return True

def is_head_residue(residue_list, res_index, residue):
    """ if residue is the head residue:
            return True
        else:
            return False
    """
    if res_index == 0:
        return True
    elif residue.get_full_id()[2] !=\
            residue_list[res_index - 1].get_full_id()[2]:
        return True
    else:
        return False


def is_tail_residue(residue_list, res_index, residue):
    """ if residue is the tail residue:
            return True
        else:
            return False
    """
    if res_index == len(residue_list) - 1:
        return True
    elif residue.get_full_id()[2] !=\
            residue_list[res_index + 1].get_full_id()[2]:
        return True
    else:
        return False


def pixelate_atoms_in_box(model, pixels):
    """ pixelate atoms in residue-centered box by BINWIDTH """

    residues = list(model.get_residues())

    ori = Vector(0, 0, 0)
    vec_x = Vector(0, 0, 0)
    vec_y = Vector(0, 0, 0)
    vec_z = Vector(0, 0, 0)

    int_coord_x = [0, 0]
    int_coord_y = [0, 0]
    int_coord_z = [0, 0]
    prob_x = [0, 0]
    prob_y = [0, 0]
    prob_z = [0, 0]

    res_index1 = 0
    for residue1 in residues:
        local_ref = [ori, vec_x, vec_y, vec_z]
        if not calc_local_reference(residue1, local_ref):
            pixels = np.delete(pixels, len(pixels)-1, 0)
            continue

        res_index2 = 0
        for residue2 in residues:
            if not is_neighbor_residue(local_ref, residue2):
                res_index2 += 1
                continue

            res_name = residue2.get_resname()[2]

            for atom in residue2:
                atom_name = atom.get_name()

                if atom_name not in ATOMS_NAME:
                    continue
                if atom_name == 'P' or atom_name == 'OP1' or atom_name == 'OP2':
                    if is_head_residue(residues, res_index2, residue2):
                        continue

                local_coordinate = Vector(0, 0, 0)
                if not is_neighbor_atom(local_ref, atom, local_coordinate):
                    continue

                lattice_3d_point(local_coordinate, int_coord_x, int_coord_y,
                                 int_coord_z, prob_x, prob_y, prob_z)

                atom_mass = DICT_MASS[atom_name]

                if atom_name in DICT_CHARGE:
                    atom_charge = DICT_CHARGE[atom_name]
                elif atom_name == 'O5\'':
                    if is_head_residue(residues, res_index2, residue2):
                        atom_charge = DICT_CHARGE_O5S['head']
                    else:
                        atom_charge = DICT_CHARGE_O5S['nonhead']
                elif atom_name == 'O3\'':
                    if is_tail_residue(residues, res_index2, residue2):
                        atom_charge = DICT_CHARGE_O3S['tail']
                    else:
                        atom_charge = DICT_CHARGE_O3S['nontail']
                elif atom_name == 'C1\'':
                    atom_charge = DICT_CHARGE_C1S[res_name]
                elif atom_name == 'N9':
                    atom_charge = DICT_CHARGE_N9[res_name]
                elif atom_name == 'C8':
                    atom_charge = DICT_CHARGE_C8[res_name]
                elif atom_name == 'N7':
                    atom_charge = DICT_CHARGE_N7[res_name]
                elif atom_name == 'C5':
                    atom_charge = DICT_CHARGE_C5[res_name]
                elif atom_name == 'C6':
                    atom_charge = DICT_CHARGE_C6[res_name]
                elif atom_name == 'N1':
                    atom_charge = DICT_CHARGE_N1[res_name]
                elif atom_name == 'C2':
                    atom_charge = DICT_CHARGE_C2[res_name]
                elif atom_name == 'N3':
                    atom_charge = DICT_CHARGE_N3[res_name]
                elif atom_name == 'C4':
                    atom_charge = DICT_CHARGE_C4[res_name]
                elif atom_name == 'O2':
                    atom_charge = DICT_CHARGE_O2[res_name]
                else:
                    continue
                
                for i in range(2):
                    if int_coord_x[i] < 0 or int_coord_x[i] >= NBINS:
                        continue
                    for j in range(2):
                        if int_coord_y[j] < 0 or int_coord_y[j] >= NBINS:
                            continue
                        for k in range(2):
                            if int_coord_z[k] < 0 or int_coord_z[k] >= NBINS:
                                continue
                            (pixels[res_index1][0][int_coord_x[i]]
                             [int_coord_y[j]][int_coord_z[k]]) +=\
                                prob_x[i] * prob_y[j] * prob_z[k]
                            (pixels[res_index1][1][int_coord_x[i]]
                             [int_coord_y[j]][int_coord_z[k]]) +=\
                                prob_x[i] * prob_y[j] * prob_z[k] * atom_mass
                            (pixels[res_index1][2][int_coord_x[i]]
                             [int_coord_y[j]][int_coord_z[k]]) +=\
                                prob_x[i] * prob_y[j] * prob_z[k] * atom_charge
                            """
                            if res_index1 == 0 and int_coord_x[i] == 23 and int_coord_y[j] == 6 and int_coord_z[k] == 3:
                                print atom_name, atom_mass, atom_charge, prob_x[i], prob_y[j], prob_z[k], residue2.get_full_id()
                            """
            res_index2 += 1
        res_index1 += 1
    
    return pixels
