CANONICAL_RES_NAME = {'  A', '  G', '  U', '  C'}

DICT_RES_NAME = {' RA': '  A', ' RG': '  G', ' RU': '  U', ' RC': '  C',
                 '2MG': '  G', 'H2U': '  U', 'M2G': '  G', 'OMG': '  G',
                 ' YG': '  G', 'PSU': '  U', '5MC': '  C', '7MG': '  G',
                 '5MU': '  U', '1MA': '  A', 'OMC': '  C', '  I': '  G',
                 '1MG': '  G', 'GDP': '  G', 'A23': '  A', '5BU': '  U',
                 '5IC': '  C', 'CB2': '  C', 'GTP': '  G', 'DHU': '  U',
                 'AET': '  A', 'G7M': '  G', '4SU': '  U', 'CCC': '  C',
                 'S4U': '  U', '  T': '  U', 'FHU': '  U', 'AVC': '  A',
                 'OMU': '  U', 'UR3': '  U', 'T6A': '  A', 'RIA': '  A',
                 'PGP': '  G', 'BRU': '  U', 'U34': '  U', 'YYG': '  G',
                 'CBR': '  C', 'A2M': '  A', 'BGM': '  G', 'UMS': '  U',
                 'CSL': '  C', ' IU': '  U', 'UD5': '  U', 'S4C': '  C',
                 'FMU': '  U', '5FU': '  U', ' DU': '  U', 'XUG': '  G',
                 'TM2': '  U', 'G46': '  G', '1SC': '  C', 'CFL': '  C',
                 'UFT': '  U', 'SUR': '  U', 'MTU': '  G', '6FC': '  C',
                 ' CH': '  C', 'U8U': '  U', 'RUS': '  U', ' IG': '  G',
                 ' IC': '  C', '6MZ': '  A', 'CM0': '  U', 'MIA': '  A',
                 ' 0C': '  C', ' 0U': '  U', ' 0G': '  G', ' DG': '  G',
                 'AP7': '  A', 'LCA': '  A', '10C': '  C', 'SSU': '  U',
                 'CBV': '  C', 'RA5': '  A', 'RG5': '  G', 'RU5': '  U',
                 'RC5': '  C', 'RA3': '  A', 'RG3': '  G', 'RU3': '  U',
                 'RC3': '  C', 'PPU': '  A', '  N': '  C', ' rA': '  A',
                 ' rG': '  G', ' rU': '  U', ' rC': '  C', '6IA': '  A'}


def modify_residue_atom_name(residues):
    """ change '*' in atom name into '\''
        change 'O1P' into 'OP1'
        change 'O2P' into 'OP2'
    """
    for residue in residues:
        if residue.resname in DICT_RES_NAME:
            residue.resname = DICT_RES_NAME[residue.resname]
        if residue.resname not in CANONICAL_RES_NAME:
            print (residue.resname, 'is not a canonical residue name.')
            exit(0)
        for atom in residue:
            if '*' in atom.name:
                del atom.parent.child_dict[atom.id]
                atom.id = atom.id.replace('*', '\'')
                atom.name = atom.name.replace('*', '\'')
                atom.parent.child_dict[atom.id] = atom
            if 'O1P' in atom.name:
                del atom.parent.child_dict[atom.id]
                atom.id = atom.id.replace('O1P', 'OP1')
                atom.name = atom.name.replace('O1P', 'OP1')
                atom.parent.child_dict[atom.id] = atom
            if 'O2P' in atom.name:
                del atom.parent.child_dict[atom.id]
                atom.id = atom.id.replace('O2P', 'OP2')
                atom.name = atom.name.replace('O2P', 'OP2')
                atom.parent.child_dict[atom.id] = atom
