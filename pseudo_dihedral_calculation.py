import numpy as np
import pandas as pd

'''User input file name'''
fileopen = open(input(), 'r')

'''Formating file for the operation'''
pdbfile = fileopen.readlines()

atoms: list[str] = []
group: list[str] = []
r_no: list[int] = []
X: list[float] = []
Y: list[float] = []
Z: list[float] = []
for lines in pdbfile[1:-1]:
    atoms.append(lines.split()[2])
    group.append(lines.split()[3])
    r_no.append(int(lines.split()[5]))
    X.append(float(lines.split()[6]))
    Y.append(float(lines.split()[7]))
    Z.append(float(lines.split()[8]))

pdb_dataframe = pd.DataFrame({'atoms': atoms,
                              'group': group,
                              'r_no': r_no,
                              'X': X,
                              'Y': Y,
                              'Z': Z},
                             columns=['atoms', 'group', 'r_no', 'X', 'Y', 'Z'])
pdb_dataframe.index += 1

'''Mass of each atom to calculate center of mass'''
atom_name: list[str] = ['C', 'N', 'O', 'P', 'H']
mass: list[float] = [12.011, 14.007, 15.999, 30.974, 1.008]

Mass = pd.DataFrame({'atom_name': atom_name,
                     'mass': mass},
                    columns=['atom_name', 'mass'])


def collecting_mass_axis_by_residue_number(residue_number: list[int]):
    X_i: list[float] = []
    Y_i: list[float] = []
    Z_i: list[float] = []
    m_i: list[float] = []
    for values in residue_number:
        in_atoms = pdb_dataframe.loc[pdb_dataframe["r_no"] == values]['atoms']
        X_i.append(pdb_dataframe.loc[pdb_dataframe["r_no"] == values]['X'])
        Y_i.append(pdb_dataframe.loc[pdb_dataframe["r_no"] == values]['Y'])
        Z_i.append(pdb_dataframe.loc[pdb_dataframe["r_no"] == values]['Z'])
        for atom_name_1 in in_atoms:
            m_1 = float(Mass.loc[Mass["atom_name"] == atom_name_1[0]]["mass"].iloc[0])
            m_i.append(m_1)

    return m_i, X_i, Y_i, Z_i


def collecting_mass_axis_by_atom_index(index_number: list[int]):
    X_i: list[float] = []
    Y_i: list[float] = []
    Z_i: list[float] = []
    m_i: list[float] = []
    for values in index_number:
        in_atoms = pdb_dataframe.loc[values]['atoms']

        X_i.append(pdb_dataframe.loc[values]['X'])
        Y_i.append(pdb_dataframe.loc[values]['Y'])
        Z_i.append(pdb_dataframe.loc[values]['Z'])
        m_1 = float(Mass.loc[Mass["atom_name"] == in_atoms[0]]["mass"].iloc[0])
        m_i.append(m_1)

    return m_i, X_i, Y_i, Z_i


def center_of_mass(no_type: str, numbers: list[int]):
    if no_type == 'r':
        center_of_mass = collecting_mass_axis_by_residue_number(numbers)
        CM_X = np.sum(np.multiply(center_of_mass[0], np.concatenate(center_of_mass[1]))) / np.sum(center_of_mass[0])
        CM_Y = np.sum(np.multiply(center_of_mass[0], np.concatenate(center_of_mass[2]))) / np.sum(center_of_mass[0])
        CM_Z = np.sum(np.multiply(center_of_mass[0], np.concatenate(center_of_mass[3]))) / np.sum(center_of_mass[0])

    if no_type == 'i':
        center_of_mass = collecting_mass_axis_by_atom_index(numbers)
        CM_X = np.sum(np.multiply(center_of_mass[0], center_of_mass[1])) / np.sum(center_of_mass[0])
        CM_Y = np.sum(np.multiply(center_of_mass[0], center_of_mass[2])) / np.sum(center_of_mass[0])
        CM_Z = np.sum(np.multiply(center_of_mass[0], center_of_mass[3])) / np.sum(center_of_mass[0])

    return CM_X, CM_Y, CM_Z


'''Input from the user for four dihedral points'''

P1 = center_of_mass('i', [89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 292, 293, 294, 295, 296, 297, 298, 299, 300, 50, 51, 52, 53, 54, 55, 56, 57, 58, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342])
P2 = center_of_mass('i', [71, 78, 78, 80, 81])
P3 = center_of_mass('i', [74, 75, 76, 77, 47])
P4 = center_of_mass('i', [59, 60, 61, 62, 63, 64, 65, 66])

'''calculating dihedral angle'''
b0 = np.subtract(P2, P1)
b1 = np.subtract(P3, P2)
b2 = np.subtract(P4, P3)

b0_C_b1 = np.cross(b0, b1)
b1_C_b2 = np.cross(b1, b2)

cos_theta = np.dot(b0_C_b1, b1_C_b2) / (np.linalg.norm(b0_C_b1)* np.linalg.norm(b1_C_b2))

sin_theta = np.dot(np.cross(b0_C_b1, b1_C_b2), b1) / (np.linalg.norm(b0_C_b1) *
                                                      np.linalg.norm(b1_C_b2) * np.linalg.norm(b1))

theta = np.arctan2(sin_theta, cos_theta) * 180/np.pi
print("{:.2f}".format(theta))
