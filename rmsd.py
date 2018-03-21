#!/usr/bin/env python
"""
Calculate the rooted mean squared distance between two protein structures.

Usage:


$ python rmsd.py file1.pdb file2.pdb >> output.txt
"""
from __future__ import print_function

__author__ = "Rodrigo Dorantes-Gilardi"
__email__ = "rodgdor@gmail.com"


import sys
import Bio.PDB
import numpy as np


def rmsd(prot1, prot2):
    """
    Return RMSD between two pdb files.

    Water molecules and het-atoms are removed from the structures.

    input
    -----
    prot1, prot2: str. The path to the two pdb files.
    
    output
    ------
    float. RMSD between prot1 and prot2.
    """
    RMSD = []
    parser = Bio.PDB.PDBParser(QUIET=True)
    try:
        prot1 = parser.get_structure(prot1, prot1)[0]
    except FileNotFoundError:
        print("Cannot find the pdb file {}".format(prot1))

    try:
        prot2 = parser.get_structure(prot2, prot2)[0]
    except FileNotFoundError:
        print("Cannot find the pdb file {}".format(prot2))

    # Remove water molecules
    def remove_water(prot):
        for residue in prot.get_residues():
            if residue.id[0] != ' ':
                residue.parent.detach_child(residue.id)
        pass

    remove_water(prot1)
    remove_water(prot2)

    for atom in prot1.get_atoms():
        chain = atom.parent.parent
        res = atom.parent
        name = atom.name
        try:
            atom2 = prot2[chain.id][res.id[1]][name]
            RMSD.append((atom - atom2) ** 2)
        except KeyError:
            print(
                ("Atom {0} in chain {1} and position {2}"
                " does not seem to be in both "
                "proteins").format(name, chain.id, res.id[1]))

    return np.sqrt(np.mean(RMSD))
    


if __name__ == '__main__':
    args = sys.argv
    if args[1][-3:] != "pdb" or args[2][-3:] != "pdb":
        raise NameError("Both files must end with '.pdb'")

    _rmsd = rmsd(args[1], args[2])
    print("RMSD value (angstroms): {}\n".format(_rmsd))
