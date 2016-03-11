#!/usr/bin/python

'''
author: Shantanu S. Bhattacharyya

This script requires Biopython to be installed beforehand.

Purpose : Detect residues that form the interface between specified chain(s)

Usage : python interface_chains.py --pdb <input pdb> --chain1 <list of chains that form the query> --chain2 <list of chains in which
            the surrounding residues are to be found>

Note : For detecting residues interfacing specific residues rather than specific chains, use interface_residues script

'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", required=True)
parser.add_argument("--chain1", nargs='+', required=True)
parser.add_argument("--chain2", nargs='+', required=True)

args = parser.parse_args()


from Bio.PDB import NeighborSearch, PDBParser, Selection

structure = PDBParser().get_structure('X', args.pdb)

center_atoms = []
pymol_command = ""

all_atom_list = [atom for atom in structure.get_atoms() if atom.name == 'CA' ]

for k in args.chain1 :
    chain_atoms = [atom for atom in structure[0][k].get_atoms() if atom.name == 'CA' ]

    center_atoms += chain_atoms

#atom_list = [x for x in all_atom_list if x not in center_atoms]

for j in args.chain2 :
    atom_list = [atom for atom in structure[0][j].get_atoms() if atom.name == 'CA' ]
    ns = NeighborSearch(atom_list)

    nearby_residues = {res for center_atom in center_atoms
                        for res in ns.search(center_atom.coord, 8.5, 'R')}

    print "\nNeighbor residues in chain ", j, ": \n"
    print sorted(res.id[1] for res in nearby_residues)

    pymol_command = "show spheres, chain " + j + " and resi "

    for m in sorted(res.id[1] for res in nearby_residues):
        pymol_command = pymol_command + str(m) + "+"

    print pymol_command[:-1] + " and name CA \n" 