from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

 
def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol

f = open('rxn.txt', 'rb')
while True:
    line = f.readline().strip('\r\n')
    if line:
        index = line.split('\t')[0]
        rxn = line.split('\t')[1]
        react = rxn.split('>>')[0]
        prod = rxn.split('>>')[1]
        R = Chem.MolFromSmiles(react)
        P = Chem.MolFromSmiles(prod)
        NewR = mol_with_atom_index(R)
        NewP = mol_with_atom_index(P)
        Draw.MolToMPL(NewR,size=(500, 500), fitImage = True)
        plt.savefig(index + '_R.jpg', dpi=72, bbox_inches = 'tight')
        Draw.MolToMPL(NewP,size=(500, 500), fitImage = True)
        plt.savefig(index + '_P.jpg', dpi=72, bbox_inches = 'tight')
        plt.close('all')
    else:
        break
##smi = "C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65"
##mol = Chem.MolFromSmiles(smi)
##new_mol = mol_with_atom_index(mol)
##
##Draw.MolToMPL(mol,size=(500, 500), fitImage = True)
##plt.savefig('mol.jpg',dpi=300, bbox_inches = 'tight')
##plt.show()
