import os
from rdkit import Chem
from rdkit.Chem import AllChem

def prepare_ligand(smiles, output_file='ligand.pdbqt'):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    
    # Save as PDB
    Chem.MolToPDBFile(mol, 'ligand.pdb')
    
    #  export PYTHONPATH="/home/deep/MGLTools-1.5.7/MGLToolsPckgs:$PYTHONPATH"

    # Convert to PDBQT (requires 'prepare_ligand' script from MGLTools in PATH)
    ret = os.system(f"python2 /home/deep/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand.py -l ligand.pdb -o {output_file}")
    if ret != 0:
        print("Error running 'prepare_ligand' script. Check MGLTools/AutoDock Tools installation.")
    else:
        print(f"Ligand prepared and saved to {output_file}")

if __name__ == "__main__":
    ritonavir_smiles = "CC(C)(C)OC1=CC=C(C=C1)C(=O)NC1CC(C(C(=O)N1CCCCC(C)C)=O)C(=O)NC(C)C"
    prepare_ligand(ritonavir_smiles)

