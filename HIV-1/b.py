import os
from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina

def prepare_receptor_openbabel(pdb_file, output_pdbqt='receptor.pdbqt'):
    """
    Convert a cleaned receptor PDB file to PDBQT using Open Babel.

    Args:
        pdb_file (str): Path to the cleaned receptor PDB file.
        output_pdbqt (str): Desired output PDBQT filename.
    """
    # Check if Open Babel is installed
    if not shutil.which("obabel"):
        raise EnvironmentError("Open Babel is not installed or not found in PATH.")
    
    cmd = f"obabel {pdb_file} -O {output_pdbqt} --partialcharge gasteiger"
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError("Error: Open Babel failed to convert receptor PDB to PDBQT.")
    else:
        print(f"Receptor PDBQT saved to {output_pdbqt}")

def prepare_ligand_openbabel(smiles, output_pdbqt='ligand.pdbqt'):
    """
    Convert a SMILES string to a PDBQT file using RDKit and Open Babel.

    Args:
        smiles (str): SMILES string of the ligand.
        output_pdbqt (str): Desired output PDBQT filename.
    """
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        raise RuntimeError("RDKit failed to embed molecule.")
    
    # Optimize geometry
    AllChem.UFFOptimizeMolecule(mol)
    
    # Save as PDB
    Chem.MolToPDBFile(mol, 'ligand.pdb')
    
    # Convert PDB to PDBQT using Open Babel
    if not shutil.which("obabel"):
        raise EnvironmentError("Open Babel is not installed or not found in PATH.")
    
    cmd = f"obabel ligand.pdb -O {output_pdbqt} --partialcharge gasteiger"
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError("Error: Open Babel failed to convert ligand PDB to PDBQT.")
    else:
        print(f"Ligand PDBQT saved to {output_pdbqt}")

def perform_docking(receptor_pdbqt, ligand_pdbqt, output_file='out.pdbqt', center=(0,0,0), size=(20,20,20)):
    """
    Perform docking using AutoDock Vina.

    Args:
        receptor_pdbqt (str): Path to the receptor PDBQT file.
        ligand_pdbqt (str): Path to the ligand PDBQT file.
        output_file (str): Output filename for docking poses.
        center (tuple): (x, y, z) coordinates of the docking box center.
        size (tuple): (x, y, z) dimensions of the docking box.
    """
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=size)
    v.dock(exhaustiveness=8, n_poses=5)
    v.write_poses(output_file, n_poses=5)
    print(f"Docking completed. Results saved to {output_file}")

def analyze_docking(output_file='out.pdbqt'):
    """
    Analyze docking results by extracting binding affinities.

    Args:
        output_file (str): Path to the docking output PDBQT file.
    """
    with open(output_file, 'r') as file:
        lines = file.readlines()
    
    scores = []
    for line in lines:
        if line.startswith('REMARK VINA RESULT:'):
            parts = line.strip().split()
            try:
                affinity = float(parts[3])
                scores.append(affinity)
            except (IndexError, ValueError):
                continue
    
    if not scores:
        print("No docking results found. Check if docking was performed correctly.")
        return
    
    # Sort scores (lower is better)
    scores_sorted = sorted(scores)
    print("Top binding affinities (kcal/mol):", scores_sorted[:5])

def main():
    # Step 1: Prepare Receptor
    receptor_pdb = './pdb_files/1hiv_with_h.pdb'  # Ensure you've manually added hydrogens
    prepare_receptor_openbabel(receptor_pdb, 'receptor.pdbqt')
    
    # Step 2: Prepare Ligand (Ritonavir)
    ritonavir_smiles = "CC(C)(C)OC1=CC=C(C=C1)C(=O)NC1CC(C(C(=O)N1CCCCC(C)C)=O)C(=O)NC(C)C"
    prepare_ligand_openbabel(ritonavir_smiles, 'ligand.pdbqt')
    
    # Step 3: Perform Docking
    docking_center = (10, 10, 10)  # Adjust based on active site coordinates
    docking_size = (20, 20, 20)     # Adjust to encompass the active site
    perform_docking('receptor.pdbqt', 'ligand.pdbqt', 'ritonavir_docked.pdbqt', center=docking_center, size=docking_size)
    
    # Step 4: Analyze Docking Results
    analyze_docking('ritonavir_docked.pdbqt')

if __name__ == "__main__":
    import shutil  # Needed for shutil.which
    main()

