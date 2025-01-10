import os
from Bio.PDB import PDBList

def download_pdb(pdb_id, save_dir='pdb_files'):
    pdbl = PDBList()
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    pdbl.retrieve_pdb_file(pdb_id, pdir=save_dir, file_format='pdb')
    print(f"PDB file for {pdb_id} downloaded to {save_dir}/pdb{pdb_id.lower()}.ent")

# Try a different PDB ID for HIV-1 protease, e.g., 5HVP
download_pdb('5HVP')

