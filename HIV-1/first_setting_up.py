import os
from Bio.PDB import PDBParser, PDBIO
from vina import Vina

def perform_docking(receptor_pdbqt, ligand_pdbqt, output_file='out.pdbqt', center=(0,0,0), size=(20,20,20)):
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=size)
    v.dock(exhaustiveness=8, n_poses=5)
    v.write_poses(output_file, n_poses=5)
    print(f"Docking completed. Results saved to {output_file}")



perform_docking('receptor.pdbqt', 'ligand.pdbqt', 'ritonavir_docked.pdbqt', center=(10, 10, 10), size=(20,20,20))

