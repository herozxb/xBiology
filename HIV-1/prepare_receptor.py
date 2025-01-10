import os
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB import Select

STANDARD_RESIDUES = [
    'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
]

class ResidueSelect(Select):
    def accept_residue(self, residue):
        if residue.id[0] == ' ' and residue.get_resname() in STANDARD_RESIDUES:
            return True
        return False

# python2 ~/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ./pdb_files/1hiv.pdb -o receptor_vina.pdbqt -A checkhydrogens -U nphs_lps_waters


def prepare_receptor(pdb_file, output_file='receptor.pdbqt'):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', pdb_file)
    
    # Save a cleaned receptor
    io = PDBIO()
    io.set_structure(structure)
    io.save('clean_receptor.pdb', select=ResidueSelect())
    
    # Convert to PDBQT (requires 'prepare_receptor' script in PATH)
    # or call 'prepare_receptor4.py' directly
    cmd = f"python2 /home/deep/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor.py -r 1hiv_with_h.pdb -o {output_file}"
    ret = os.system(cmd)
    if ret != 0:
        print("Error: 'prepare_receptor' script did not run successfully!")
    else:
        print(f"Clean receptor saved as clean_receptor.pdb; PDBQT: {output_file}")

path = "./pdb_files/1hiv.pdb"
print(path)

prepare_receptor(path)

