import os
import torch
import numpy as np
import csv
import hashlib
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser

# ── Constants ──────────────────────────────────────────────────────────────────

# 14 heavy atoms per residue (standard RoseTTAFold atom ordering)
ATOM14_NAMES = [
    'N', 'CA', 'C', 'O',           # backbone
    'CB',                            # beta carbon
    'CG', 'CG1', 'CG2',
    'CD', 'CD1', 'CD2',
    'CE', 'CE1', 'CE2',
    'NZ'
]

AA3_TO_1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','SEP':'S','TPO':'T','CSO':'C','PTR':'Y',  # common mods
}

# Per-residue atom14 ordering (which atoms go in which slot)
# This follows the standard used in RoseTTAFold
RESIDUE_ATOM14 = {
    'ALA': ['N','CA','C','O','CB','','','','','','','','',''],
    'ARG': ['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2','','',''],
    'ASN': ['N','CA','C','O','CB','CG','OD1','ND2','','','','','',''],
    'ASP': ['N','CA','C','O','CB','CG','OD1','OD2','','','','','',''],
    'CYS': ['N','CA','C','O','CB','SG','','','','','','','',''],
    'GLN': ['N','CA','C','O','CB','CG','CD','OE1','NE2','','','','',''],
    'GLU': ['N','CA','C','O','CB','CG','CD','OE1','OE2','','','','',''],
    'GLY': ['N','CA','C','O','','','','','','','','','',''],
    'HIS': ['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2','','','',''],
    'ILE': ['N','CA','C','O','CB','CG1','CG2','CD1','','','','','',''],
    'LEU': ['N','CA','C','O','CB','CG','CD1','CD2','','','','','',''],
    'LYS': ['N','CA','C','O','CB','CG','CD','CE','NZ','','','','',''],
    'MET': ['N','CA','C','O','CB','CG','SD','CE','','','','','',''],
    'PHE': ['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','','',''],
    'PRO': ['N','CA','C','O','CB','CG','CD','','','','','','',''],
    'SER': ['N','CA','C','O','CB','OG','','','','','','','',''],
    'THR': ['N','CA','C','O','CB','OG1','CG2','','','','','','',''],
    'TRP': ['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'],
    'TYR': ['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH','',''],
    'VAL': ['N','CA','C','O','CB','CG1','CG2','','','','','','',''],
}
# fallback for modified residues
DEFAULT_ATOM14 = ['N','CA','C','O','CB','','','','','','','','','']


# ── Atom14 extraction ──────────────────────────────────────────────────────────

def residue_to_atom14(residue):
    """
    Returns:
        xyz  : (14, 3) float32
        mask : (14,)   bool
        bfac : (14,)   float32
        occ  : (14,)   float32
    """
    resname = residue.resname.strip()
    atom_order = RESIDUE_ATOM14.get(resname, DEFAULT_ATOM14)

    xyz  = np.zeros((14, 3), dtype=np.float32)
    mask = np.zeros(14,      dtype=bool)
    bfac = np.zeros(14,      dtype=np.float32)
    occ  = np.zeros(14,      dtype=np.float32)

    for i, atom_name in enumerate(atom_order):
        if not atom_name:
            continue
        if atom_name in residue:
            atom = residue[atom_name]
            xyz[i]  = atom.coord
            mask[i] = True
            bfac[i] = atom.bfactor
            occ[i]  = atom.occupancy
        # if atom missing, xyz stays 0, mask stays False

    return xyz, mask, bfac, occ


# ── Chain-level .pt ────────────────────────────────────────────────────────────

def extract_chain_pt(chain):
    """
    Returns dict with:
        seq  - str of 1-letter AA codes
        xyz  - FloatTensor [L,14,3]
        mask - BoolTensor  [L,14]
        bfac - FloatTensor [L,14]
        occ  - FloatTensor [L,14]
    """
    seq_chars = []
    all_xyz   = []
    all_mask  = []
    all_bfac  = []
    all_occ   = []

    for residue in chain.get_residues():
        # skip HETATM (waters, ligands) — keep only standard residues
        if residue.id[0] != ' ':
            continue
        resname = residue.resname.strip()
        if 'CA' not in residue:
            continue

        aa1 = AA3_TO_1.get(resname, 'X')
        seq_chars.append(aa1)

        xyz, mask, bfac, occ = residue_to_atom14(residue)
        all_xyz.append(xyz)
        all_mask.append(mask)
        all_bfac.append(bfac)
        all_occ.append(occ)

    if not seq_chars:
        return None

    return {
        'seq':  ''.join(seq_chars),
        'xyz':  torch.tensor(np.stack(all_xyz),  dtype=torch.float32),
        'mask': torch.tensor(np.stack(all_mask), dtype=torch.bool),
        'bfac': torch.tensor(np.stack(all_bfac), dtype=torch.float32),
        'occ':  torch.tensor(np.stack(all_occ),  dtype=torch.float32),
    }


# ── Biological assembly parsing ────────────────────────────────────────────────

def parse_remark350(pdb_path):
    """
    Parse REMARK 350 from a PDB file to extract biological assembly info.
    Returns lists: asmb_ids, asmb_details, asmb_method, asmb_chains, asmb_xforms
    """
    asmb_ids     = []
    asmb_details = []
    asmb_method  = []
    asmb_chains  = []
    asmb_xforms  = []

    current_id      = None
    current_chains  = []
    current_xforms  = []
    current_detail  = ''
    current_method  = ''
    current_matrix_rows = []
    current_matrix_chains = []

    def flush_assembly():
        if current_id is not None:
            asmb_ids.append(current_id)
            asmb_details.append(current_detail)
            asmb_method.append(current_method)
            asmb_chains.append(','.join(current_chains))
            # stack xforms: [n,4,4]
            if current_xforms:
                asmb_xforms.append(torch.tensor(current_xforms, dtype=torch.float32))
            else:
                asmb_xforms.append(torch.eye(4).unsqueeze(0))

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('REMARK 350'):
                continue
            content = line[11:].strip()

            if content.startswith('BIOMOLECULE:'):
                flush_assembly()
                current_id      = content.split('BIOMOLECULE:')[1].strip()
                current_chains  = []
                current_xforms  = []
                current_detail  = ''
                current_method  = ''
                current_matrix_rows = []
                current_matrix_chains = []

            elif 'AUTHOR DETERMINED BIOLOGICAL UNIT' in content:
                current_detail = 'author'
            elif 'SOFTWARE DETERMINED QUATERNARY STRUCTURE' in content:
                current_detail = 'software'
            elif 'SOFTWARE USED:' in content:
                current_method = content.split('SOFTWARE USED:')[1].strip()

            elif 'APPLY THE FOLLOWING TO CHAINS:' in content:
                chain_str = content.split('CHAINS:')[1].strip().rstrip(',')
                current_matrix_chains = [c.strip() for c in chain_str.split(',')]
                current_chains.extend(current_matrix_chains)

            elif content.startswith('AND CHAINS:'):
                chain_str = content.split('AND CHAINS:')[1].strip().rstrip(',')
                extra = [c.strip() for c in chain_str.split(',')]
                current_matrix_chains.extend(extra)
                current_chains.extend(extra)

            elif content.startswith('BIOMT'):
                parts = content.split()
                row_idx = int(parts[0][5]) - 1   # BIOMT1/2/3
                row = [float(x) for x in parts[2:6]]  # 3 rotation + 1 translation
                current_matrix_rows.append((row_idx, row))

                # every 3 rows = one 4x4 xform
                if len(current_matrix_rows) == 3:
                    mat = np.eye(4, dtype=np.float32)
                    for ridx, r in current_matrix_rows:
                        mat[ridx, :3] = r[:3]   # rotation
                        mat[ridx,  3] = r[3]    # translation
                    current_xforms.append(mat)
                    current_matrix_rows = []

    flush_assembly()
    return asmb_ids, asmb_details, asmb_method, asmb_chains, asmb_xforms


# ── Metadata .pt ───────────────────────────────────────────────────────────────

def extract_metadata(structure, pdb_path, valid_chain_ids):
    """Build the PDBID.pt metadata dict."""
    header = structure.header

    asmb_ids, asmb_details, asmb_method, asmb_chains, asmb_xforms = \
        parse_remark350(pdb_path)

    # Placeholder TM-score matrix (requires running TM-align externally)
    n = len(valid_chain_ids)
    tm = torch.zeros(n, n, 3, dtype=torch.float32)
    for i in range(n):
        tm[i, i] = torch.tensor([1.0, 1.0, 0.0])  # self-similarity

    return {
        'method':       header.get('structure_method', 'unknown').upper(),
        'date':         header.get('deposition_date', 'unknown'),
        'resolution':   float(header.get('resolution', 0.0) or 0.0),
        'chains':       valid_chain_ids,
        'tm':           tm,
        'asmb_ids':     asmb_ids,
        'asmb_details': asmb_details,
        'asmb_method':  asmb_method,
        'asmb_chains':  asmb_chains,
        'asmb_xforms':  asmb_xforms,  # list of [n,4,4] tensors, one per assembly
    }


# ── Hashing & CSV ──────────────────────────────────────────────────────────────

def seq_hash(seq: str) -> str:
    """6-character hash of the sequence."""
    return hashlib.md5(seq.encode()).hexdigest()[:6].upper()


def write_list_csv(records, output_dir):
    """Write list.csv with one row per chain."""
    csv_path = Path(output_dir) / 'list.csv'
    fieldnames = ['CHAINID', 'DEPOSITION', 'RESOLUTION', 'HASH', 'CLUSTER', 'SEQUENCE']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print(f"  Wrote {csv_path}")


# ── Main conversion ────────────────────────────────────────────────────────────

def convert_pdb(pdb_path, output_dir, cluster_map=None):
    """
    Convert a single PDB file into:
        PDBID_CHAINID.pt  (one per chain)
        PDBID.pt          (metadata)
    and return a list of CSV row dicts.

    cluster_map: optional dict {seq_hash -> cluster_id}
    """
    pdb_path   = Path(pdb_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_id = pdb_path.stem.upper()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))
    model = structure[0]

    header     = structure.header
    date       = header.get('deposition_date', 'unknown')
    resolution = float(header.get('resolution', 0.0) or 0.0)

    valid_chain_ids = []
    csv_rows = []

    for chain in model:
        chain_data = extract_chain_pt(chain)
        if chain_data is None:
            print(f"  Skipping empty chain {chain.id}")
            continue

        chain_label = f"{pdb_id}_{chain.id}"
        valid_chain_ids.append(chain.id)

        out_path = output_dir / f"{chain_label}.pt"
        torch.save(chain_data, out_path)
        print(f"  Saved {out_path.name}  (L={len(chain_data['seq'])})")

        h = seq_hash(chain_data['seq'])
        csv_rows.append({
            'CHAINID':    chain_label,
            'DEPOSITION': date,
            'RESOLUTION': resolution,
            'HASH':       h,
            'CLUSTER':    cluster_map.get(h, '') if cluster_map else '',
            'SEQUENCE':   chain_data['seq'],
        })

    # metadata
    meta = extract_metadata(structure, str(pdb_path), valid_chain_ids)
    meta_path = output_dir / f"{pdb_id}.pt"
    torch.save(meta, meta_path)
    print(f"  Saved {meta_path.name}")

    return csv_rows


def convert_all_pdbs(pdb_dir, output_dir, cluster_map=None):
    """Convert all PDB files in a directory and write list.csv."""
    pdb_dir  = Path(pdb_dir)
    out_dir  = Path(output_dir)
    pdb_files = sorted(pdb_dir.glob("*.pdb")) + sorted(pdb_dir.glob("*.ent"))

    print(f"Found {len(pdb_files)} PDB files\n")
    all_rows = []
    for pdb_file in pdb_files:
        print(f"Processing {pdb_file.name} ...")
        try:
            rows = convert_pdb(pdb_file, out_dir, cluster_map)
            all_rows.extend(rows)
        except Exception as e:
            print(f"  ERROR: {e}")

    write_list_csv(all_rows, out_dir)
    print(f"\nDone. {len(all_rows)} chains written to {out_dir}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        convert_pdb(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python convert.py <input.pdb> <output_dir>")
        print("       python convert.py  (edit __main__ for batch)")