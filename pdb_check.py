import os
import shutil
import pandas as pd

def check_pdb_files(pdb_dir):
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    empty_list = []
    if not pdb_files:
        print("No PDB files found in the directory.")
        return

    print(f"Found {len(pdb_files)} PDB files in the directory.")

    empty_dir = os.path.join(pdb_dir, "empty_files")
    empty_count = 0

    for fname in pdb_files:
        fpath = os.path.join(pdb_dir, fname)
        try:
            file_size = os.stat(fpath).st_size
            if file_size == 0:
                print(f"Warning: Empty PDB file found: {fname}")
                os.makedirs(empty_dir, exist_ok=True)
                shutil.move(fpath, os.path.join(empty_dir, fname))
                empty_count += 1
                empty_list.append(fname)
        except OSError as e:
            print(f"Error checking {fname}: {e}")

    if empty_count:
        print(f"Moved {empty_count} empty file(s) to {empty_dir}")
        empty_records = []
        for empty_file in empty_list:
            pdb_name = empty_file.split("-")[1]  # extract UniProt ID
            empty_records.append({'filename': empty_file, 'uniprot_id': pdb_name})
        
        pd.DataFrame(empty_records).to_csv("empty_pdb_files.csv", index=False)  # ✅ outside loop
    else:
        print("All PDB files are non-empty.")


def find_length(empty_list, check_csv):
    """
    Looks up the sequence length for each empty PDB file.

    Args:
        empty_list : list of dicts with 'filename' and 'uniprot_id'
                     OR a path to empty_pdb_files.csv (str)
        check_csv  : DataFrame with 'acc' and 'length' columns
                     OR a path to the CSV/TSV file (str)

    Returns:
        DataFrame with 'filename', 'uniprot_id', 'length'
    """
    # accept file paths as well as already-loaded data
    if isinstance(empty_list, str):
        empty_list = pd.read_csv(empty_list).to_dict('records')
    if isinstance(check_csv, str):
        check_csv = pd.read_csv(check_csv)

    # sanity checks
    required_keys = {'filename', 'uniprot_id'}
    if empty_list and not required_keys.issubset(empty_list[0].keys()):
        raise ValueError(f"empty_list records must have keys: {required_keys}. "
                         f"Got: {set(empty_list[0].keys())}")
    if 'acc' not in check_csv.columns or 'length' not in check_csv.columns:
        raise ValueError(f"check_csv must have 'acc' and 'length' columns. "
                         f"Got: {check_csv.columns.tolist()}")

    # ensure both ID columns are strings for safe comparison
    check_csv = check_csv.copy()
    check_csv['acc'] = check_csv['acc'].astype(str)

    length_records = []
    for row in empty_list:
        uid = str(row['uniprot_id'])
        match = check_csv[check_csv['acc'] == uid]
        length = match['length'].values[0] if not match.empty else 'Not found'
        length_records.append({
            'filename':   row['filename'],
            'uniprot_id': uid,
            'length':     length,
        })

    df = pd.DataFrame(length_records)
    df.to_csv("empty_pdb_files_with_length.csv", index=False)
    print(f"Wrote {len(df)} records to empty_pdb_files_with_length.csv")
    return df


if __name__ == "__main__":
    # Option A — chain directly
    empty_records = check_pdb_files("/path/to/pdbs/")
    check_csv     = pd.read_csv("mobidb_filtered_data.csv")
    find_length(empty_records, check_csv)

    # Option B — load from saved CSV
    # find_length("empty_pdb_files.csv", "mobidb_filtered_data.csv")
    