import pandas as pd
import os
import numpy as np
from Bio.PDB import PDBList, MMCIFParser, MMCIF2Dict
import time

# Parses ATOM records from an mmCIF file and extracts atomic coordinates and metadata
# Input: file_path (string) – path to the mmCIF file
# Output: tuple of DataFrames (full atom data, backbone-only atoms ["N", "O", "C"])

def atoms (file_path):

    columns = ["group_PDB","model","chain", "id", "type_symbol", "label_atom_id", "label_comp_id","res num", "Cartn_x", "Cartn_y", "Cartn_z"]
    data = []
    
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("ATOM "):  # Filter only ATOM lines
                chunk = line.split()      # Split line by whitespace
                atom_info = [
                    chunk[0],          # ATOM
                    int(chunk[20]),    # Model
                    chunk[6],          # chain
                    int(chunk[1]),     # Atom serial number
                    chunk[2],          # Atom name 
                    chunk[3],          
                    chunk[5],          # Residue name (e.g., VAL, LEU, SER, etc.)
                    int(chunk[8]),     # Residual num
                    float(chunk[10]),  # X coordinate
                    float(chunk[11]),  # Y coordinate
                    float(chunk[12]),  # Z coordinate
                ]  
                data.append(atom_info)
    
    df = pd.DataFrame(data, columns=columns)
    df_backbone = df[df["label_atom_id"].isin(["N", "O", "C"])]

    return(df, df_backbone)

# Parses HELX_P records from an mmCIF file and extracts helix metadata
# Input: file_path (string) – path to the mmCIF file
# Output: DataFrame containing helix information (start/end positions, chain, IDs)

def helix (file_path):

    data = []
    columns = ['helix','Helix id','chain start','start','start res','chain end','end','end res']
    atom_info = []
    
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("HELX_P "):  # Filter only ATOM lines
                chunk = line.split()  # Split line by whitespace
                atom_info = [
                    chunk[1],       #Helix
                    chunk[2],       #Helix id
                    chunk[4],       #Chain 
                    chunk[5],       #Start
                    chunk[3],       #Start res
                    chunk[8],       #Chain
                    chunk[9],       #end
                    chunk[7],       #end res
                ]
                data.append(atom_info)
    
    helix_df = pd.DataFrame(data, columns=columns)

    return(helix_df)

# Parses sheet section from an mmCIF file and extracts strand metadata
# Input: file_path (string) – path to the mmCIF file
# Output: DataFrame containing sheet/strand information (start/end positions, chain, IDs)

def sheet (file_path):

    sheet_lines = []
    found_header = False
    
    with open(file_path) as file:
        for line in file:
            if not found_header:
                if "_struct_sheet_range.end_auth_seq_id" in line:
                    found_header = True
                continue
            if line.startswith("#"):
                break
            sheet_lines.append(line.strip())
    
    data = []
    
    columns = ['Sheet id','segment id','chain start','start','start res','chain','end','end res']
    atom_info = []
    
    for line in sheet_lines:  # Filter only ATOM lines
        chunk = line.split()  # Split line by whitespace
        atom_info = [
            chunk[0],       #sheet-level identifier - groups together segments (or strands)
            chunk[1],       #id segment within sheet
            chunk[3],       #Chain 
            chunk[4],       #Start
            chunk[2],       #start res
            chunk[7],       #Chain 
            chunk[8],       #Start
            chunk[6],       #start res
        ]
        data.append(atom_info)
    
    sheet_df = pd.DataFrame(data, columns=columns)

    return(sheet_df)

# Creates= secondary structure for each residue in a chain
# Input: helix_df (DataFrame), sheet_df (DataFrame), chain (string), chain_len (int)
# Output: list of secondary structure labels per residue (e.g., "loop", "sheet", "helix1")

def SS(helix_df, sheet_df, chain, chain_len):

    helix_df = helix_df[helix_df['chain start'] == chain].reset_index(drop=True)
    sheet_df = sheet_df[sheet_df['chain start'] == chain].reset_index(drop=True)
    helix_names = list(helix_df['helix'])

    #print(helix_df)
    #print(sheet_df)
    
    helix_start_end = [] 
    for i in range(len(helix_df)):
        start = helix_df.loc[i, 'start']
        end = helix_df.loc[i, 'end']
        helix_start_end.append((start, end))

    # Extract sheet start and end positions as tuples
    sheet_start_end = []
    for i in range(len(sheet_df)):
        start = sheet_df.loc[i, 'start']
        end = sheet_df.loc[i, 'end']
        sheet_start_end.append((start, end))

    # Create a dictionary mapping each residue in a helix to its helix name
    expanded_helix_dict = {}
    for idx, (start, end) in enumerate(helix_start_end):
        helix_name = helix_names[idx]  # Assume the order of names matches the ranges
        for pos in range(int(start), int(end) + 1):
            expanded_helix_dict[pos] = helix_name

    # Create a list of all residues that are part of sheets
    expanded_sheet_list = []
    for start, end in sheet_start_end:
        expanded_sheet_list.extend(range(int(start), int(end) + 1))

    # Build the final list of secondary structure annotations for each residue
    helix_sheet_loop = []
    for i in range(1, chain_len + 1):  # Using 1-based indexing for residues
        if (i in expanded_helix_dict) and (i in expanded_sheet_list):
            # Residue in both helix and sheet: include the helix name along with "and sheet"
            helix_sheet_loop.append(expanded_helix_dict[i] + " and sheet")
        elif i in expanded_helix_dict:
            # Residue in helix: include the helix name
            helix_sheet_loop.append(expanded_helix_dict[i])
        elif i in expanded_sheet_list:
            # Residue in sheet only
            helix_sheet_loop.append("sheet")
        else:
            # Residue is a loop (neither helix nor sheet)
            helix_sheet_loop.append("loop")
    
    return helix_sheet_loop

# Computes Euclidean distance between two atoms (3D coordinates)
# Input: atom1 (Pandas Series), atom2 (Pandas Series or NaN)
# Output: float distance or NaN if atom2 is missing

def dist(atom1, atom2):
    
    if isinstance(atom2, float) and np.isnan(atom2):
        return np.nan
    
    x1, y1, z1 = atom1['Cartn_x'], atom1['Cartn_y'], atom1['Cartn_z']
    x2, y2, z2 = atom2['Cartn_x'], atom2['Cartn_y'], atom2['Cartn_z']
    
    distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance

# Retrieves an atom (row) at a given position from a DataFrame
# Input: pos (int), df (DataFrame)
# Output: Pandas Series representing the atom or NaN if out of range

def atom(pos, df):
    try:
        return df.iloc[pos]
    except IndexError: 
        return np.nan

# Computes a set of pairwise distances and differences for a residue and its neighbors
# Input: i (int) – residue index, df (DataFrame) – backbone atom data
# Output: list [index, residue, SS, distances and differences for step 2 to 5]

def distance_metrics(i, df):

    BD2 = dist(atom(3*i+2, df), atom(3*i+6, df))
    AD2 = dist(atom(3*i+1, df), atom(3*i+6, df))
    AD2_BD2 = AD2 - BD2
    
    BD3 = dist(atom(3*i+2, df), atom(3*i+9, df))
    AD3 = dist(atom(3*i+1, df), atom(3*i+9, df))
    AD3_BD3 = AD3 - BD3
    
    BD4 = dist(atom(3*i+2, df), atom(3*i+12, df))
    AD4 = dist(atom(3*i+1, df), atom(3*i+12, df))
    AD4_BD4 = AD4 - BD4
    
    BD5 = dist(atom(3*i+2, df), atom(3*i+15, df))
    AD5 = dist(atom(3*i+1, df), atom(3*i+15, df))
    AD5_BD5 = AD5 - BD5

    info = atom(3*i, df)
    res = info['label_comp_id']
    SS = info['PDB secondary structure']
    
    return([i+1, res, SS , BD2, AD2_BD2, BD3, AD3_BD3, BD4, AD4_BD4, BD5, AD5_BD5])

# Highlights cells in a DataFrame based on value and column type
# Input: val (any), col_name (string)
# Output: string – cell style (e.g., "background-color: yellow") or empty string

def highlight_cells(val, col_name):
    if val == "NaN":
        return ""
    try:
        v = float(val)
    except ValueError:
        return ""
    if "distance" in col_name:
        if v < 4.0:
            return "background-color: yellow"
    elif "difference" in col_name:
        if v >= 0.9:
            return "background-color: yellow"
    return ""

# Applies conditional formatting to a DataFrame column
# Input: s (Pandas Series)
# Output: list of cell styles for the column

def highlight_column(s):
    # s is a Series; s.name gives the column name.
    return [highlight_cells(v, s.name) for v in s]

# Main function to compute backbone H-bond distance metrics and save results to Excel
# Input: pdb_id (string), pdb_file_loc (string), output_loc (string)
# Output: None (saves styled Excel file with distance metrics)

def distances (pdb_id, pdb_file_loc, output_loc):

    pdbl = PDBList()
    
    #pdb_id = "1hho"
    pdb_id = pdb_id.upper()
    file_path = pdbl.retrieve_pdb_file(pdb_id, file_format="mmCif", pdir=pdb_file_loc)
    print(f"File saved at: {file_path}")
    
    cleaned = pd.read_csv("PDB707K_cleaned_chains_sequences_19Feb2025.csv")
    
    df, df_backbone = atoms(file_path)
    helix_df = helix(file_path)
    sheet_df = sheet(file_path)

    models = df_backbone['model'].unique()
    chains = df_backbone['chain'].unique()
    models_chains_cleaned = cleaned[cleaned['pdb_id'] == pdb_id]

    print(f"models: {models} chains: {chains}")
    
    for model in models:
    
        model_df = df_backbone[df_backbone['model'] == model]
        
        for chain in chains:
    
            models_chains_cleaned = cleaned[cleaned['pdb_id'] == pdb_id]
            cleaned_models = list(models_chains_cleaned['model_id'].unique()) #unique models is from cleaned data
            print(f"model: {model}, chain: {chain}")
    
            if model in cleaned_models:
                # Get the list of chains for this model_id from the cleaned data
                cleaned_chains = list(models_chains_cleaned.loc[models_chains_cleaned['model_id'] == model, 'chain_id'])
                if chain in cleaned_chains:
                    print(f"✅ Model {model} and chain {chain} in cleaned data -- Compute bonds")
                    compute_distances = True
                else:
                    print(f"❌Model {model} but not chain {chain} in cleaned data")
                    compute_distances = False
            else:
                print(f"❌ Model {model} is not in cleaned data")
                compute_distances = False
    
    
            if compute_distances == True:
    
                model_chain_df = model_df[model_df['chain'] == chain]
                chain_len = int(len(model_chain_df)/3)
                SS_list = SS(helix_df, sheet_df, chain, chain_len+1)
                rep_SS_list = [ss for ss in SS_list for _ in range(3)]
                rep_SS_list = rep_SS_list[:len(model_chain_df)]
                model_chain_df.insert(2, "PDB secondary structure", rep_SS_list)
                chain_start = '1'
                entity_id = models_chains_cleaned.loc[(models_chains_cleaned['model_id'] == model) & (models_chains_cleaned['chain_id'] == chain), 'entity_id'].values[0]
        
                data = []
                
                for i in range (0,chain_len):
                    data.append(distance_metrics(i, model_chain_df))
        
                columns = [
                    "index",
                    "residue",
                    "PDB secondary structure",
                    "distance |O N+2|",
                    "difference |C N+2| - |O N+2|",
                    "distance |O N+3|",
                    "difference |C N+3| - |O N+3|",
                    "distance |O N+4|",
                    "difference |C N+4| - |O N+4|",
                    "distance |O N+5|",
                    "difference |C N+5| - |O N+5|"
                ]
                
                model_chain_dist = pd.DataFrame(data, columns=columns)
                
                excel_filename = f"{output_loc}/{pdb_id}/{pdb_id}-{entity_id}-{model}-{chain}-{chain_start}-{chain_len}_short_hydro_bonds.xlsx"
                protein_dir = f"{output_loc}/{pdb_id}"
        
                os.makedirs(protein_dir, exist_ok=True)
                
                print(excel_filename)
        
                # Determine which columns should be styled (those containing "distance" or "difference")
                highlight_cols = [col for col in model_chain_dist.columns if "distance" in col or "difference" in col]
                
                # Apply the highlighting function to each of these columns
                styled_df = model_chain_dist.style.apply(highlight_column, subset=highlight_cols)
                
                # Write the styled DataFrame to an Excel file using an ExcelWriter.
                with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
                    styled_df.to_excel(writer, index=False, sheet_name='Sheet1')
                
                print(f"Excel file saved at: {excel_filename}")


df = pd.read_csv("PDB707K_cleaned_chains_sequences_19Feb2025.csv")

df_9 = df[df["pdb_id"].astype(str).str.startswith("9")]
proteins = df_9['pdb_id'].unique() 

print(f"Number of unique proteins {len(proteins)} number of models/chains {len(df_9)}")

print(proteins)
print("---------------------------")

input_path = "mmcif_files9"
output_path = "new_output_files_9"

failed_proteins = []

start_time = time.time() 

for i, pdb_id in enumerate(proteins, start=1):
    print(f"Protein {i}/{len(proteins)}: {pdb_id}")
    try:
        distances(pdb_id, input_path, output_path)
    except:
        print(f"❌❌❌could not compute distances for {pdb_id}")
        failed_proteins.append(pdb_id)
        
    print()

end_time = time.time()  # End the timer
elapsed_time = end_time - start_time

print("------")
print(f"Time taken to run: {elapsed_time:.2f} seconds")
print(f"failed : {failed_proteins}")
print("------")
