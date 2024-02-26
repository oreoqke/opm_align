import requests
import concurrent.futures
import time

def download_fasta_from_pdb(pdb_ids, logger):
    fasta_sequences = {}
    def download(pdb_id):
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            fasta_sequences[pdb_id] = response.text
        elif response.status_code == 404:
            logger.warning(f"FASTA not found for PDB ID: {pdb_id}")
            print(f"FASTA not found for PDB ID: {pdb_id}")
            fasta_sequences[pdb_id] = None
        elif response.status_code == 429:
            logger.warning(f"Access to FASTA for PDB ID: {pdb_id}" 
                           " is forbidden due to too many requests")
            print(f"Access to FASTA for PDB ID: {pdb_id}" 
                           " is forbidden due to too many requests")
            time.sleep(5)
            download(pdb_id)
        else:
            logger.warning(f"Failed to download FASTA for PDB ID: {pdb_id}")
            print(f"Failed to download FASTA for PDB ID: {pdb_id}")
            fasta_sequences[pdb_id] = None

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        executor.map(download, pdb_ids)

    return fasta_sequences

def get_fasta(input_file, output_path, logger):
    with open(input_file, "r") as f:
        pdb_ids = f.read().splitlines()
    fasta_sequences = download_fasta_from_pdb(pdb_ids, logger)
    with open(output_path, "w") as f:
        for pdb_id, fasta_content in fasta_sequences.items():
            if fasta_content:
                f.write(str(fasta_content))

