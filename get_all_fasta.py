import requests
import concurrent.futures
import time

'''This file contains functions to download FASTA files for a list of PDBids from RCSB. 
Unortunately there is a limit to how often the client can make API requests, so there is 
a sleep timer included and the maximum number of workers is set to 5. If you find a batch download method
please consider modifying this part.'''


def download_fasta_from_pdb(pdb_ids, logger):
    '''Downloads FASTA for a pdbid and logs some warnings'''
    
    # Consider modifying adding more warnings for the logger if needed
    fasta_sequences = {}
    def download(pdb_id):
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            fasta_sequences[pdb_id] = response.text
        elif response.status_code == 404:
            logger.warning(f"FASTA not found for PDB ID: {pdb_id}, status code: " 
                           f"{response.status_code} \n Response: {response.text}")
            # print(f"FASTA not found for PDB ID: {pdb_id}")
            fasta_sequences[pdb_id] = None
        elif response.status_code == 429:
            logger.warning(f"Access to FASTA for PDB ID: {pdb_id}" 
                           " is forbidden due to too many requests")
            # print(f"Access to FASTA for PDB ID: {pdb_id}" 
            #                " is forbidden due to too many requests")
            # This can be adjusted if we need to increase the pause
            time.sleep(5)
            download(pdb_id)
        else:
            logger.warning(f"Failed to download FASTA for PDB ID: {pdb_id}, status code: " 
                           f"{response.status_code} \n Response: {response.text}")
            # print(f"Failed to download FASTA for PDB ID: {pdb_id}")
            fasta_sequences[pdb_id] = None
            
    # max_workers can be further reduced in case the server is overwhelmed and responds with 429
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        executor.map(download, pdb_ids)

    return fasta_sequences

def get_fasta(input_file, output_path, logger):
    '''reads in a list of pdb ids and downloads their Fasta in a single file'''
    with open(input_file, "r") as f:
        pdb_ids = f.read().splitlines()
    fasta_sequences = download_fasta_from_pdb(pdb_ids, logger)
    with open(output_path, "w") as f:
        for pdb_id, fasta_content in fasta_sequences.items():
            if fasta_content:
                f.write(str(fasta_content))

