""" Module Providing the SuperpositionOPM class. """
import time
import os
import sys
import pandas as pd
from tqdm import tqdm
import urllib.request
import re
import concurrent.futures
import psutil
import queue
import threading
import argparse


import warnings
warnings.filterwarnings("ignore", message="The behavior of DataFrame concatenation with empty or all-NA entries is deprecated.*")


def input_parser():
    parser = argparse.ArgumentParser(
        description='Superposition OPM structures')
    # Lists of reference and target structures as pdbids
    parser.add_argument('--reference_list', type=str, required=False, const=None,
                        help='List of reference structures', default='reference_pdb.csv')
    parser.add_argument('--target_list', type=str, required=False, const=None,
                        help='List of target structures', default='list_pdb')
    # Custom path to USalign executable
    parser.add_argument('--usalign', type=str, required=False, const=None,
                        help='Path to USalign executable', default='script/USalign')
    # Parameters for classification
    parser.add_argument('--P_over_protein', type=float, required=False, const=None,
                        help='Minimal percentage of overlapped residues in the superposition (same protein)', default=95)
    parser.add_argument('--P_ide_protein', type=float, required=False, const=None,
                        help='Minimal percentage of identity (same protein)', default=98)
    parser.add_argument('--P_over_family', type=float, required=False, const=None,
                        help='Minimal percentage of overlapped residues in the superposition (same family)', default=25)
    parser.add_argument('--P_ide_family', type=float, required=False, const=None,
                        help='Minimal percentage of identity (same family)', default=50)
    parser.add_argument('--rmsd_same',  type=float, required=False, const=None,
                        help='Maximal RMSD of the superposition (A)', default=1.0)
    parser.add_argument('--mode', '-m', type=str, required=False, const=None, 
                        help='Mode of superposition\n 1: -ter 1 -byresi 1\n 2: -ter 1 -byresi 0\n 3:-ter 1 -mm 1\n Else: set by user', default='3')
    parser.add_argument('--rmsd_max', type=float, required=False, const=None,
                        help='Maximal RMSD of the superposition (A)', default=5.0)
    
    return parser.parse_args()


# Perform superposition using USalign
# Iterate over the elements of target_list
def run_align(pdb_folder, options, target, ref):
    # Run USalign and get the output
    pipe = os.popen(f'./USalign {pdb_folder}/{target}.cif {pdb_folder}/{ref}.cif -outfmt 1 {options} 2>/dev/null')
    output_splits = pipe.read()
    pipe.close()
    result_df_1 = pd.DataFrame(columns=['PDB_ID_1', 'chain_1', 'PDB_ID_2', 
                                  'chain_2', 'length_1', 'length_2', 
                                  'TM-score_1', 'TM-score_2', 'd0_1', 'd0_2', 
                                  'RMSD', 'Lali', 'seqid_1', 'seqid_2','seqid_ali', 
                                  'N_over', 'N_ide', 'P_over', 'P_seqide', 'P_ide_1', 
                                  'P_ide_2', 'Nsub_over'])
    # Define a pattern to search for the desired information
    pattern = fr'>{pdb_folder}/{target}.cif:(\S+)\s+L=(\S+)\s+d0=(\d+\.\d+)\s+seqID=(\d+\.\d+)\s+TM-score=(\d+\.\d+)\n(\S+)\n>{pdb_folder}/{ref}.cif:(\S+)\s+L=(\S+)\s+d0=(\d+\.\d+)\s+seqID=(\d+\.\d+)\s+TM-score=(\d+\.\d+)\n(\S+)\n# Lali=(\d+)\s+RMSD=(\d+\.\d+)\s+seqID_ali=(\d+\.\d+)'
    
    # Search for the pattern in the output of USalign
    for result in re.finditer(pattern, output_splits):
        if result is None:
            print(f'US align failure for {target} and {ref}:')
            print(output_splits)
            if not sys.stdout.isatty():
                sys.stdout.flush()
            continue

        # Extract the information from the match
        chain_1 = result.group(1)
        l_1 = int(result.group(2))
        d0_1 = float(result.group(3))
        seqid_1 = float(result.group(4))
        score_1 = float(result.group(5))
        seq_1 = result.group(6)
        chain_2 = result.group(7)
        l_2 = int(result.group(8))
        d0_2 = float(result.group(9))
        seqid_2 = float(result.group(10))
        score_2 = float(result.group(11))
        seq_2 = result.group(12)
        Lali = int(result.group(13))
        rmsd = float(result.group(14))
        seqid_ali = float(result.group(15))
            
        # Leave if the length of the sequence is less than 20
        if min(l_1, l_2) < 20:
            continue
            
        N_over = 0
        N_identical = 0
        Nsub_over = 0
        star = True
        
        # Calculate the number of residues that are overlapped
        for s1, s2 in zip(seq_1, seq_2):
            if s1 != '*' and s1 != '-' and s2 != '-':
                N_over += 1
                if s1 == s2:
                    N_identical += 1
            if s1 == '*':
                star = True
            if star == True and s1 != '-' and s2 != '-' and s1 != '*': 
                Nsub_over += 1
                star = False
                    
        P_over = round(100 * N_over / min(l_1, l_2))
        P_seqide = round(100 * N_identical / N_over)
        rmsd = round(rmsd, 2)
        P_identical_1 = round(100 * N_identical / l_1)
        P_identical_2 = round(100 * N_identical / l_2)                        
        
        # Append the extracted information to the data frame as a new row
        to_add = pd.DataFrame([[target, chain_1, ref, chain_2, l_1, l_2, score_1, score_2, d0_1, d0_2, rmsd, Lali, seqid_1, seqid_2, seqid_ali, N_over, N_identical, P_over, P_seqide, P_identical_1, P_identical_2, Nsub_over]], columns=result_df_1.columns)
        result_df_1 = pd.concat([result_df_1, to_add], ignore_index=True)
    
    return result_df_1

# Shared queue for results
results_queue = queue.Queue()

# Flag to indicate when all tasks are done
all_tasks_done = threading.Event()

# Function to write results to file
def write_results_to_file(output_dir):
    while not all_tasks_done.is_set() or not results_queue.empty():
        try:
            #dframe = pd.DataFrame(['PDB_tar', 'chain_tar', 'PDB_ref', 'chain_ref', 'type', 'N_over', 'P_over', 'P_seqide', 'RMSD', 'P_ide_1', 'P_ide_2'])
            resultat = results_queue.get(timeout=1)  # Timeout to avoid blocking indefinitely
            #dframe = pd.concat([dframe, resultat])
            resultat.to_csv(f'{output_dir}/results.csv', mode='a', header=False, index=False)

            results_queue.task_done()
        except queue.Empty:
            continue


def align_structures(output_dir, logger):
    """Align the structures using the superposition algorithm."""
    args = input_parser()
    exe = "USalign"
    P_over_protein = args.P_over_protein
    P_over_family = args.P_over_family
    P_ide_protein = args.P_ide_protein
    P_ide_family = args.P_ide_family
    rmsd_same = args.rmsd_same
    rmsd_max = args.rmsd_max
    mode = args.mode
    
    
    start_time = time.time()
    # Check if Executable USalign is present
    if os.path.isfile(exe) is False:
        print(f"Error: {exe} is not found.")
        logger.warning(f"Error: {exe} is not found.")
        logger.error(f"Error: {exe} is not found.")
        sys.exit(1)
    
    # Set USalign options
    if mode == '1':
        options = '-ter 1 -byresi 1'
    elif mode == '2':
        options = '-ter 1 -byresi 0'
    elif mode == '3':
        options = '-ter 1 -mm 1'
    else:
        options = mode

    # Create a list of all PDB IDs that we will use to align
    ref_list = pd.read_csv(os.path.join(output_dir, "new_pdb.txt"),
                           header=None, sep=" ")[0].tolist()
    target_list = pd.read_csv(os.path.join(output_dir, "to_align.txt"),
                              header=None, sep=" ")[1].tolist()
    all_pdb = list(set().union(ref_list, target_list))
    #print(all_pdb)
    pdb_folder = os.path.join(output_dir, "pdb")
    if os.path.isdir(pdb_folder) is False:
        os.mkdir(pdb_folder)
    #print(target_list)
    #print(ref_list)
    
    # Fetch mmCIF files from RCSB
    for i in tqdm(range(len(all_pdb)), desc='\033[94m'+'Downloading mmCIF files from RCSB' + '\033[0m', disable=None, leave=False):
        pdb = all_pdb[i]
        if os.path.isfile(f'{pdb_folder}/{pdb}.cif') is False:
            try:
                urllib.request.urlretrieve('https://files.rcsb.org/download/'+pdb+'.cif', f'{pdb_folder}/{pdb}.cif')
            except:
                print('\033[91m'+'Error: '+pdb+' does not exist on RCSB'+'\033[0m')
                target_list = target_list[target_list != pdb]
                ref_list = ref_list[ref_list != pdb]
            continue
    print('\033[92m'+'All mmCIF files are downloaded.'+'\033[0m')
    
    # Create an empty data frame to store the results
    result_df = pd.DataFrame(columns=['PDB_ID_1', 'chain_1', 'PDB_ID_2', 
                                  'chain_2', 'length_1', 'length_2', 
                                  'TM-score_1', 'TM-score_2', 'd0_1', 'd0_2', 
                                  'RMSD', 'Lali', 'seqid_1', 'seqid_2','seqid_ali', 
                                  'N_over', 'N_ide', 'P_over', 'P_seqide', 'P_ide_1', 
                                  'P_ide_2', 'Nsub_over', 'seq_ide_blast', 'length_blast'])
    
    # if os.path.isfile(os.path.join(output_dir, "result.csv")):
    #     result_df = pd.read_csv(os.path.join(output_dir, "result.csv"), header=None, sep="\t")

    to_skip = []
    target_list = open(os.path.join(output_dir, "to_align.txt"), "r").read().splitlines()
    
    
    # This is a single thread version that may take longer than your PhD to compute
    # for target in tqdm(target_list, desc='\033[94m'+'Performing superposition using USalign'+ '\033[0m',disable=None, postfix='targets', leave = False):        # Iterate over the elements of ref_list
    #     ref, target = target.split()
    #     target_time = time.time()
    #     # for ref in tqdm(ref_list, leave = False, disable=None, postfix='references'):
    #     result_df = run_align(pdb_folder, output_dir, options, result_df, to_skip, target, ref)  
    #     if not sys.stdout.isatty(): sys.stdout.flush()



    # Start the writer thread
    writer_thread = threading.Thread(target=write_results_to_file, args=(output_dir,))
    writer_thread.start()
    
    # Create a single tqdm progress bar to track overall progress
    overall_progress = tqdm(total=len(target_list), desc='Processing references', leave=False)

    # Create a single tqdm progress bar to track overall progress
    if os.path.isfile(f'{output_dir}/results.csv'):
        cmd = f"rm -f {output_dir}/results.csv"
        os.system(cmd)
    result_df.to_csv(f'{output_dir}/results.csv', mode='a', index=False)

    # Count the total number of references
    def process_reference(target, ref, seq_ide_blast, length_blast):
        resultat = run_align(pdb_folder, options, target, ref)
        # add the seq_ide_blast to the result
        resultat['seq_ide_blast'] = seq_ide_blast
        resultat['length_blast'] = length_blast
        results_queue.put(resultat)

    # This creates a pool of workers, consider using a smaller number of workers if you need to run other tasks
    with concurrent.futures.ThreadPoolExecutor(max_workers=psutil.cpu_count(logical=False)) as executor:
        futures = []
        for target in target_list:
            target, ref, seq_ide_blast, length_blast = target.split()
            futures.append(executor.submit(process_reference, target, ref, seq_ide_blast, length_blast))
        
        for future in concurrent.futures.as_completed(futures):
            overall_progress.update(1)
    concurrent.futures.wait(futures)
    all_tasks_done.set()
    writer_thread.join() 
            
    print('\033[92m'+f'Finished superpositioning using USalign in {time.time() - start_time:.2f} seconds.'+'\033[0m')

    sys.stdout.flush()

    classifier_time = time.time()
    def classify_match(group):
        def classify_match_helper(row):
            if row['RMSD'] > rmsd_max:
                return 0
            elif row['P_over'] >= P_over_protein and row['P_ide_1'] >= P_ide_protein and row['RMSD'] <= rmsd_same:
                return 1
            elif row['P_ide_1'] >= P_ide_protein:
                return 2
            elif row['P_over'] >= P_over_family and row['P_ide_1'] >= P_ide_family:
                return 3
            elif  row['P_over'] >= P_over_family:
                return 4
            else:
                return 0

        group['type'] = group.apply(classify_match_helper, axis=1)
        return group

    result_df = pd.concat([result_df,pd.read_csv(f'{output_dir}/results.csv')], ignore_index=True)

    selection_df = result_df.groupby(['PDB_ID_1', 'chain_1', 'PDB_ID_2', 'chain_2'], group_keys=False).apply(classify_match)

    # drop the rows that classified to type 0
    try:
        selection_df = selection_df[selection_df['type'] != 0]
    except:
        print('No matches found. Check the parameters and results.csv file.')

    # for each chain_1 in the selection_df, delete the ':' and get length of the string
    if mode == '3':
        for row in selection_df.itertuples():
            # count how many there is string between ':' and ':' in the chain_1
            count_1 = [x for x in row.chain_1.split(':') if x]
            selection_df.at[row.Index, 'chain_1'] = sum(1 for x in count_1)
            count_2 = [x for x in row.chain_2.split(':') if x]
            selection_df.at[row.Index, 'chain_2'] = sum(1 for x in count_2)

        # rename to format
        selection_df.rename(columns = {'PDB_ID_1': 'PDB_tar','PDB_ID_2': 'PDB_ref', 'chain_1': 'Nsub_tar', 'chain_2': 'Nsub_ref'}, inplace=True)
        selection_df = selection_df[['PDB_tar', 'PDB_ref', 'type', 'N_over', 'P_over', 'P_seqide', 'RMSD', 'P_ide_1', 'P_ide_2', 'Nsub_tar', 'Nsub_ref', 'Nsub_over']]
    
    else: 
        selection_df.rename(columns = {'PDB_ID_1': 'PDB_tar','PDB_ID_2': 'PDB_ref', 'chain_1': 'chain_tar', 'chain_2': 'chain_ref'}, inplace=True)
        selection_df = selection_df[['PDB_tar', 'chain_tar', 'PDB_ref', 'chain_ref', 'type', 'N_over', 'P_over', 'P_seqide', 'RMSD', 'P_ide_1', 'P_ide_2']]
        
    selection_df.to_csv(f'{output_dir}/selection.csv', index=False)
    if mode == '3':
        selection_df.groupby(['PDB_tar', 'PDB_ref', 'type'])['RMSD'].idxmin()
        selection_df = selection_df.drop_duplicates(subset=['PDB_tar', 'type'], keep='first')
    else: 
        selection_df.groupby(['PDB_tar', 'chain_tar', 'PDB_ref', 'type'])['RMSD'].idxmin()
        selection_df = selection_df.drop_duplicates(subset=['PDB_tar', 'chain_tar', 'type'], keep='first')
    # Drop values that are not the minimum RMSD
    selection_df.to_csv(f'{output_dir}/top_matches.csv', index=False)
    print('\033[92m'+f'Finished classifying matches in {time.time() - classifier_time:.2f} seconds.'+'\033[0m')
        
#align_structures('results', 'logger')