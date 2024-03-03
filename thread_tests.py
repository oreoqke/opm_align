import concurrent.futures
import queue
import time
import threading
import argparse
import os
import psutil
import sys
from psutil._common import bytes2human
from tqdm import tqdm
import pandas as pd
import glob
import urllib
import re

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
    
    # directories
    parser.add_argument('-i', '--input_dir', type=str, required=False, const=None,
                        help = 'Input directory', default='run_1')
    parser.add_argument('-o', '--output_dir', type=str, required=False, const=None,
                        help='Output directory', default='run_1/results')
    return parser.parse_args()



def run_align(pdb_folder, options, target, ref):
    # Run USalign and get the output
    pipe = os.popen(f'./script/USalign {pdb_folder}/{target}.cif {pdb_folder}/{ref}.cif -outfmt 1 {options} 2>/dev/null')
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
                    
        P_over = round(100 * N_over / min(l_1, l_2), 2)
        P_seqide = round(100 * N_identical / N_over, 2)
        P_identical_1 = round(100 * N_identical / l_1, 2)
        P_identical_2 = round(100 * N_identical / l_2, 2)                        
        
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


# Main function
def main():
    start_time = time.time()
    testing = False
    exact_matching  = False 
    pdb_folder = 'pdb'
    
    # Parse arguments and assign them to variables
    args = input_parser()
    ref_list = args.reference_list
    target_list = args.target_list
    exe = args.usalign
    P_over_protein = args.P_over_protein
    P_over_family = args.P_over_family
    P_ide_protein = args.P_ide_protein
    P_ide_family = args.P_ide_family
    rmsd_same = args.rmsd_same
    rmsd_max = args.rmsd_max
    mode = args.mode

        # Input and output directories
    if args.input_dir == '':
        input_dir = ''
    elif args.input_dir[-1] != '/':
        input_dir = args.input_dir + '/'
    else:
        input_dir = args.input_dir
    
    output_dir = args.output_dir
    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    if not sys.stdout.isatty():
        sys.stdout = open(f'{output_dir}/output.log', "w")
        print('Output redirected to output.log')
        print('PID of current proccess: ', os.getpid())
        print('Job started at: ', time.strftime('%X %x %Z'))
        print('Current working directory: ', os.getcwd())
        print('System platform: ', sys.platform)
        print('Python version: ', sys.version)
        print('OS name: ', os.name)

        print('Number of CPU cores: ', os.cpu_count())
        print('Total RAM: ', bytes2human(psutil.virtual_memory().total))
        print('Total disk space: ', bytes2human(psutil.disk_usage('/').total))

        print('Number of processes: ', len(psutil.pids()))
        print('Number of threads: ', psutil.cpu_count(logical=False))
        sys.stdout.flush()
        
    # Check USalign executable
    if os.path.isfile(exe) is False:
        print('Error: USalign executable not found')
        sys.exit()
    
    # Set USalign options  
    if mode == '1':
        options = '-ter 1 -byresi 1'
    elif mode == '2':
        options = '-ter 1 -byresi 0'
    elif mode == '3':
        options = '-ter 1 -mm 1'
    else:
        options = mode
        
        # Read the reference and target lists (testing mode or full version)
    if testing is False and exact_matching is False:
        ref_list = pd.read_csv(input_dir+ref_list)
        ref_list.columns = ['Family', 'Species', 'Name', 'PDB', 'Superfamily']
        target_list = pd.read_csv(input_dir+target_list, header=None)
        target_list.columns = ['PDB']
    elif testing is True and exact_matching is False:
        ref_list = pd.read_csv('test/ref_pdb', header=None)
        ref_list.columns = ['PDB']
        target_list = pd.read_csv('test/target_pdb', header=None)
        target_list.columns = ['PDB']
        print('\033[93m' + 'Important Notice: ' + '\033[0m' +
            'It\'s a test case, reading from test folder. To run the full version, please set testing to False in the script.')
    elif exact_matching is True:
        data = pd.read_csv(input_dir+'structures_matched.txt', sep = '|', header=None)
        target_list = data.iloc[:,1].to_frame()
        target_list.columns = ['PDB']
        ref_list = data.iloc[:,2].to_frame()
        ref_list.columns = ['PDB']
        count = 0

    # Remove blank PDB ID, duplicates, spaces in names. Make PDB IDs lowercase for target_list and ref_list
    target_list = target_list.dropna()
    target_list['PDB'] = target_list['PDB'].str.strip()
    target_list['PDB'] = target_list['PDB'].str.lower()

    ref_list = ref_list.dropna()
    ref_list['PDB'] = ref_list['PDB'].str.strip()
    ref_list['PDB'] = ref_list['PDB'].str.lower()

        # drop duplicates is only used if we don't do 1-1 matching
    if exact_matching is False:
        ref_list = ref_list.drop_duplicates()
        target_list = target_list.drop_duplicates()

    # Create a list of all PDB IDs
    all_pdb = list(set().union(list(target_list['PDB']), list(ref_list['PDB'])))
    downloaded = glob.glob(f'{pdb_folder}/*.cif')
    all_pdb = list(set(all_pdb)-set(downloaded))
    if os.path.isdir('pdb') is False:
        os.mkdir('pdb')

        # Fetch mmCIF files from RCSB
    for i in tqdm(range(len(all_pdb)), desc='\033[94m'+'Downloading mmCIF files from RCSB' + '\033[0m', disable=None, leave=False):
        pdb = all_pdb[i]
        if os.path.isfile(f'{pdb_folder}/{pdb}.cif') is False:
            try:
                urllib.request.urlretrieve('https://files.rcsb.org/download/'+pdb+'.cif', f'{pdb_folder}/{pdb}.cif')
            except:
                print('\033[91m'+'Error: '+pdb+' does not exist on RCSB'+'\033[0m')
                target_list = target_list[target_list['PDB'] != pdb]
                ref_list = ref_list[ref_list['PDB'] != pdb]
                continue
    print('\033[92m'+'All mmCIF files are downloaded.'+'\033[0m')
    if not sys.stdout.isatty(): sys.stdout.flush()
    
    
    # Create an empty data frame for results
    result_df = pd.DataFrame(columns=['PDB_ID_1', 'chain_1', 'PDB_ID_2', 
                                    'chain_2', 'length_1', 'length_2', 
                                    'TM-score_1', 'TM-score_2', 'd0_1', 'd0_2', 
                                    'RMSD', 'Lali', 'seqid_1', 'seqid_2','seqid_ali', 
                                    'N_over', 'N_ide', 'P_over', 'P_seqide', 'P_ide_1', 
                                    'P_ide_2', 'Nsub_over'])

    # Read 'results.csv' using pandas (if it exists) and read unique pairs of PDB_ID_1 and PDB_ID_2 into a list
    if os.path.isfile(f'{output_dir}/results.csv') is True:
        result_df = pd.read_csv(f'{output_dir}/results.csv')
        to_skip = result_df[['PDB_ID_1', 'PDB_ID_2']].to_records(index=False).tolist()
        to_skip = list(set(to_skip))
    else:
        to_skip = []
        result_df.to_csv(f'{output_dir}/results.csv', index=False, mode='a')
    # Create an empty data frame for results with header
    
    
    
    # Start the writer thread
    writer_thread = threading.Thread(target=write_results_to_file, args=(output_dir,))
    writer_thread.start()
    
    
    total_references = len(target_list['PDB']) * len(ref_list['PDB'])

    # Create a single tqdm progress bar to track overall progress
    overall_progress = tqdm(total=total_references, desc='Processing references', leave=False)

    # Count the total number of references
    def process_reference(ref, target):
        resultat = run_align(pdb_folder, options, target, ref)
        results_queue.put(resultat)

    # Submit tasks for processing
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # tasks = [1,2,3,4,5,6,7,8,9]  # List of tasks to process
        # futures = [executor.submit(process_task, task) for task in tasks]
        futures = []
        count = 1  
        for target in target_list['PDB']:
            # Iterate over the elements of ref_list
            target_time = time.time()
            if exact_matching is True:
                print(count)
                ref = ref_list['PDB'][count]
                print(ref)
                count += 1
                if (target, ref) in to_skip:
                    print(f"{count} skipping {target} and {ref}")
                    continue
                futures.append(executor.submit(process_reference, ref))
            else:
                for ref in ref_list['PDB']:
                    count +=1
                    if (target,ref) in to_skip or target == ref:
                        #print(f"{count} skipping {target} and {ref}")
                        continue
                    futures.append(executor.submit(process_reference, ref, target))
                    with open(f'{output_dir}/missing.txt', 'a') as f:
                        f.write(f'{target},{ref}\n')
                                    
        for future in concurrent.futures.as_completed(futures):
            overall_progress.update(1)  # Update the overall progress bar
    # Notify the writer thread that all tasks are done
    all_tasks_done.set()
    writer_thread.join()  # Wait for the writer thread to finish
    
    
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

    result_df = pd.read_csv(f'{output_dir}/results.csv')
    selection_df = result_df.groupby(['PDB_ID_1', 'chain_1', 'PDB_ID_2', 'chain_2'], group_keys=False).apply(classify_match)

    # drop the rows that classified to type 0
    try:
        selection_df = selection_df[selection_df['type'] != 0]
    except:
        print('No matches found. Check the parameters and results.csv file.')
        return 0

    # for each chain_1 in the selection_df, delete the ':' and get length of the string
    if mode == '3':
        for row in selection_df.itertuples():
            # count how many there is string between ':' and ':' in the chain_1
            count_1 = [x for x in row.chain_1.split(':') if x]
            selection_df.loc[row.Index, 'chain_1'] = sum(1 for x in count_1)
            count_2 = [x for x in row.chain_2.split(':') if x]
            selection_df.loc[row.Index, 'chain_2'] = sum(1 for x in count_2)

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

if __name__ == "__main__":
    main()
