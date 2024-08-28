import click
import os
import logging
import get_all_fasta
import get_all_pdbid
import subprocess
import pandas as pd
import align_structures
import transmembrane_region as tr
'''The part that does USalign was made posssible by the generous contribution of Stanislav Cherepanov'''
def write_bad_match(output_dir, info):
    with open(os.path.join(output_dir, "bad_matches.txt"), "a") as f:
        f.write(f"{info}\n")
        
def criteria_for_usalign(output_dir, filtered_df, i):
    ''' accepts valid output directory, index i, and filtered data frame
    returns whether the pair of pdbids should be used for structure alignment'''
    # if the best sequence alignment is less than 15% don't bother aligning
    # write down this entry
    if filtered_df.iloc[0, 2] < 15:
        write_bad_match(output_dir, filtered_df.iloc[0])
        return False
    # if the best sequence alignment is better than 30% align only the best value
    elif filtered_df.iloc[0, 2] > 30:
        return (filtered_df.iloc[i, 2] == filtered_df.iloc[0, 2])
    # if the best value is between 15 and 30% align the best values within 5% range
    else:
        return filtered_df.iloc[i, 2] >= filtered_df.iloc[0, 2]-5
    

def align_sequences(output_dir, logger):
    # read in the data frame from the alignment
    file_path = os.path.join(output_dir, "result.csv")
    df = pd.read_csv(file_path, header=None, sep="\t")
    # I'm using aligment score as the similarity measure, so I know which 3d structures
    # are likely to be most similar. I will try to find the best match from any subunit to any subunit
    # of the target protein.
    df = df.map(lambda x: x.upper() if isinstance(x, str) else x)
    
    # cuts off everything except for the pdb id and the subunit
    df = df.map(lambda x: x.split('|')[0] if isinstance(x, str) else x)

    # For each target value, get the best match
    to_align = []
    targets = open(os.path.join(output_dir, "new_pdb.txt"), "r").read().splitlines()
    
    # Convert all targets to uppercase to match the dictionary
    upper_targets = [target.upper() for target in targets]
    
    # Get the transmembrane subunit with the most transmembrane regions
    longest_subs = tr.longest_sub(targets)        
    
    # Filter the data frame for the best match for each target
    for target in upper_targets:
        # If the target does not have any transmembrane regions, we only filter by the 
        # sequence alignment score
        # TODO: I am not sure if this would ever be the case. Need to be tested
        if target not in longest_subs:
            filtered_df = df.map(lambda x: x[:4] if isinstance(x, str) else x)
            filtered_df = filtered_df[filtered_df[0] == target]
        # Take the subunit with the most transmembrane regions and drop other ones from the data frame
        else:
            memb_subunit = longest_subs[target]["subunit"]
            if memb_subunit is not None:
                filtered_df = df[df[0] == f"{target}_{memb_subunit}"]
                filtered_df = filtered_df.map(lambda x: x[:4] if isinstance(x, str) else x)
                print(f"Filtering for {target}_{memb_subunit}")
                print(filtered_df)
            # If the target does not have any transmembrane regions, we only filter by the 
            # sequence alignment score
            else:
                filtered_df = df.map(lambda x: x[:4] if isinstance(x, str) else x)
                print(filtered_df)
                filtered_df = filtered_df[filtered_df[0] == target]
                print(f"Filtering for {target}_{memb_subunit}")
                print(filtered_df)
                
        # Drop duplicates if there are any
        # This might happen if we don't have any information about the subunit
        filtered_df = filtered_df.drop_duplicates()
        filtered_df = filtered_df.sort_values(by=2, ascending=False)
        
        # Write a helpful message that no matches were found
        if filtered_df.empty:
            print(f"No match found for {target} by BLAST")
            logger.warning(f"No match found for {target} by BLAST")
            write_bad_match(output_dir, f"{target} No match found by BLAST")
            continue
        i = 0
        # This loop will determine the pairs for which to run the alignment
        while filtered_df.shape[0] !=i and criteria_for_usalign(output_dir, filtered_df, i):
            print(filtered_df.iloc[i])
            to_align.append(filtered_df.iloc[i])
            i += 1
    
    with open(os.path.join(output_dir, "to_align.txt"), "w") as f:
        for line in to_align:
            f.write(f"{line[0]} {line[1]} {line[2]} {line[3]}\n")
           


@click.command()
# Specify the input for the new pdb ids, the default is new_pdb.txt
@click.option("--input", "-i", default="new_pdb.txt", 
              type=click.Path(exists=True), help="Specify input file, default is new_pdb.txt", required=False)
@click.option("--output", "-o", default="result.csv", 
              type=click.Path(), help="Specify output file, default is result.csv", required=False)
@click.option("--output-dir", "-d", default="results",
              type=click.Path(), help="Specify output directory, default is results", required=False)
def main(input, output, output_dir):
    """This funciton accepts a list of pdb ids and finds the best match for each one from the database.
    It first runs the alignment by sequence and then by structure for the best match.
    It then writes the results into a results directory."""
      
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        click.echo(f"Directory '{output_dir}' created")
        
    # Configure logging
    logging.basicConfig(filename= os.path.join(output_dir, "find_best_match.log"), level=logging.DEBUG)
    # Create a file handler for warnings
    warning_handler = logging.FileHandler(filename= os.path.join(output_dir, "warnings.log"))
    warning_handler.setLevel(logging.WARNING)
    warning_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger('').addHandler(warning_handler)
    
    logger = logging.getLogger()
    
    output_path = os.path.join(output_dir, output)
    click.echo(f"Output file path: {output_path}")
    
    # go into get_all_fasta.py and download all fasta files for
    # the new pdb ids
    new_fasta = os.path.join(output_dir, "new_fasta.txt")
    get_all_fasta.get_fasta(input, new_fasta, logger)
    # Copy the input file to the output directory
    cmd = f"cp {input} {os.path.join(output_dir, 'new_pdb.txt')}"
    subprocess.run(cmd, shell=True)
    
    # if you don't provide a reference list, it will take all the ids from
    # the existing opm database
    opm_pdb_ids = os.path.join(output_dir, "opm_pdbid.txt")
    if not os.path.exists(os.path.join(output_dir, "opm_pdbid.txt")):
        # Get the updataed list of pdb ids from the database
        get_all_pdbid.get_all_pdbid(output_file=opm_pdb_ids)
    
    
    # Now check if the db fasta files are already in the output
    # directory. This means we have already downloaded them
    # in previous runs. If not, download them
    db_fasta = os.path.join(output_dir, "db_fasta.txt")
    if not os.path.exists(db_fasta):
        get_all_fasta.get_fasta(opm_pdb_ids, db_fasta, logger)
    
    # Now we can construct the database for the alignment
    db_path = os.path.join(output_dir, "db")
    command = f"makeblastdb -in {db_fasta} -dbtype prot -out {db_path}"
    res = subprocess.run(command, shell=True, capture_output=True, text=True)
    # Log everything and capture any errors
    print(res.stdout)
    logger.info(res.stdout)
    if res.stderr:
        print(res.stderr)
        logger.warning(res.stderr)
        
    # Now we can run the sequence alignment
    # for more information on the command see the blastp documentation
    # https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/
    command = f'blastp -query {new_fasta} -db {db_path} -outfmt "6 qseqid sseqid pident length" -out {output_path}'
    res = subprocess.run(command, shell=True, capture_output=True, text=True)
    # Print logs and capture any errors
    print(res.stdout)
    logger.info(res.stdout)
    if res.stderr:
        print(res.stderr)
        logger.warning(res.stderr)
    
    # Finally, we can run the structure alignment
    align_sequences(output_dir, logger)
    align_structures.align_structures(output_dir, logger)
    
    
    
    
    
if __name__ == "__main__":
    main()