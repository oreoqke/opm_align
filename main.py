import click
import os
import logging
import get_all_fasta
import get_all_pdbid
import subprocess
import pandas as pd
import align_structures

def align_sequences(output_dir, logger):
    # read in the data frame from the alignment
    file_path = os.path.join(output_dir, "result.csv")
    df = pd.read_csv(file_path, header=None, sep="\t")
    # I'm using aligment score as the similarity measure, so I know which 3d structures
    # are likely to be most similar. I will try to find the best match from any subunit to any subunit
    # of the target protein.
    df = df.map(lambda x: x.lower() if isinstance(x, str) else x)
    df = df.map(lambda x: x[:4] if isinstance(x, str) else x)
    # multisubunit proteins will result in duplicate matches
    # This is needed because we don't care about different subunits
    df = df.drop_duplicates()

    # For each target value, get the best match
    to_align = []
    targets = open(os.path.join(output_dir, "new_pdb.txt"), "r").read().splitlines()
    for target in targets:
        filtered_df = df[df[0] == target]
        filtered_df = filtered_df.sort_values(by=2, ascending=False)
        if filtered_df.empty:
            print(f"No match found for {target}")
            logger.warning(f"No match found for {target}")
            continue
        i = 0
        while filtered_df.shape[0] != i and filtered_df.iloc[i, 2] >= filtered_df.iloc[0, 2]-10:
            to_align.append(filtered_df.iloc[i])
            i += 1
        #print(target)
        #print(filtered_df)
    with open(os.path.join(output_dir, "to_align.txt"), "w") as f:
        for line in to_align:
            f.write(f"{line[0]} {line[1]} {line[2]}\n")


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
    command = f'blastp -query {new_fasta} -db {db_path} -outfmt "6 qseqid sseqid pident" -out {output_path}'
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