import psycopg2

def get_all_pdbid(dbname="lomize-opm", output_file="opm_pdbid.txt"):
    # Connect to the database
    conn = psycopg2.connect(
        dbname=dbname,
    )

    # Create a cursor
    cursor = conn.cursor()
    query = "SELECT pdbid FROM primary_structures"
    cursor.execute(query)
    results = cursor.fetchall()
    # Close the cursor and connection
    cursor.close()
    conn.close()
    # Write the results to a file
    with open(output_file, 'w') as f:
        for row in results:
            f.write(str(row[0]) + '\n')
