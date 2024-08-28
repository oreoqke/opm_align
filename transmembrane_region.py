import requests
import sys
import time

def graphQL_query(pdb_ids):
    """
    This function queries the RCSB GraphQL API to get the uniprot ids 
    for all the pdb ids in the list at once instead of querying them one by one.
    It can be reverted to querying a different API if needed, but it will slow down the code.
    """
    
    # Template for the GraphQL query
    query = """
    query getSubunitsWithUniProtAccession($ids: [String!]!) {
        entries(entry_ids: $ids) {
            rcsb_id
            rcsb_entry_container_identifiers {
            polymer_entity_ids
            }
            polymer_entities {
            entity_poly {
                pdbx_strand_id
            }
            rcsb_polymer_entity_container_identifiers {
                uniprot_ids
            }
            }
        }
    }
    """
    
    variables = {
        "ids": pdb_ids
    }
    url = "https://data.rcsb.org/graphql"
    payload = {
        "query": query,
        "variables": variables
    }
    response = requests.post(url, json=payload)
    
    # if the response is not successful the data will not be returned
    if response.status_code == 200:
        data = response.json()
        print(data)
        return data
    

def fetch_uniprot_entry(uniprot_id):
    url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.xml"
    response = requests.get(url)
    try:
        response.raise_for_status()
        return response.text
    except requests.exceptions.HTTPError as error:
        with open("error_log", "a") as f:
            f.write(f"Error occurred when fetching UniProt entry: {error}\n")
            
def parse_uniprot_entry(uniprot_xml):
    from xml.etree import ElementTree as ET
    root = ET.fromstring(uniprot_xml)
    subcellular_location = []
    keywords = []
    
    for comment in root.findall(".//{http://uniprot.org/uniprot}comment[@type='subcellular location']"):
        location = comment.find("{http://uniprot.org/uniprot}subcellularLocation/{http://uniprot.org/uniprot}location")
        if location is not None and location.text is not None:
            subcellular_location.append(location.text.lower())
    
    for keyword in root.findall(".//{http://uniprot.org/uniprot}keyword"):
        if keyword.text is not None:
            keywords.append(keyword.text.lower())
    
    return subcellular_location, keywords

def is_membrane(subcellular_locations, keywords):
    # Miight need to add more terms to this list later
    transmembrane_terms = ["membrane", "cell outer membrane", "transmembrane"]
    for term in transmembrane_terms:
        if any(term in location for location in subcellular_locations) or any(term in keyword for keyword in keywords):
               return True
    return False
        

def get_uniprot_entry(accession):
    """Fetch UniProt entry data suing UniProt API."""
    
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{accession}"
    print("requesting: "+accession)
    response = requests.get(url)
    try:
        response.raise_for_status()
        return response.json()
    # catch the error if the response is not successful
    # sometimes the pdb might not have a valid uniprot code
    except requests.exceptions.HTTPError as error:
        with open("error_log", "a") as f:
            f.write(f"Error occurred when fetching UniProt entry: {error}\n")
    # Sometimes the query might time out. 
    # To prevent the code from crashing retry the request in 5 seconds
    except requests.exceptions.RequestException as error:
        with open("error_log", "a") as f:
            f.write(f"Request error occurred when fetching UniProt entry: {error}\n")
        time.sleep(5)
        get_uniprot_entry(accession)
    
def extract_transmembrane_regions(uniprot_entry):
    """Extract transmembrane regions from UniProt entry. This gets the features 
    information and looks for any transmembrane labels."""
    transmembrane_regions = []
    for feature in uniprot_entry.get('features', []):
        if feature['type'] == 'TRANSMEM':
            transmembrane_regions.append({"begin": feature["begin"],
                                          "end": feature["end"]})
    
    return transmembrane_regions


def longest_sub(pdb_ids):
    # Read in the list of pdb ids from a file
    # with open("list_pdb", "r") as f:
    #     pdb_ids = f.read().splitlines()
    pdb_uniprot = {}
    
    # pdb_ids = ["8A9Y"]
    
    # Retrieve the uniprot codes for the pdb ids, including polymer entity ids
    data = graphQL_query(pdb_ids)['data']['entries']
    
    for entry in data:
        entities = {}
        for i,entt in enumerate(entry['polymer_entities']):
            entities[i+1] = { "uniprot": entt['rcsb_polymer_entity_container_identifiers']['uniprot_ids'],
                            "chains": entt['entity_poly']['pdbx_strand_id'],
                            "is_membrane": False,
                            "num_transmembrane": 0,
                            "transmembrane_regions": []}
        pdb_uniprot[entry['rcsb_id']] = entities
        
    
    longest_subs = {}
    
    # Retrieve the uniprot data for each uniprot code to determine whether a part of 
    # the chain is in membrane
    for key, value in pdb_uniprot.items():
        # this will loop over the polymer entities
        for i,entity in value.items():
            # some entries do not have uniprot codes, just skip them, assuming non-membrane
            if entity['uniprot'] == None:
                continue
            # Some entries have multiple uniprot codes, check if any of them have transmembrane regions
            # (I think multiple uniprots indicate chimeric protein but I am not sure, need to check)
            for uniprot in entity['uniprot']:                
                uniprot_data = get_uniprot_entry(uniprot)
                # print(uniprot_data)
                if uniprot_data is None:
                    continue
                transmembrane_regions = extract_transmembrane_regions(uniprot_data)
                # print(transmembrane_regions)
                if transmembrane_regions:
                    pdb_uniprot[key][i]["is_membrane"] = True
                    pdb_uniprot[key][i]["num_transmembrane"] = len(transmembrane_regions)
                    pdb_uniprot[key][i]["transmembrane_regions"] = transmembrane_regions
                    break
                else:
                    uni_data = fetch_uniprot_entry(uniprot)
                    subcellular_location, keywords = parse_uniprot_entry(uni_data)
                    if is_membrane(subcellular_location, keywords):
                        print(f"Membrane protein found: {key}")
                        pdb_uniprot[key][i]["is_membrane"] = True
                        break  
        
        max_num_transmem = 0
        max_uniprot_code = None
        max_subunit = None
        for k,i in pdb_uniprot[key].items():
            if i["num_transmembrane"] > max_num_transmem:
                max_num_transmem = i["num_transmembrane"]
                max_uniprot_code = i["uniprot"]
                max_subunit = k
        longest_subs[key] = {"uniprot": max_uniprot_code, "num_transmembrane": max_num_transmem, "subunit": max_subunit}
        print(f"Longest subunit for {key} is {max_uniprot_code} with {max_num_transmem} transmembrane regions in subunit {max_subunit}")
        
        # Write the uniprot codes to a file for debugging purposes
        with open("test_Uniprot", "a") as f:
            f.write(f"{key}: {value}")
            f.write("\n")
    
    return longest_subs

def main(pdb_id):
    #pdb_entry = get_pdb_entry(pdb_id)
    
    # for element in pdb_entry:
    #     print("----")
    #     print(element)
    #     print(pdb_entry[element])
    
    # try:
    #     entity_ids = pdb_entry['rcsb_entry_container_identifiers']['polymer_entity_ids']
    # except KeyError:
    #     # Skip if the polymer entity IDs are not available
    #     with open("error_log", "a") as f:
    #         f.write(f"Polymer entity IDs not found for PDB ID: {pdb_id}\n")
    #     return []
    
    # transmembrane_polypeptides = []
    
    # for entity_id in entity_ids:
    #     uniprot_entries = get_uniprot_accessions(pdb_id, entity_id)
    #     for element in uniprot_entries:
    #         print("----")
    #         print(element)  
    #         print(uniprot_entries[element])
        
    #     if uniprot_entries is None:
    #         continue
        
    #     for uniprot_entry in uniprot_entries:
    #         if uniprot_entry is None:
    #             # Skip if the UniProt accession is not available
    #             # TODO: Add this into a log file later
    #             continue
    #         print(uniprot_entry)
    #         accession = uniprot_entry['rcsb_id']
    #         uniprot_data = get_uniprot_entry(accession)
    #         transmembrane_regions = extract_transmembrane_regions(uniprot_data)
    #         if transmembrane_regions:
    #             transmembrane_polypeptides.append({
    #                 'entity_id': entity_id,
    #                 'uniprot_accession': accession,
    #                 "transmembrane_regions": transmembrane_regions
    #             })
    
    return transmembrane_polypeptides


if __name__ == "__main__":
    # Read in the list of pdb ids from a file
    with open("list_pdb", "r") as f:
        pdb_ids = f.read().splitlines()
    pdb_uniprot = {}
    
    # pdb_ids = ["8A9Y"]
    
    # Retrieve the uniprot codes for the pdb ids, including polymer entity ids
    data = graphQL_query(pdb_ids)['data']['entries']
    
    for entry in data:
        entities = {}
        for i,entt in enumerate(entry['polymer_entities']):
            entities[i+1] = { "uniprot": entt['rcsb_polymer_entity_container_identifiers']['uniprot_ids'],
                            "chains": entt['entity_poly']['pdbx_strand_id'],
                            "is_membrane": False,
                            "num_transmembrane": 0,
                            "transmembrane_regions": []}
        pdb_uniprot[entry['rcsb_id']] = entities
        
    
    longest_subs = {}
    
    # Retrieve the uniprot data for each uniprot code to determine whether a part of 
    # the chain is in membrane
    for key, value in pdb_uniprot.items():
        # this will loop over the polymer entities
        for i,entity in value.items():
            # some entries do not have uniprot codes, just skip them, assuming non-membrane
            if entity['uniprot'] == None:
                continue
            # Some entries have multiple uniprot codes, check if any of them have transmembrane regions
            # (I think multiple uniprots indicate chimeric protein but I am not sure, need to check)
            for uniprot in entity['uniprot']:                
                uniprot_data = get_uniprot_entry(uniprot)
                # print(uniprot_data)
                if uniprot_data is None:
                    continue
                transmembrane_regions = extract_transmembrane_regions(uniprot_data)
                # print(transmembrane_regions)
                if transmembrane_regions:
                    pdb_uniprot[key][i]["is_membrane"] = True
                    pdb_uniprot[key][i]["num_transmembrane"] = len(transmembrane_regions)
                    pdb_uniprot[key][i]["transmembrane_regions"] = transmembrane_regions
                    break
                else:
                    uni_data = fetch_uniprot_entry(uniprot)
                    subcellular_location, keywords = parse_uniprot_entry(uni_data)
                    if is_membrane(subcellular_location, keywords):
                        print(f"Membrane protein found: {key}")
                        pdb_uniprot[key][i]["is_membrane"] = True
                        break  
        
        max_num_transmem = 0
        max_uniprot_code = None
        max_subunit = None
        for k,i in pdb_uniprot[key].items():
            if i["num_transmembrane"] > max_num_transmem:
                max_num_transmem = i["num_transmembrane"]
                max_uniprot_code = i["uniprot"]
                max_subunit = k
        longest_subs[key] = {"uniprot": max_uniprot_code, "num_transmembrane": max_num_transmem, "subunit": max_subunit}
        print(f"Longest subunit for {key} is {max_uniprot_code} with {max_num_transmem} transmembrane regions in subunit {max_subunit}")
        
        # Write the uniprot codes to a file for debugging purposes
        with open("test_Uniprot", "a") as f:
            f.write(f"{key}: {value}")
            f.write("\n")
        with open("longest_subs", "a") as f:
            f.write(f"{key}: {longest_subs[key]}")
            f.write("\n")
            
        

            

    