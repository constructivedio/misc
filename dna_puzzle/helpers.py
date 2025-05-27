import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation, SimpleLocation
from Bio import SeqIO
from copy import deepcopy
import os

def read_order_csv(order_csv_path):
    """
    Read the order CSV file and return a dictionary mapping filenames to their ordered sources.
    
    Args:
        order_csv_path (str): Path to the CSV file containing the element order
        
    Returns:
        dict: Dictionary with filenames as keys and lists of sources as values
    """
    # Read CSV without headers since columns are unnamed
    df = pd.read_csv(order_csv_path, header=None, na_filter=False)

    # Create dictionary of filename -> sources
    order_dict = {}
    for idx, row in df.iterrows():
        # Remove filename and any NaN values to get clean list of sources
        sources = []
        raw_sources = row[1:].dropna().tolist()
        for source in raw_sources:
            this_source = {"value":source.split(":")[0]}
            if ":annotate" in source:
                this_source["annotate"] = True
            if ":previous_overlap" in source:
                this_source["previous_overlap"] = True
            sources.append(this_source)
        order_dict[row[0]] = sources
        
    return order_dict

def parse_dna_info(dna_csv_path):
    """
    Parse DNA information CSV file and extract sequence information.
    
    Args:
        dna_csv_path (str): Path to the CSV file containing DNA information
        
    Returns:
        dict: Dictionary with element IDs as keys and tuples of (type, info) as values
            where type is one of: 'sequence', 'file', 'feature', 'region'
            and info contains relevant data for each type
    """
    df = pd.read_csv(dna_csv_path, header=None)
    
    if len(df.columns) != 2:
        raise ValueError(f"DNA info CSV must have exactly 2 columns. Found {len(df.columns)} columns.")
    
    dna_dict = {}
    
    for _, row in df.iterrows():
        element_id = row[0]
        dna_info = row[1]
        
        if ':' not in dna_info and not dna_info.endswith('.dna'):
            # Raw DNA sequence
            dna_dict[element_id] = ('raw_sequence', dna_info)
        elif dna_info.endswith('.dna') and ':' not in dna_info:
            # Full SnapGene file
            dna_dict[element_id] = ('file', dna_info)
        elif dna_info.count(':') == 1 and dna_info.endswith('.dna'):
            # Feature from SnapGene file
            feature_id, snapgene_file = dna_info.split(':')
            dna_dict[element_id] = ('snapgene_feature', (feature_id, snapgene_file))
        elif dna_info.count(':') == 1 and dna_info.count('-') == 1:
            # Region from SnapGene file
            snapgene_file, coordinates = dna_info.split(':')
            start, end = coordinates.split('-')
            dna_dict[element_id] = ('snapgene_region', (snapgene_file, int(start), int(end)))
            
    return dna_dict

def create_record_from_raw_sequence(element_id, sequence_info):
    """
    Create a SeqFeature object from a raw sequence.
    
    Args:
        element_id (str): ID of the element to use as feature label
        sequence_info (tuple): Tuple containing ('raw_sequence', sequence)
        
    Returns:
        SeqRecord: A SeqRecord with a misc feature with the element_id as label
    """
    record = SeqIO.SeqRecord(sequence_info[1])
    feature = SeqFeature(
        location=FeatureLocation(0, len(sequence_info[1])),
        type="misc_feature",
        id=element_id,
        qualifiers={"label": [element_id]}
    )
    
    record.features.append(feature)
    
    return record

def create_record_from_snapgene_file(element_id, file_path):
    """
    Create a SeqRecord object from a SnapGene file.
    
    Args:
        element_id (str): ID of the element to use as feature label
        file_path (str): Path to the SnapGene file
        
    Returns:
        SeqRecord: A SeqRecord with misc feature on its full length with the element_id as label
    """
    
    records = list(SeqIO.parse(file_path, "snapgene"))
    if len(records) != 1:
        raise ValueError(f"Expected exactly one record in {file_path}, found {len(records)}")
        
    record = records[0]
    feature = SeqFeature(
        location=FeatureLocation(0, len(record.seq)),
        type="misc_feature",
        id=element_id,
        qualifiers={"label": [element_id]}
    )
    record.features.append(feature)
    
    new_record = record[feature.location.start:feature.location.end]

    return new_record

def create_record_from_snapgene_region(file_path, start, end):
    """
    Create a SeqRecord object from a region in a SnapGene file.
    
    Args:
        file_path (str): Path to the SnapGene file
        start (int): Start coordinate of the region
        end (int): End coordinate of the region
        If start>end (circular plasmid case), the region will be concatenated from start to the end of the sequence and from the beginning to end
        
    Returns:
        SeqRecord: A SeqRecord with the region from start to end
    """
    records = list(SeqIO.parse(file_path, "snapgene"))
    if len(records) != 1:
        raise ValueError(f"Expected exactly one record in {file_path}, found {len(records)}")
        
    record = records[0]
    if start < 0 or end > len(record.seq):
        raise ValueError(f"Region coordinates ({start}, {end}) out of bounds for sequence length {len(record.seq)}")
    

    if start > end:
        chunk1 = deepcopy(record[start:])
        chunk2 = deepcopy(record[:end])
        chunk = chunk1 + chunk2
    else:
        chunk = record[start:end]
    return chunk

def create_record_from_searching_snapgene_with_feature_label(file_path, feature_name):
    """
    Searches all the features in a SnapGene file for a feature with the specified label and returns the
    whole region of that feature including all the features in the region.
    
    Args:
        file_path (str): Path to the SnapGene file
        feature_name (str): Name of the feature to search for
        
    Returns:
        SeqRecord: A SeqRecord with the feature matching the label and all other features in the region
        
    Raises:
        ValueError: If no feature or multiple features with the given name are found
    """
    records = list(SeqIO.parse(file_path, "snapgene"))
    
    all_matching_features_records = []

    for record in records:
        matching_features = [{"record":record,"feature":f} for f in record.features if "label" in f.qualifiers and feature_name in f.qualifiers["label"]]
        all_matching_features_records.extend(matching_features)

    if len(all_matching_features_records) == 0:
        raise ValueError(f"No feature with label '{feature_name}' found in {file_path}")
    if len(all_matching_features_records) > 1:
        raise ValueError(f"Multiple features with label '{feature_name}' found in {file_path}")
    
    matching_feature_record = all_matching_features_records[0]
    recordoi = matching_feature_record["record"]
    featureoi = matching_feature_record["feature"]

    return(recordoi[featureoi.location.start:featureoi.location.end])

def dna_info_to_dict(dna_info_key, dna_info, source_files_path):
    """
    Convert a tuple of DNA information to a dictionary.
    
    Args:
        dna_info (tuple): Tuple containing (type, info)
        
    Returns:
        dict: Dictionary with 'type' and 'info' keys
    """
    if dna_info[0] == "raw_sequence":
        new_feature = create_record_from_raw_sequence(dna_info_key, dna_info)
    elif dna_info[0] == "file":
        new_feature = create_record_from_snapgene_file(dna_info_key, os.path.join(source_files_path,dna_info[1]))
    elif dna_info[0] == "snapgene_feature":
        new_feature = create_record_from_searching_snapgene_with_feature_label(os.path.join(source_files_path,dna_info[1][1]), dna_info[1][0])
    elif dna_info[0] == "snapgene_region":
        new_feature = create_record_from_snapgene_region(os.path.join(source_files_path,dna_info[1][0]), dna_info[1][1], dna_info[1][2])
    else:
        raise ValueError(f"Invalid DNA info type: {dna_info[0]}")
    return({dna_info_key: new_feature})
    
def parse_dna_info_file_to_feature_dict(dna_csv_path, source_files_path="."):
    """
    Parse DNA information CSV file and convert to a dictionary of SeqFeature objects.
    
    Args:
        dna_csv_path (str): Path to the CSV file containing DNA information
        
    Returns:
        dict: Dictionary with element IDs as keys and SeqFeature objects as values
    """
    dna_info = parse_dna_info(dna_csv_path)
    dna_info_dict = {}
    for key, value in dna_info.items():
        if key in dna_info_dict:
            raise ValueError(f"Duplicate element ID found: {key}. Element IDs must be unique.")
        dna_info_dict.update(dna_info_to_dict(key, value, source_files_path))
    return dna_info_dict

def find_overlap_index_for_appending_sequence(original_dna, appending_dna):
    """
    Merge two DNA strings by finding overlapping regions and avoiding duplication.
    
    Args:
        original_dna (str): The first DNA string
        appending_dna (str): The second DNA string to append
        
    Returns:
        str: Merged DNA string with overlapping region included only once
    """
    # Edge cases
    if not original_dna:
        return appending_dna
    if not appending_dna:
        return original_dna
    # Find the longest overlap
    max_overlap = 0
    for i in range(min(len(original_dna), len(appending_dna)), 0, -1):
        if original_dna[-i:].lower() == appending_dna[:i].lower():
            max_overlap = i
            break
    return max_overlap

def find_and_annotate_sequence_on_seqrecord(record, sequence_to_annoate, feature_type="misc_features", feature_name=None):
    """
    Find a sequence in a SeqRecord and annotate it as a feature.
    
    Args:
        record (SeqRecord): The SeqRecord to search
        sequence_to_annoate (SeqRecord): The sequence to search for and annotate
        feature_type (str): The type of feature to create
        feature_name (str): The name of the feature
    
    Returns:
        SeqRecord: The SeqRecord with the annotated feature
    """
    if not feature_name:
        feature_name = "annotated_feature"
    element_to_add_start = record.seq.lower().find(sequence_to_annoate.seq.lower())
    element_to_add_feature = SeqFeature(
        location=FeatureLocation(element_to_add_start, element_to_add_start+len(sequence_to_annoate.seq)),
        type=feature_type,
        id=feature_name,
        qualifiers={"label": feature_name}
    )
    record.features.append(element_to_add_feature)
    return record

def merge_two_overlapping_seqrecords(original_record, appending_record, feature_name=None):
    """
    Merge two SeqRecord objects by finding overlapping regions and avoiding duplication.
    
    Args:
        original_record (SeqRecord): The first SeqRecord
        appending_record (SeqRecord): The second SeqRecord to append
        
    Returns:
        SeqRecord: Merged SeqRecord with overlapping region included only once
    """
    # Edge cases
    if not original_record:
        return appending_record
    if not appending_record:
        return original_record
    # Find the overlap index
    overlap_index = find_overlap_index_for_appending_sequence(str(original_record.seq), str(appending_record.seq))
    # Append the second record to the first
    merged_record = original_record + appending_record[overlap_index:]
    annotated_record = find_and_annotate_sequence_on_seqrecord(merged_record, appending_record, feature_name=feature_name)
    return annotated_record
    

def construct_order_of_elements_to_dna(order, source_dict, circular=False, molecule_type="DNA"):
    """Takes in the order of DNA elements and the source dictionary and constructs the DNA sequence

    Args:
        order (list): List of dictionaries, specifying the order and type of DNA elements to construct the SeqRecord object
                      Example: [{"value": "100k02_LSup"}, {"value": "100k21_chonks"}, {"value": "1A1_sC"}]
        sources_dict (dict): Dictionary with element IDs as keys and SeqRecord objects as values (output of parse_dna_info_file_to_feature_dict)
        
    Returns:
        SeqRecord: A SeqRecord with the constructed DNA sequence
    """
    annotate_elements = []
    constructed_dna_fragment = None
    for idx, element in enumerate(order):
        if element["value"] not in source_dict:
            raise ValueError(f"Element '{element}' not found in source CSV")
        if "annotate" in element and element["annotate"]:
            annotate_elements.append(element["value"])
            continue
        if constructed_dna_fragment is None:
            constructed_dna_fragment = source_dict[element["value"]]
        else:
            if "previous_overlap" in element and element["previous_overlap"]:
                constructed_dna_fragment = merge_two_overlapping_seqrecords(constructed_dna_fragment, source_dict[element["value"]], feature_name=element["value"])
            else:
                constructed_dna_fragment += source_dict[element["value"]]
        
        if annotate_elements != []:
            for element_to_annotate in annotate_elements:
                constructed_dna_fragment = find_and_annotate_sequence_on_seqrecord(constructed_dna_fragment, source_dict[element_to_annotate], feature_name=element_to_annotate)
        
        constructed_dna_fragment.annotations = {
            "molecule_type": molecule_type,
            "topology": "circular" if circular else "linear"
        }
    
    return constructed_dna_fragment

def write_seqrecord_to_genbank(seqrecord, outfile):
    """
    Write a SeqRecord object to a GenBank file.
    
    Args:
        seqrecord (SeqRecord): The SeqRecord to write
        outfile (str): Path to the output GenBank file
    """
    if not outfile.endswith('.gb') or outfile.endswith('.gbk'):
        outfile = f'{outfile}.gb'
    SeqIO.write(seqrecord, outfile, "genbank")

def write_seqrecord_to_file(seqrecord, outtype, output_dir, outfile):
    """Write a SeqRecord object to a file."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    full_output_file = os.path.join(output_dir, outfile)
    if outtype == 'genbank':
        write_seqrecord_to_genbank(seqrecord, full_output_file)
    else:
        raise ValueError(f"Output format '{outtype}' not supported")
    return(full_output_file)