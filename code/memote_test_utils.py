import memote.suite.api as memote_api
import cobra
import pandas as pd
import re
import os
import glob
import logging

# Function to find model file path
def find_model_file_path(model_dir):
    """Finds an XML or YML file in the model directory, prioritizing XML."""
    xml_files = glob.glob(os.path.join(model_dir, "*.xml"))
    yml_files = glob.glob(os.path.join(model_dir, "*.yml"))
    
    # Convert to forward slashes
    xml_files = [file.replace("\\", "/") for file in xml_files]  
    yml_files = [file.replace("\\", "/") for file in yml_files]
    
    # Prioritizes XML
    if xml_files:
        logging.info(f"Returned path to XML model file: {xml_files[0]}")
        return xml_files[0]
    elif yml_files:
        logging.info(f"Returned path to YML model file: {yml_files[0]}")
        return yml_files[0]
    else:
        logging.error("No XML or YML model file found in the directory.")
        raise FileNotFoundError("No XML or YML model file found in the directory.")

# Function to check and load TSV files if they exist
def load_annotation_tsv_files(model_dir):
    """Loads genes.tsv, mets.tsv, and rxns.tsv if all exist, otherwise raises an error."""
    required_tsv_files = ["genes.tsv", "mets.tsv", "rxns.tsv"]
    tsv_paths = {tsv: os.path.join(model_dir, tsv) for tsv in required_tsv_files}

    missing_files = [tsv for tsv, path in tsv_paths.items() if not os.path.exists(path)]
    if missing_files:
        logging.error(f"Missing required annotation file(s): {', '.join(missing_files)}")
        raise FileNotFoundError(f"Missing required annotation file(s): {', '.join(missing_files)}")

    logging.info("Annotation TSV information loaded.")
    return {tsv: pd.read_csv(path, sep="\t") for tsv, path in tsv_paths.items()}

# Function to parse annotation field
def parse_annotation_field(field_str):
    """Parses a structured annotation field (e.g., NAME, MIRIAM) into a dictionary."""
    annotation_dict = {}
    if pd.notna(field_str):
        for entry in field_str.split(";"):
            key_value = entry.split("/")
            if len(key_value) == 2:
                key, value = key_value
                annotation_dict.setdefault(key, []).append(value)
    return annotation_dict

# Function to append annotation information to the model
def append_annotations(model, genes_df, mets_df, rxns_df):
    """Adds annotation data to the model from TSV files."""
    # Ensure all genes have SBO:0000243 (standard annotation that every gene is expected to have as default)
    for gene in model.genes:
        gene.annotation.setdefault("sbo", "SBO:0000243")
        
    # Update gene annotations
    for _, row in genes_df.iterrows():
        gene_id_clean = row[genes_df.columns[0]].removeprefix("G_")  # Cleaned ID
        if gene_id_clean in model.genes:
            gene = model.genes.get_by_id(gene_id_clean)
            annotations = row.drop(genes_df.columns[0]).dropna().to_dict()
    
            # Expand NAME field
            if "NAME" in annotations:
                name_dict = parse_annotation_field(annotations.pop("NAME"))
                annotations.update(name_dict)
    
            gene.annotation.update(annotations)
    
    # Update metabolite annotations
    for _, row in mets_df.iterrows():
        met_id_clean = row[mets_df.columns[8]].removeprefix("M_")  # Clean ID
        if met_id_clean in model.metabolites:
            met = model.metabolites.get_by_id(met_id_clean)
            annotations = row.drop(mets_df.columns[8]).dropna().to_dict()
    
            # Expand MIRIAM field
            if "MIRIAM" in annotations:
                miriam_dict = parse_annotation_field(annotations.pop("MIRIAM"))
                annotations.update(miriam_dict)
    
            met.annotation.update(annotations)
    
    # Update reaction annotations
    for _, row in rxns_df.iterrows():
        rxn_id_clean = row[rxns_df.columns[1]].removeprefix("R_")  # Clean ID
        if rxn_id_clean in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id_clean)
            annotations = row.drop(rxns_df.columns[1]).dropna().to_dict()
    
            # Update EC-NUMBER key (required by MEMOTE)
            if "EC-NUMBER" in annotations:
                annotations["ec-code"] = annotations.pop("EC-NUMBER")
    
            # Expand MIRIAM field
            if "MIRIAM" in annotations:
                miriam_dict = parse_annotation_field(annotations.pop("MIRIAM"))
                annotations.update(miriam_dict)
            
            # Extract SBO term from MIRIAM if present
            miriam = annotations.get("MIRIAM", "")
            sbo_match = re.search(r"sbo/SBO:\d+", miriam)
            if sbo_match:
                sbo_term = sbo_match.group().split("/")[-1]  # Extract SBO code
                annotations["sbo"] = sbo_term
            
            rxn.annotation.update(annotations)
            
