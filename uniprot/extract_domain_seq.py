#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function to extract protein domain sequences from UniProt IDs.
Outputs sequences in FASTA format with UniProt ID as the annotation.

Can process a file with UniProt IDs and domain names to extract multiple sequences.
"""

import requests
import json
import os
import sys
from typing import Dict, List, Optional, Union, Tuple


def extract_domain_sequence(uniprot_id: str, domain_name: Optional[str] = None) -> Union[Dict[str, str], str, None]:
    """
    Extract the sequence of a protein domain given a UniProt ID.
    
    Args:
        uniprot_id (str): The UniProt ID of the protein.
        domain_name (str, optional): The name of the specific domain to extract.
                                    If None, all domains will be returned.
    
    Returns:
        If domain_name is None: Dictionary mapping domain names to their sequences.
        If domain_name is provided: String containing the sequence of the specified domain.
        If no domains are found or an error occurs: None
    
    Example:
        >>> # Get all domains for a protein
        >>> domains = extract_domain_sequence("P00533")
        >>> # Get a specific domain
        >>> kinase_domain = extract_domain_sequence("P00533", "Protein kinase")
    """
    # UniProt API URL
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        # Make request to UniProt API
        response = requests.get(url, timeout=30)
        response.raise_for_status()  # Raise exception for HTTP errors
        
        # Parse JSON response
        data = response.json()
        
        # Extract full protein sequence
        full_sequence = data.get("sequence", {}).get("value", "")
        if not full_sequence:
            print(f"No sequence found for UniProt ID: {uniprot_id}")
            return None
        
        # Extract domain information
        domains = {}
        
        # Look for features that represent domains
        for feature in data.get("features", []):
            feature_type = feature.get("type", "")
            
            # Check if the feature is a domain
            if feature_type.lower() in ["domain", "region"]:
                domain_description = feature.get("description", "Unknown domain")
                
                # Get domain boundaries
                start = feature.get("location", {}).get("start", {}).get("value")
                end = feature.get("location", {}).get("end", {}).get("value")
                
                if start is not None and end is not None:
                    # Extract domain sequence (UniProt positions are 1-based)
                    domain_seq = full_sequence[start-1:end]
                    domains[domain_description] = domain_seq
        
        # If no domains found
        if not domains:
            print(f"No domain information found for UniProt ID: {uniprot_id}")
            return None
        
        # Return specific domain if requested
        if domain_name:
            # Try exact match first
            if domain_name in domains:
                return domains[domain_name]
            
            # Try case-insensitive partial match
            domain_name_lower = domain_name.lower()
            for key, seq in domains.items():
                if domain_name_lower in key.lower():
                    return seq
            
            print(f"Domain '{domain_name}' not found for UniProt ID: {uniprot_id}")
            print(f"Available domains: {', '.join(domains.keys())}")
            return None
        
        # Return all domains
        return domains
    
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from UniProt: {e}")
        return None
    except (json.JSONDecodeError, KeyError) as e:
        print(f"Error parsing UniProt response: {e}")
        return None


def write_fasta_file(uniprot_id: str, sequence: str, output_file: Optional[str] = None, domain_label: Optional[str] = None, mode: str = 'w') -> str:
    """
    Write a sequence to a FASTA format file.
    
    Args:
        uniprot_id (str): The UniProt ID to use as the sequence annotation.
        sequence (str): The sequence to write.
        output_file (str, optional): The path to the output file. If None, a default name will be used.
        domain_label (str, optional): The domain label to include in the header. If provided, 
                                     the header will be >uniprot_id_domain_label.
        mode (str): File opening mode ('w' for write, 'a' for append).
    
    Returns:
        str: The path to the created FASTA file.
    """
    if output_file is None:
        output_file = f"{uniprot_id}_domain.fasta"
    
    with open(output_file, mode) as f:
        # Format header based on whether domain_label is provided
        if domain_label:
            # Clean domain label for use in FASTA header (remove spaces, special chars)
            clean_label = domain_label.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
            header = f">{uniprot_id}_{clean_label}"
        else:
            header = f">{uniprot_id}"
        
        f.write(f"{header}\n")
        
        # Write sequence in lines of 60 characters (standard FASTA format)
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
    
    return output_file


def append_to_fasta_file(uniprot_id: str, sequence: str, output_file: str, domain_label: Optional[str] = None) -> str:
    """
    Append a sequence to an existing FASTA file.
    
    Args:
        uniprot_id (str): The UniProt ID to use as the sequence annotation.
        sequence (str): The sequence to write.
        output_file (str): The path to the output file.
        domain_label (str, optional): The domain label to include in the header.
    
    Returns:
        str: The path to the updated FASTA file.
    """
    return write_fasta_file(uniprot_id, sequence, output_file, domain_label, mode='a')


def get_kinase_domain_sequences(uniprot_id: str, output_file: Optional[str] = None) -> Dict[str, str]:
    """
    Extract all kinase domain sequences from a protein.
    
    Args:
        uniprot_id (str): The UniProt ID of the protein.
        output_file (str, optional): If provided, write the sequences to this file in FASTA format.
    
    Returns:
        Dict[str, str]: Dictionary mapping domain names to their sequences.
                       Empty dictionary if no kinase domains are found.
    
    Example:
        >>> kinase_domains = get_kinase_domain_sequences("P00533")
        >>> # Write to FASTA file
        >>> kinase_domains = get_kinase_domain_sequences("P00533", "P00533_kinase.fasta")
    """
    # Common names for kinase domains in UniProt
    kinase_domain_names = [
        "protein kinase", 
        "kinase", 
        "catalytic domain", 
        "tyrosine-protein kinase", 
        "serine/threonine-protein kinase"
    ]
    
    # Get all domains
    all_domains = extract_domain_sequence(uniprot_id)
    
    if not all_domains:
        return {}
    
    # Find all kinase domains
    kinase_domains = {}
    
    for domain_name, sequence in all_domains.items():
        for kinase_name in kinase_domain_names:
            if kinase_name.lower() in domain_name.lower():
                kinase_domains[domain_name] = sequence
                break
    
    # Write to FASTA file if requested
    if output_file is not None and kinase_domains:
        first_domain = True
        for domain_name, sequence in kinase_domains.items():
            if first_domain:
                write_fasta_file(uniprot_id, sequence, output_file, domain_name)
                first_domain = False
            else:
                append_to_fasta_file(uniprot_id, sequence, output_file, domain_name)
    
    # If no specific kinase domain found, print available domains
    if not kinase_domains:
        print(f"No kinase domains found for {uniprot_id}. Available domains:")
        for domain_name in all_domains.keys():
            print(f"  - {domain_name}")
    
    return kinase_domains


def get_kinase_domain_sequence(uniprot_id: str, output_file: Optional[str] = None) -> Optional[str]:
    """
    Extract the first kinase domain sequence from a protein (for backward compatibility).
    
    Args:
        uniprot_id (str): The UniProt ID of the protein.
        output_file (str, optional): If provided, write the sequence to this file in FASTA format.
    
    Returns:
        str or None: The sequence of the first kinase domain if found, None otherwise.
    """
    kinase_domains = get_kinase_domain_sequences(uniprot_id, output_file)
    
    if not kinase_domains:
        return None
    
    # Return the first kinase domain found
    domain_name = next(iter(kinase_domains))
    return kinase_domains[domain_name]


def process_uniprot_file(input_file: str, output_file: str, default_domain: str = "kinase") -> Dict[str, Dict[str, str]]:
    """
    Process a file containing UniProt IDs and domain names.
    
    The input file should have one entry per line in the format:
    UNIPROT_ID [DOMAIN_NAME]
    
    If DOMAIN_NAME is not provided, the default domain will be used.
    
    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output FASTA file.
        default_domain (str): Default domain to use if not specified.
        
    Returns:
        Dict[str, Dict[str, str]]: Dictionary mapping UniProt IDs to their domain sequences.
    """
    results = {}
    success_count = 0
    failure_count = 0
    first_entry = True
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return {}
    
    # Read input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Process each line
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        
        # Parse line
        parts = line.split()
        uniprot_id = parts[0]
        domain_name = " ".join(parts[1:]) if len(parts) > 1 else default_domain
        
        print(f"Processing {i+1}/{len(lines)}: {uniprot_id} (domain: {domain_name})")
        
        # Extract sequences
        if domain_name.lower() == "kinase":
            # Get all kinase domains
            domains = get_kinase_domain_sequences(uniprot_id)
            
            if domains:
                # Write all domains to output file
                for domain_label, sequence in domains.items():
                    if first_entry:
                        write_fasta_file(uniprot_id, sequence, output_file, domain_label)
                        first_entry = False
                    else:
                        append_to_fasta_file(uniprot_id, sequence, output_file, domain_label)
                
                results[uniprot_id] = domains
                success_count += 1
            else:
                print(f"  Warning: No kinase domains found for {uniprot_id}")
                failure_count += 1
        else:
            # Get specific domain
            sequence = extract_domain_sequence(uniprot_id, domain_name)
            
            if sequence:
                if first_entry:
                    write_fasta_file(uniprot_id, sequence, output_file, domain_name)
                    first_entry = False
                else:
                    append_to_fasta_file(uniprot_id, sequence, output_file, domain_name)
                
                results[uniprot_id] = {domain_name: sequence}
                success_count += 1
            else:
                print(f"  Warning: No sequence found for {uniprot_id} with domain '{domain_name}'")
                failure_count += 1
    
    # Count total domains extracted
    total_domains = sum(len(domains) for domains in results.values())
    
    print(f"\nProcessing complete: {success_count} proteins processed, {total_domains} domains extracted, {failure_count} failed")
    print(f"Results written to: {output_file}")
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Extract protein domain sequences from UniProt IDs")
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(dest="mode", help="Operation mode")
    
    # Single protein mode
    single_parser = subparsers.add_parser("single", help="Extract domain from a single protein")
    single_parser.add_argument("uniprot_id", help="UniProt ID of the protein")
    single_parser.add_argument("--domain", "-d", default="kinase", help="Domain name to extract (default: kinase)")
    single_parser.add_argument("--output", "-o", help="Output FASTA file path (default: <uniprot_id>_domain.fasta)")
    single_parser.add_argument("--list", "-l", action="store_true", help="List all available domains")
    
    # Batch mode
    batch_parser = subparsers.add_parser("batch", help="Extract domains from multiple proteins listed in a file")
    batch_parser.add_argument("input_file", help="Input file with UniProt IDs and optional domain names")
    batch_parser.add_argument("--output", "-o", default="domain_sequences.fasta", 
                             help="Output FASTA file path (default: domain_sequences.fasta)")
    batch_parser.add_argument("--default-domain", "-d", default="kinase", 
                             help="Default domain to use when not specified (default: kinase)")
    
    args = parser.parse_args()
    
    # If no mode specified, default to single mode for backward compatibility
    if args.mode is None:
        if len(sys.argv) > 1 and not sys.argv[1].startswith('-'):
            args.mode = "single"
            args.uniprot_id = sys.argv[1]
            args.domain = "kinase"
            args.output = None
            args.list = False
            
            # Check for additional arguments
            if len(sys.argv) > 2:
                if not sys.argv[2].startswith('-'):
                    args.domain = sys.argv[2]
                elif sys.argv[2] in ['-l', '--list']:
                    args.list = True
                elif sys.argv[2] in ['-o', '--output'] and len(sys.argv) > 3:
                    args.output = sys.argv[3]
        else:
            parser.print_help()
            sys.exit(1)
    
    # Process based on mode
    if args.mode == "single":
        if args.list:
            # List all domains
            domains = extract_domain_sequence(args.uniprot_id)
            if domains:
                print(f"Available domains for {args.uniprot_id}:")
                for name in domains.keys():
                    print(f"  - {name}")
            else:
                print(f"No domains found for {args.uniprot_id}")
        elif args.domain.lower() == "kinase":
            # Get all kinase domains
            domains = get_kinase_domain_sequences(args.uniprot_id, args.output)
            if domains:
                output_file = args.output if args.output else f"{args.uniprot_id}_domain.fasta"
                print(f"Found {len(domains)} kinase domains for {args.uniprot_id}, written to {output_file}")
                for domain_name, sequence in domains.items():
                    print(f"  - {domain_name}: {sequence[:50]}..." if len(sequence) > 50 else f"  - {domain_name}: {sequence}")
            else:
                print(f"No kinase domains found for {args.uniprot_id}")
        else:
            # Get specific domain
            result = extract_domain_sequence(args.uniprot_id, args.domain)
            if result:
                if args.output:
                    output_file = write_fasta_file(args.uniprot_id, result, args.output, args.domain)
                    print(f"Domain '{args.domain}' sequence written to {output_file}")
                else:
                    output_file = write_fasta_file(args.uniprot_id, result, domain_label=args.domain)
                    print(f"Domain '{args.domain}' sequence written to {output_file}")
                print(f"Sequence: {result[:50]}..." if len(result) > 50 else f"Sequence: {result}")
            else:
                print(f"Domain '{args.domain}' not found for {args.uniprot_id}")
    
    elif args.mode == "batch":
        # Process batch file
        process_uniprot_file(args.input_file, args.output, args.default_domain)

## Example usage
# python extract_domain_seq.py single P00533 --domain "protein kinase" --output P00533_kinase.fasta
# python extract_domain_seq.py batch uniprot_ids.txt --output domain_sequences.fasta --default-domain "kinase"