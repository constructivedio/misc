import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Remove a specific sequence from a SnapGene file if found exactly once. "
                    "The output is saved in GenBank format."
    )
    parser.add_argument(
        "--snapgene", 
        required=True, 
        help="Path to the input SnapGene file (.dna)."
    )
    parser.add_argument(
        "--sequence", 
        required=True, 
        help="The DNA sequence string to search for and remove."
    )
    parser.add_argument(
        "--output_file", 
        required=True, 
        help="Path for the output modified file (will be in GenBank format)."
    )

    args = parser.parse_args()

    if not args.sequence:
        print("Error: The --sequence argument cannot be empty.", file=sys.stderr)
        sys.exit(1)

    try:
        # Read the SnapGene file. SeqIO.read expects exactly one record.
        record = SeqIO.read(args.snapgene, "snapgene")
    except FileNotFoundError:
        print(f"Error: Input file '{args.snapgene}' not found.", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error reading SnapGene file '{args.snapgene}': {e}. "
              "Please ensure it's a valid SnapGene file containing exactly one sequence record.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading '{args.snapgene}': {e}", file=sys.stderr)
        sys.exit(1)

    sequence_to_remove = args.sequence.lower()

    # Count occurrences of the sequence to remove
    count = record.lower().count(sequence_to_remove)

    if count == 1:
        # Find the start index of the sequence
        start_index = record.lower().seq.find(sequence_to_remove)
        
        # Create the new sequence by removing the target sequence
        modified_record = record[0:start_index] + record[start_index + len(sequence_to_remove):]

        # make the modified record circular
        modified_record.annotations["topology"] = "circular"
        modified_record.annotations["molecule_type"] = "DNA"
        
        try:
            # Write the modified record to the output file in GenBank format.
            # Biopython does not write native SnapGene files.
            # SnapGene can open GenBank files.
            SeqIO.write(modified_record, args.output_file, "genbank")
            print(f"Sequence '{sequence_to_remove}' found once and removed. "
                  f"Modified sequence saved to '{args.output_file}' (GenBank format).")
        except Exception as e:
            print(f"Error writing output file '{args.output_file}': {e}", file=sys.stderr)
            sys.exit(1)
            
    elif count == 0:
        print(f"Sequence '{sequence_to_remove}' not found in '{args.snapgene}'. No changes made.")
    else:  # count > 1
        print(f"Sequence '{sequence_to_remove}' found {count} times in '{args.snapgene}'. "
              "No changes made as it was not found exactly once.")

if __name__ == "__main__":
    main()