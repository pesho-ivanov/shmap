#!/usr/bin/env python3
import re
import sys
import argparse
import subprocess
import os


def create_bed_from_fasta(fasta_file, bed_file):
    """Create a BED file containing coordinates for all reads in the FASTA file."""
    total_descriptions = 0
    sequences = {}  # Store sequences for later reconstruction
    
    # Handle stdin if fasta_file is '-'
    fasta = sys.stdin if fasta_file == '-' else open(fasta_file)
    with fasta, open(bed_file, "w") as bed:
        current_id = None
        for line in fasta:
            if line.startswith(">"):
                total_descriptions += 1
                match = re.match(r">(.*?)!(chr\S+?)_(\S+?)!(\d+?)!(\d+?)!(\S+)", line)
                if match:
                    prefix, chrom, haplotype, start, end, strand = match.groups()
                    current_id = f"{total_descriptions}"  # Use counter as unique ID
                    # Store original header info for reconstruction
                    sequences[current_id] = {
                        'header': line.strip(),
                        'sequence': '',
                        'prefix': prefix,
                        'strand': strand
                    }
                    # Write BED entry with unique ID as name
                    bed.write(f"{chrom}_{haplotype}\t{start}\t{end}\t{current_id}\t0\t{strand}\n")
            elif current_id is not None:
                sequences[current_id]['sequence'] += line.strip()
    
    return sequences, total_descriptions


def run_liftover(input_bed, chain_file, output_bed, unmapped_bed):
    """Run liftOver on the entire BED file at once."""
    try:
        # Get absolute path to liftOver relative to script directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)  # Project root is one level up from evals/
        liftover_path = os.path.join(project_root, "ext", "liftOver")
        
        # Make liftover executable
        os.chmod(liftover_path, 0o755)
        
        # Get absolute paths for all files
        input_bed = os.path.abspath(input_bed)
        chain_file = os.path.abspath(chain_file)
        output_bed = os.path.abspath(output_bed)
        unmapped_bed = os.path.abspath(unmapped_bed)
        
        cmd = [liftover_path, input_bed, chain_file, output_bed, unmapped_bed]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"liftOver stdout: {result.stdout}", file=sys.stderr)
        print(f"liftOver stderr: {result.stderr}", file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running liftOver: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Unexpected error running liftOver: {e}", file=sys.stderr)
        return False


def create_fasta_from_lifted_bed(lifted_bed, unmapped_bed, sequences):
    """Convert lifted BED coordinates back to FASTA format."""
    mapped_ids = set()
    unmapped_count = 0
    
    # Process successfully lifted coordinates
    try:
        with open(lifted_bed) as f:
            for line in f:
                chrom, start, end, read_id, _, strand = line.strip().split()[:6]
                if read_id in sequences:
                    mapped_ids.add(read_id)
                    seq_info = sequences[read_id]
                    # Create new header with lifted coordinates
                    new_header = f">{seq_info['prefix']}!{chrom}!{start}!{end}!{strand}\n"
                    sys.stdout.write(new_header)
                    sys.stdout.write(f"{seq_info['sequence']}\n")
    except FileNotFoundError:
        print(f"Warning: Lifted BED file not found: {lifted_bed}", file=sys.stderr)
    
    # Process unmapped reads
    try:
        with open(unmapped_bed) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    read_id = fields[3]
                    if read_id in sequences and read_id not in mapped_ids:
                        unmapped_count += 1
                        #seq_info = sequences[read_id]
                        #output.write(f"{seq_info['header']}\n")
                        #output.write(f"{seq_info['sequence']}\n")
    except FileNotFoundError:
        print(f"Warning: Unmapped BED file not found: {unmapped_bed}", file=sys.stderr)
    
    return len(mapped_ids), unmapped_count


def process_fasta(fasta_file, chain_file):
    """Convert all FASTA descriptions from v1.1 to CHM13 and report unmapped percentage."""
    # Create temporary filenames in script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    temp_bed = os.path.join(script_dir, "temp_all.bed")
    lifted_bed = os.path.join(script_dir, "lifted_all.bed")
    unmapped_bed = os.path.join(script_dir, "unmapped_all.bed")
    
    try:
        # Step 1: Create BED file for all reads
        print("Creating BED file from FASTA...", file=sys.stderr)
        sequences, total_descriptions = create_bed_from_fasta(fasta_file, temp_bed)
        
        # Step 2: Run liftOver on entire BED file
        print("Running liftOver on all reads...", file=sys.stderr)
        if not run_liftover(temp_bed, chain_file, lifted_bed, unmapped_bed):
            print("LiftOver failed", file=sys.stderr)
            return
        
        # Step 3: Convert lifted BED back to FASTA
        print("Converting lifted coordinates back to FASTA...", file=sys.stderr)
        mapped_count, unmapped_count = create_fasta_from_lifted_bed(
            lifted_bed, unmapped_bed, sequences
        )
        
        # Report statistics
        print(f"Total descriptions: {total_descriptions}", file=sys.stderr)
        print(f"Mapped descriptions: {mapped_count}", file=sys.stderr)
        print(f"Unmapped descriptions: {unmapped_count}", file=sys.stderr)
        if total_descriptions > 0:
            unmapped_percentage = (unmapped_count / total_descriptions) * 100
            print(f"Percentage unmapped: {unmapped_percentage:.2f}%", file=sys.stderr)
            
    finally:
        # Clean up temporary files
        for f in [temp_bed, lifted_bed, unmapped_bed]:
            try:
                os.remove(f)
            except OSError:
                pass


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert FASTA descriptions from v1.1 to CHM13 using chain files.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-c", "--chain", required=True, help="Chain file for v1.1 to CHM13")
    args = parser.parse_args()

    print(f"Processing FASTA file: {args.fasta}", file=sys.stderr)
    process_fasta(args.fasta, args.chain)


if __name__ == "__main__":
    main()
