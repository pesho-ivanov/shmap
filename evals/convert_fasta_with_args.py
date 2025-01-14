import re
import sys
import argparse
import subprocess
import os


def create_bed_from_fasta(fasta_file, bed_file):
    """Create a BED file containing coordinates for all reads in the FASTA file."""
    total_descriptions = 0
    sequences = {}  # Store sequences for later reconstruction
    
    with open(fasta_file) as fasta, open(bed_file, "w") as bed:
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
                    bed.write(f"{chrom}_{haplotype}\t{start}\t{end}\t{current_id}\n")
            elif current_id is not None:
                sequences[current_id]['sequence'] += line.strip()
    
    return sequences, total_descriptions


def run_liftover(input_bed, chain_file, output_bed, unmapped_bed):
    """Run liftOver on the entire BED file at once."""
    try:
        cmd = ["ext/liftOver", input_bed, chain_file, output_bed, unmapped_bed]
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running liftOver: {e}", file=sys.stderr)
        return False


def create_fasta_from_lifted_bed(lifted_bed, unmapped_bed, sequences, output_file):
    """Convert lifted BED coordinates back to FASTA format."""
    mapped_ids = set()
    unmapped_count = 0
    
    with open(output_file, "w") as output:
        # Process successfully lifted coordinates
        try:
            with open(lifted_bed) as f:
                for line in f:
                    chrom, start, end, read_id = line.strip().split()[:4]
                    if read_id in sequences:
                        mapped_ids.add(read_id)
                        seq_info = sequences[read_id]
                        # Create new header with lifted coordinates
                        new_header = f">{seq_info['prefix']}!{chrom}!{start}!{end}!{seq_info['strand']}\n"
                        output.write(new_header)
                        output.write(f"{seq_info['sequence']}\n")
        except FileNotFoundError:
            pass  # Handle case where no reads were mapped
        
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
            pass  # Handle case where no reads were unmapped
    
    return len(mapped_ids), unmapped_count


def process_fasta(fasta_file, chain_file, output_file):
    """Convert all FASTA descriptions from v1.1 to CHM13 and report unmapped percentage."""
    # Create temporary filenames
    temp_bed = "temp_all.bed"
    lifted_bed = "lifted_all.bed"
    unmapped_bed = "unmapped_all.bed"
    
    try:
        # Step 1: Create BED file for all reads
        print("Creating BED file from FASTA...")
        sequences, total_descriptions = create_bed_from_fasta(fasta_file, temp_bed)
        
        # Step 2: Run liftOver on entire BED file
        print("Running liftOver on all reads...")
        if not run_liftover(temp_bed, chain_file, lifted_bed, unmapped_bed):
            print("LiftOver failed", file=sys.stderr)
            return
        
        # Step 3: Convert lifted BED back to FASTA
        print("Converting lifted coordinates back to FASTA...")
        mapped_count, unmapped_count = create_fasta_from_lifted_bed(
            lifted_bed, unmapped_bed, sequences, output_file
        )
        
        # Report statistics
        print(f"Total descriptions: {total_descriptions}")
        print(f"Mapped descriptions: {mapped_count}")
        print(f"Unmapped descriptions: {unmapped_count}")
        if total_descriptions > 0:
            unmapped_percentage = (unmapped_count / total_descriptions) * 100
            print(f"Percentage unmapped: {unmapped_percentage:.2f}%")
            
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
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    print(f"Processing FASTA file: {args.fasta}")
    process_fasta(args.fasta, args.chain, args.output)
    print(f"Conversion complete. Output written to {args.output}")


if __name__ == "__main__":
    main()
