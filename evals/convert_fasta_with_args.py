import re
import sys
import argparse

#def parse_chain_file(chain_file):
#    """Parse chain file into a dictionary for fast lookup, skipping extra arguments."""
#    chain_dict = {}
#    with open(chain_file) as f:
#        for line in f:
#            parts = line.strip().split()
#            if len(parts) < 12:  # Adjusted to account for the skipped argument
#                continue
#            src_chrom, src_start, src_end, _, _, _, _, dest_chrom, _, dest_start, dest_end, _ = parts[:12]
#            key = (src_chrom, int(src_start), int(src_end))
#            chain_dict[key] = (dest_chrom, int(dest_start), int(dest_end))
#            #print(key, ' -> ', chain_dict[key])
#    return chain_dict

#def liftover_coordinates(chain_dict, chrom, start, end):
#    """Map v1.1 coordinates to CHM13 using the chain dictionary."""
#    for (src_chrom, src_start, src_end), (dest_chrom, dest_start, dest_end) in chain_dict.items():
#        if src_chrom == chrom and src_start <= start < src_end:
#            new_start = dest_start + (start - src_start)
#            new_end = dest_start + (end - src_start)
#            return dest_chrom, new_start, new_end
#    return None, None, None  # Unmapped

def liftover_coordinates(chain_file, chrom, start, end):
    """Map coordinates using liftOver tool."""
    # Create temporary BED file with coordinates
    with open("temp.bed", "w") as f:
        f.write(f"{chrom}\t{start}\t{end}\n")
    
    # Run liftOver command
    try:
        import subprocess
        cmd = ["ext/liftOver", "temp.bed", chain_file, "mapped.bed", "unmapped.bed"]
        subprocess.run(cmd, check=True, capture_output=False)
        
        # Read mapped coordinates
        with open("mapped.bed") as f:
            line = f.readline().strip()
            if line:
                new_chrom, new_start, new_end = line.split()[:3]
                return new_chrom, int(new_start), int(new_end)
                
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Error running liftOver: {e}", file=sys.stderr)
    finally:
        # Clean up temp files
        import os
        for f in ["temp.bed", "mapped.bed", "unmapped.bed"]:
            try:
                os.remove(f)
            except OSError:
                pass
                
    return None, None, None  # Return None if unmapped or error


def process_fasta(fasta_file, chain_file, output_file):
    """Convert all FASTA descriptions from v1.1 to CHM13 and report unmapped percentage."""
    total_descriptions = 0
    unmapped_descriptions = 0

    with open(fasta_file) as fasta, open(output_file, "w") as output:
        for line in fasta:
            if line.startswith(">"):
                total_descriptions += 1
                # Extract chrom, start, end from description
                match = re.match(r">(.*?)!(chr\S+?)_(\S+?)!(\d+?)!(\d+?)!(\S+)", line)
                if match:
                    prefix, chrom, haplotype, start, end, strand = match.groups()
                    start, end = int(start), int(end)
                    # Perform liftover
                    new_chrom, new_start, new_end = liftover_coordinates(chain_file, f"{chrom}_{haplotype}", start, end)
                    if new_chrom:
                        # Write updated description
                        new_desc = f">{prefix}!{new_chrom}!{new_start}!{new_end}!{strand}\n"
                        output.write(new_desc)
                    else:
                        unmapped_descriptions += 1
                        #print(f"Unmapped: {line.strip()}")
                        output.write(line)  # Write original if unmapped
                else:
                    output.write(line)  # Write original if no match
            else:
                output.write(line)  # Write sequence lines unchanged

    # Report unmapped percentage
    unmapped_percentage = (unmapped_descriptions / total_descriptions) * 100 if total_descriptions > 0 else 0
    print(f"Total descriptions: {total_descriptions}")
    print(f"Unmapped descriptions: {unmapped_descriptions}")
    print(f"Percentage unmapped: {unmapped_percentage:.2f}%")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert FASTA descriptions from v1.1 to CHM13 using chain files.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-c", "--chain", required=True, help="Chain file for v1.1 to CHM13")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    # Parse chain file and process FASTA
    print("Parsing chain file...")
    #chain_dict = parse_chain_file(args.chain)
    print(f"Processing FASTA file: {args.fasta}")
    process_fasta(args.fasta, args.chain, args.output)
    print(f"Conversion complete. Output written to {args.output}")

if __name__ == "__main__":
    main()

