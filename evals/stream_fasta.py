#!/usr/bin/env python3
import sys
import subprocess
import os

def read_fasta_batch(fh, batch_size=1000):
    """Read FASTA sequences in batches."""
    batch = []
    sequence = []
    header = None
    
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            if header is not None:
                # Store the previous sequence
                batch.append((header, ''.join(sequence)))
                sequence = []
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
            header = line
        else:
            sequence.append(line)
            
    # Don't forget the last sequence
    if header is not None and sequence:
        batch.append((header, ''.join(sequence)))
    
    # Don't forget the last batch
    if batch:
        yield batch

def write_fasta_batch(batch):
    """Convert batch to FASTA format string."""
    return '\n'.join(f"{header}\n{seq}" for header, seq in batch)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.fa chain_file", file=sys.stderr)
        sys.exit(1)
        
    input_file = sys.argv[1]
    chain_file = sys.argv[2]
    
    # Get absolute paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    lift_fasta = os.path.join(script_dir, "convert_fasta_with_args.py")
    chain_file = os.path.abspath(chain_file)
    
    with open(input_file) as fh:
        for batch in read_fasta_batch(fh, 1000):
            # Create the lift_fasta command
            cmd = ["python", lift_fasta, "-f", "-", "-c", chain_file]
            
            # Create the subprocess
            proc = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=script_dir  # Run in evals directory for relative paths
            )
            
            # Write batch to stdin and get output
            batch_fasta = write_fasta_batch(batch) + "\n"  # Add newline
            stdout, stderr = proc.communicate(batch_fasta)
            
            # Write lifted sequences to stdout
            sys.stdout.write(stdout)
            
            # Write any errors to stderr
            if stderr:
                sys.stderr.write(stderr)
            
            if proc.returncode != 0:
                sys.stderr.write(f"Error processing batch\n")
                sys.exit(1)

if __name__ == "__main__":
    main()
