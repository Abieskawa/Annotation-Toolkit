#!/usr/bin/env python3
import argparse
import os

def remove_trailing_star(input_file, output_file, log_file, command):
    """Removes trailing '*' from amino acid sequences in the input file."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Write header lines unchanged
                outfile.write(line)
            else:
                # Remove trailing '*' and write the sequence
                outfile.write(line.rstrip().rstrip('*') + '\n')

    # Write the command to the log file
    with open(log_file, 'a') as logfile:
        logfile.write(f"Command used: {command}\n")
        logfile.write(f"Processed input file: {input_file} -> output file: {output_file}\n")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Remove trailing '*' from amino acid sequences.")
    parser.add_argument('-I', '--input', required=True, help="Input .aa file")
    parser.add_argument('-O', '--output', required=True, help="Output .aa file")
    parser.add_argument('-L', '--log', required=True, help="Log file to store the command")

    args = parser.parse_args()

    # Construct the command used for debugging and logging
    command_used = f"python3 remove_star.py -I {os.path.basename(args.input)} -O {os.path.basename(args.output)} -L {os.path.basename(args.log)}"
    print(f"Command used: {command_used}")

    # Remove trailing '*'
    remove_trailing_star(args.input, args.output, args.log, command_used)
