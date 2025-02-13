import sys
import re

def rename_proteins(input_file, output_file, new_prefix):
    """Rename protein sequences in a FASTA file by replacing the 'g' prefix and updating the suffix."""
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    header = line.strip()
                    # Update to match '>g<number>.t<number>' and replace appropriately.
                    header = re.sub(r'^>g(\d+)\.t(\d+)', fr'>{new_prefix}_g\1.p\2', header)
                    outfile.write(header + '\n')
                else:
                    outfile.write(line)
        print(f"Renaming complete. Updated file saved as {output_file}.")
    except FileNotFoundError:
        print(f"Error: File {input_file} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_proteins.py <input_file> <output_file> <new_prefix>")
        print("<new_prefix>_gx.tx")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    new_prefix = sys.argv[3]

    rename_proteins(input_file, output_file, new_prefix)
