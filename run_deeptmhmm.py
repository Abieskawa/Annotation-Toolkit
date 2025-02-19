import biolib
import sys,stat,os,time,subprocess

def format_time(elapsed_time):
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    return f"{hours} hr {minutes} mins {seconds} seconds"

def check_nvidia_gpu():
    try:
        # Run the nvidia-smi command to check for NVIDIA GPU
        subprocess.run(["nvidia-smi"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print("NVIDIA GPU detected.")
        if os.getenv("BIOLIB_DOCKER_RUNTIME") != "nvidia":
            print("Environment variable BIOLIB_DOCKER_RUNTIME is not set to 'nvidia'. Setting it now...")
            os.environ["BIOLIB_DOCKER_RUNTIME"] = "nvidia"
            print("Environment variable BIOLIB_DOCKER_RUNTIME set to 'nvidia'.")
    except subprocess.CalledProcessError:
        print("Warning: NVIDIA GPU not detected. The program will run on CPU. For optimal performance, ensure an NVIDIA GPU is available.")

def make_files_read_only(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            try:
                # Set file to read-only
                os.chmod(file_path, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)
                print(f"Set read-only permissions for: {file_path}")
            except Exception as e:
                print(f"Error setting permissions for {file_path}: {e}")

def run_deeptmhmm(input_file, output_dir, basename):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Record start time
    start_time = time.time()

    # Enable real-time streaming of progress
    biolib.utils.STREAM_STDOUT = True

    # Load DeepTMHMM model
    deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24')

    # Prepare CLI arguments
    args = f"--fasta {input_file}"

    # Run the job locally
    deeptmhmm_job = deeptmhmm.cli(args=args, machine='local', machine_count=4)

    # Save result files
    deeptmhmm_job.save_files(output_dir)

    try:
        # Adjust the paths if necessary (if "predicted_topologies.3line" is inside output_dir)
        parser_input_3line = os.path.join(output_dir, "predicted_topologies.3line")
        parser_output_tsv  = os.path.join(output_dir, f"{basename}_deeptmhmm.tsv")
        parser_tmr_gff3  = os.path.join(output_dir, "TMRs.gff3")

        subprocess.run([
            "DeepTMHMM_parser.py", 
            parser_input_3line, 
            parser_output_tsv, 
            "-topology", 
            parser_tmr_gff3,
            "-delimiter", "1"
        ], check=True)

        print("DeepTMHMM_parser.py ran successfully.")
    except FileNotFoundError:
        print("Error: DeepTMHMM_parser.py not found. Make sure it's in your PATH or provide the full path.")
    except subprocess.CalledProcessError as e:
        print(f"Error running DeepTMHMM_parser.py: {e}")

    make_files_read_only(output_dir)

    # Record end time
    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"DeepTMHMM results saved to {output_dir}")
    print(f"Time taken: {format_time(elapsed_time)}")

if __name__ == "__main__":
    check_nvidia_gpu()

    if len(sys.argv) != 4:
        print("Usage: python run_deeptmhmm.py <input_file> <output_dir> <output_basename>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    output_basename = sys.argv[3]
    run_deeptmhmm(input_file, output_dir, output_basename)
