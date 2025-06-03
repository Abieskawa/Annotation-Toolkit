#!/usr/bin/env python3
"""
PacBio RNA-seq Processing Pipeline
Handles cutadapt -> picard -> isoseq -> lordec workflow with conda environment management
Supports multiple samples with corresponding short reads for polishing
"""

import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path
import time

class PacBioPipeline:
    def __init__(self, base_dir=".", threads=128, five_prime_adapters=None, three_prime_adapters=None, anywhere_adapters=None):
        self.base_dir = Path(base_dir)
        self.threads = threads
        self.five_prime_adapters = five_prime_adapters or []
        self.three_prime_adapters = three_prime_adapters or []
        self.anywhere_adapters = anywhere_adapters or []
        self.setup_logging()
        self.setup_directories()
        
    def setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('pacbio_pipeline.log'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def setup_directories(self):
        """Create necessary output directories"""
        dirs = [
            "05_preprocessed_RNA_PacBio"
        ]
        for d in dirs:
            (self.base_dir / d).mkdir(exist_ok=True)
            
    def run_conda_command(self, command, env_name, description):
        """Run command in specified conda environment"""
        self.logger.info(f"Starting: {description}")
        self.logger.info(f"Environment: {env_name}")
        self.logger.info(f"Command: {command}")
        
        # Create the full conda command
        full_command = f"conda run -n {env_name} {command}"
        
        try:
            start_time = time.time()
            result = subprocess.run(
                full_command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                cwd=self.base_dir
            )
            
            duration = time.time() - start_time
            self.logger.info(f"✓ Completed: {description} ({duration:.1f}s)")
            
            # Log stdout and stderr to the main log file
            if result.stdout:
                self.logger.info(f"STDOUT from {description}:")
                for line in result.stdout.strip().split('\n'):
                    if line.strip():  # Only log non-empty lines
                        self.logger.info(f"  {line}")
                        
            if result.stderr:
                self.logger.info(f"STDERR from {description}:")
                for line in result.stderr.strip().split('\n'):
                    if line.strip():  # Only log non-empty lines
                        self.logger.info(f"  {line}")
                
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"✗ Failed: {description}")
            self.logger.error(f"Exit code: {e.returncode}")
            
            if e.stdout:
                self.logger.error(f"STDOUT from failed {description}:")
                for line in e.stdout.strip().split('\n'):
                    if line.strip():
                        self.logger.error(f"  {line}")
                        
            if e.stderr:
                self.logger.error(f"STDERR from failed {description}:")
                for line in e.stderr.strip().split('\n'):
                    if line.strip():
                        self.logger.error(f"  {line}")
            return False
            
    def step1_cutadapt(self, input_fastq, sample_id):
        """Step 1: Cutadapt primer removal and poly-A trimming"""
        if not self.five_prime_adapters and not self.three_prime_adapters and not self.anywhere_adapters:
            self.logger.error(f"No adapters specified for {sample_id}. Use -g, -a, or --anywhere options or skip cutadapt step.")
            return False
            
        flnc_fastq = f"05_preprocessed_RNA_PacBio/{sample_id}_flnc.fastq"
        final_fastq = f"05_preprocessed_RNA_PacBio/{sample_id}_no_internal_primers.fastq"
        
        # Build first cutadapt command
        cmd1_parts = ["cutadapt"]
        
        # Add 5' adapters
        for adapter in self.five_prime_adapters:
            cmd1_parts.append(f'-g "{adapter}"')
            
        # Add 3' adapters  
        for adapter in self.three_prime_adapters:
            cmd1_parts.append(f'-a {adapter}')
            
        cmd1_parts.extend([
            "--poly-a", "--trimmed-only",
            f"--cores {self.threads}",
            f"-o {flnc_fastq}",
            input_fastq
        ])
        
        # Build second cutadapt command for anywhere adapters (only if specified)
        if self.anywhere_adapters:
            cmd2_parts = ["cutadapt"]
            
            for adapter in self.anywhere_adapters:
                cmd2_parts.append(f"--anywhere {adapter}")
                
            cmd2_parts.extend([
                "--minimum-length 5", "--discard-trimmed",
                f"--cores {self.threads}",
                f"-o {final_fastq}",
                flnc_fastq
            ])
            
            # Combine both commands
            command = f'bash -c \'{" ".join(cmd1_parts)} && {" ".join(cmd2_parts)}\''
        else:
            # No anywhere adapters - just rename the output
            command = f'bash -c \'{" ".join(cmd1_parts)} && mv {flnc_fastq} {final_fastq}\''
        
        return self.run_conda_command(command, "annotation", f"Cutadapt primer removal - {sample_id}")
        
    def step2_picard(self, sample_id):
        """Step 2: Convert FASTQ to unaligned BAM using Picard"""
        input_fastq = f"05_preprocessed_RNA_PacBio/{sample_id}_no_internal_primers.fastq"
        output_bam = f"05_preprocessed_RNA_PacBio/{sample_id}_no_internal_primers.bam"
        
        command = (
            f'picard FastqToSam '
            f'-F1 {input_fastq} '
            f'-O {output_bam} '
            f'-SM {sample_id} '
            f'-PL PACBIO'
        )
        
        return self.run_conda_command(command, "picard", f"Picard FastqToSam conversion - {sample_id}")
        
    def step3_isoseq(self, sample_id):
        """Step 3: IsoSeq clustering"""
        input_bam = f"05_preprocessed_RNA_PacBio/{sample_id}_no_internal_primers.bam"
        output_bam = f"05_preprocessed_RNA_PacBio/{sample_id}_clustered.bam"
        
        command = (
            f'isoseq cluster {input_bam} {output_bam} '
            f'--verbose --use-qvs -j {self.threads}'
        )
        
        return self.run_conda_command(command, "py27", f"IsoSeq clustering - {sample_id}")
        
    
    def step4_lordec(self, sample_id, short_reads):
        """Step 4: LoRDEC error correction"""
        # IsoSeq generates compressed FASTA files
        input_fasta_gz = f"05_preprocessed_RNA_PacBio/{sample_id}_clustered.hq.fasta.gz"
        input_fasta = f"05_preprocessed_RNA_PacBio/{sample_id}_clustered.hq.fasta"
        output_fasta = f"05_preprocessed_RNA_PacBio/{sample_id}_polished.hq.fasta"
        
        # Check if compressed file exists
        if not self.check_file_exists(input_fasta_gz, f"IsoSeq HQ FASTA (compressed) - {sample_id}"):
            return False
        
        # Decompress the FASTA file in place (removes .gz file)
        self.logger.info(f"Decompressing {input_fasta_gz} in place")
        decompress_cmd = f"gunzip {input_fasta_gz}"
        if not self.run_conda_command(decompress_cmd, "annotation", f"Decompress HQ FASTA - {sample_id}"):
            return False
        
        if not short_reads:
            self.logger.warning(f"No short reads specified for {sample_id} - skipping LoRDEC polishing")
            # Copy the decompressed file
            command = f"cp {input_fasta} {output_fasta}"
            return self.run_conda_command(command, "annotation", f"Copy unpolished FASTA - {sample_id}")
        
        # Remove duplicates from short reads list
        unique_reads = list(set(short_reads))
        reads_str = ",".join(unique_reads)
        
        command = (
            f'lordec-correct '
            f'-i {input_fasta} '
            f'-2 {reads_str} '
            f'-o {output_fasta} '
            f'-T {self.threads} '
            f'-k 19 -s 3'
        )
        
        return self.run_conda_command(command, "annotation", f"LoRDEC error correction - {sample_id}")
        
    def check_file_exists(self, filepath, description):
        """Check if required input file exists"""
        if not os.path.exists(filepath):
            self.logger.error(f"Required file not found: {filepath} ({description})")
            return False
        return True
        
    def run_pipeline(self, sample_id=None, input_fastq=None, short_reads=None, skip_steps=None):
        """Run the complete pipeline for a single sample"""
        skip_steps = skip_steps or []
        short_reads = short_reads or []
        
        self.logger.info(f"Starting PacBio pipeline for sample: {sample_id}")
        self.logger.info(f"Input FASTQ: {input_fastq}")
        self.logger.info(f"Short reads: {short_reads}")
        
        # Check input files
        if not self.check_file_exists(input_fastq, f"PacBio long reads - {sample_id}"):
            return False
            
        for read_file in short_reads:
            if not self.check_file_exists(read_file, f"Illumina short reads - {sample_id}"):
                return False
        
        # Run pipeline steps
        steps = [
            (1, "cutadapt", lambda: self.step1_cutadapt(input_fastq, sample_id)),
            (2, "picard", lambda: self.step2_picard(sample_id)),
            (3, "isoseq", lambda: self.step3_isoseq(sample_id)),
            (4, "lordec", lambda: self.step4_lordec(sample_id, short_reads))
        ]
        
        for step_num, step_name, step_func in steps:
            if step_num in skip_steps:
                self.logger.info(f"Skipping step {step_num}: {step_name} for {sample_id}")
                continue
                
            self.logger.info(f"=" * 60)
            self.logger.info(f"STEP {step_num}: {step_name.upper()} - {sample_id}")
            self.logger.info(f"=" * 60)
            
            if not step_func():
                self.logger.error(f"Pipeline failed at step {step_num}: {step_name} for {sample_id}")
                return False
                
        self.logger.info(f"Sample {sample_id} completed successfully!")
        return True


def parse_input_files(input_args):
    """Parse input files with colon syntax: input.fastq:short1.fastq,short2.fastq"""
    samples = {}
    
    for i, input_spec in enumerate(input_args):
        if ':' in input_spec:
            # Format: input.fastq:short1.fastq,short2.fastq
            input_fastq, short_reads_str = input_spec.split(':', 1)
            short_reads = [f.strip() for f in short_reads_str.split(',') if f.strip()]
        else:
            # Just input file, no short reads
            input_fastq = input_spec
            short_reads = []
        
        # Generate sample ID from filename
        sample_id = Path(input_fastq).stem
        
        samples[sample_id] = {
            'input_fastq': input_fastq,
            'short_reads': short_reads
        }
    
    return samples


def main():
    parser = argparse.ArgumentParser(
        description="PacBio RNA-seq processing pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single sample with polishing reads
  python pacbio_pipeline.py \\
    -i "A.fastq:B_1.fastq,B_2.fastq" \\
    -g "AAGCAGTGGTATCAACGCAGAGTACT{30}" \\
    -a AAGCAGTGGTATCAACGCAGAGTAC
    
  # Multiple samples
  python pacbio_pipeline.py \\
    -i "sample1.fastq:short1_1.fastq,short1_2.fastq" \\
    -i "sample2.fastq:short2_1.fastq,short2_2.fastq" \\
    -g "ADAPTER{30}" \\
    -a ADAPTER
    
  # Single sample without polishing reads
  python pacbio_pipeline.py \\
    -i "A.fastq" \\
    -g "ADAPTER{30}" \\
    -a ADAPTER
    
  # Multiple adapters
  python pacbio_pipeline.py \\
    -i "A.fastq:B_1.fastq,B_2.fastq" \\
    -g "ADAPTER1{30}" \\
    -g "ADAPTER2{25}" \\
    -a ADAPTER3 \\
    -a ADAPTER4
    
  # With internal primer removal
  python pacbio_pipeline.py \\
    -i "A.fastq:B_1.fastq,B_2.fastq" \\
    -g "ADAPTER{30}" \\
    -a ADAPTER \\
    --anywhere INTERNAL_PRIMER
        """
    )
    
    # Input files (required)
    parser.add_argument('-i', '--input-fastq', action='append', default=[], required=True,
                       help='Input PacBio FASTQ file with optional short reads: "input.fastq:short1.fastq,short2.fastq" (use multiple -i for multiple samples)')
    
    # General arguments
    parser.add_argument('-t', '--threads', type=int, default=128,
                        help='Number of threads (default: 128)')
    parser.add_argument('-d', '--base-dir', default='.',
                        help='Base directory for output (default: current directory)')
    parser.add_argument('-g', '--five-prime-adapter', action='append', default=[],
                        help='5\' adapter sequence for cutadapt -g (can be used multiple times)')
    parser.add_argument('-a', '--three-prime-adapter', action='append', default=[],
                        help='3\' adapter sequence for cutadapt -a (can be used multiple times)')
    parser.add_argument('--anywhere', action='append', default=[],
                        help='Adapter sequence for cutadapt --anywhere (can be used multiple times). Optional - only runs if specified.')
    parser.add_argument('--skip', type=int, nargs='+',
                        help='Skip specific steps (1=cutadapt, 2=picard, 3=isoseq, 4=lordec)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show commands that would be run without executing')
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = PacBioPipeline(
        base_dir=args.base_dir, 
        threads=args.threads,
        five_prime_adapters=args.five_prime_adapter,
        three_prime_adapters=args.three_prime_adapter,
        anywhere_adapters=args.anywhere
    )
    
    if args.dry_run:
        pipeline.logger.info("DRY RUN MODE - Commands will be displayed but not executed")
        # You could add dry-run logic here
        return
    
    # Parse input files and run pipeline
    samples = parse_input_files(args.input_fastq)
    
    success_count = 0
    total_count = len(samples)
    
    for sample_id, sample_data in samples.items():
        pipeline.logger.info(f"\n{'='*80}")
        pipeline.logger.info(f"PROCESSING SAMPLE {sample_id} ({success_count + 1}/{total_count})")
        pipeline.logger.info(f"{'='*80}")
        
        if pipeline.run_pipeline(
            sample_id=sample_id,
            input_fastq=sample_data['input_fastq'],
            short_reads=sample_data['short_reads'],
            skip_steps=args.skip
        ):
            success_count += 1
        else:
            pipeline.logger.error(f"Failed to process sample {sample_id}")
    
    pipeline.logger.info(f"\n{'='*80}")
    pipeline.logger.info(f"PIPELINE SUMMARY: {success_count}/{total_count} samples completed successfully")
    pipeline.logger.info(f"{'='*80}")
    
    success = success_count == total_count
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()