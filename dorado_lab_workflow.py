#!/usr/bin/env python3
"""
Dorado Lab Workflow Manager
==========================

Handles the complete Dorado workflow from basecalling to post-processing,
specifically designed for your lab's common use cases.

This script automates:
1. Initial basecalling with your typical parameters
2. Demultiplexing 
3. Post-processing (BAM->FASTQ, NanoTel, alignment)

Based on your actual workflow patterns.
"""

import os
import re
import sys
import json
import logging
import argparse
import subprocess
from pathlib import Path
from typing import Optional
from dataclasses import dataclass, asdict
from datetime import datetime
from dorado_post_processing import DoradoPostProcessor, LabConfig


@dataclass
class LabWorkflowConfig:
    """Configuration for typical lab workflows"""
    # Your common settings
    default_model: str = "/home/tzfati/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
    default_kit: str = "SQK-NBD114-24"
    mouse_reference: str = "/home/tzfati/Documents/Mouse-assembly/mouse.241018.v1.1.0.-40k_edges-telo-trimmed-Y-included.fasta"
    human_reference: str = "/home/tzfati/Documents/Human/T2T-CHM13-150KchrHeadTail-Yincluded-telo-trimmed.fasta"
    min_qscore: int = 9
    modified_bases: str = "5mCG_5hmCG"

    # Common paths
    nanotel_script: str = "/home/tzfati/Telomere-Analyzer/NanoTel.R"
    telomere_pattern: str = "CCCTAA"
    min_density: float = 0.5

    # Processing options
    no_trim: bool = True
    recursive: bool = True
    sort_bam: bool = True
    emit_summary: bool = True

    # R Analysis Configuration
    run_r_analysis: bool = True
    run_nanotel_analysis: bool = True
    run_mapping_analysis: bool = True
    run_methylation_analysis: bool = True
    stop_on_error: bool = True

    # R Analysis Parameters
    nanotel_density_threshold: float = 0.75
    nanotel_max_telomere_start: int = 150
    mapping_min_mapq: int = 10
    mapping_head_max_start: int = 5000
    mapping_tail_min_end: int = 35000
    methylation_head_max_pos: int = 5000
    methylation_tail_min_pos: int = 145000
    methylation_create_plots: bool = True
    methylation_create_shiny_app: bool = False


class DoradoLabWorkflow:
    """Main workflow manager for lab-specific Dorado processing"""
    
    def __init__(self, trial_name: str, base_output_dir: str, config: LabWorkflowConfig = None):
        self.logger = None
        self.trial_name = trial_name
        self.base_output_dir = Path(base_output_dir)
        self.config = config or LabWorkflowConfig()
        
        # Create trial-specific directory structure
        self.trial_dir = self.base_output_dir / trial_name
        self.setup_directories()
        self.setup_logging()
    
    def setup_directories(self):
        """Create the standard directory structure for a trial"""
        self.dirs = {
            'base': self.trial_dir,
            'rebasecalled': self.trial_dir / 'rebasecalled',
            'demuxed': self.trial_dir / 'demuxed',
            'fastqs': self.trial_dir / 'fastqs',
            'nanotel_output': self.trial_dir / 'nanotel_output',
            'aligned': self.trial_dir / 'aligned',
            'logs': self.trial_dir / 'logs'
        }
        
        # Create all directories
        for dir_path in self.dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def setup_logging(self):
        """Setup logging for the workflow"""
        log_file = self.dirs['logs'] / f"{self.trial_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Lab workflow initialized for trial: {self.trial_name}")
    
    def build_basecall_command(self, pod5_input: str) -> str:
        """Build the basecalling command with your typical parameters"""
        timestamp = datetime.now().strftime('%Y-%m-%d_T%H-%M-%S')
        output_file = self.dirs['rebasecalled'] / f"calls_{timestamp}.bam"
        
        cmd_parts = [
            "dorado", "basecaller",
            f"--min-qscore {int(self.config.min_qscore)}",
            "-r" if self.config.recursive else "",
            f"--modified-bases {self.config.modified_bases}",
            "--no-trim" if self.config.no_trim else "",
            f"--kit-name {self.config.default_kit}",
            f"--reference {self.config.mouse_reference}",
            f"--output-dir {self.dirs['rebasecalled']}",
            self.config.default_model,
            f'"{pod5_input}"'
        ]
        
        # Remove empty strings
        cmd_parts = [part for part in cmd_parts if part]
        command = " ".join(cmd_parts)
        
        self.logger.info(f"Basecall command: {command}")
        return command, output_file
    
    def build_demux_command(self, basecalled_bam: str) -> str:
        """Build the demultiplexing command"""
        cmd_parts = [
            "dorado", "demux",
            f"--output-dir {self.dirs['demuxed']}",
            # "--no-classify",
            f"--kit-name {self.config.default_kit} "
            "--no-trim" if self.config.no_trim else "",
            "--sort-bam" if self.config.sort_bam else "",
            "--emit-summary",
            f'"{basecalled_bam}"'
        ]
        
        # Remove empty strings
        cmd_parts = [part for part in cmd_parts if part]
        command = " ".join(cmd_parts)
        
        self.logger.info(f"Demux command: {command}")
        return command

    def organize_demuxed_files(self) -> bool:
        """Organize demuxed BAM files into barcode subdirectories"""
        self.logger.info("Organizing demuxed files into barcode directories...")

        demuxed_dir = self.dirs['demuxed']

        # Find all BAM files in the demuxed directory
        bam_files = list(demuxed_dir.glob("*.bam"))
        bai_files = list(demuxed_dir.glob("*.bam.bai"))

        if not bam_files:
            self.logger.warning("No BAM files found to organize")
            return True

        # Group files by barcode or unclassified
        barcode_files = {}

        for file_path in bam_files + bai_files:
            # Check if it's an unclassified file
            if 'unclassified' in file_path.name:
                folder_name = 'unclassified'
            else:
                # Extract barcode from filename (e.g., "barcode01" from "something_barcode01.bam")
                barcode_match = re.search(r'barcode\d+', file_path.name, re.IGNORECASE)
                if barcode_match:
                    folder_name = barcode_match.group().lower()
                else:
                    self.logger.warning(f"Could not determine barcode for file: {file_path.name}")
                    continue

            if folder_name not in barcode_files:
                barcode_files[folder_name] = []
            barcode_files[folder_name].append(file_path)

        # Create directories and move files
        for folder_name, files in barcode_files.items():
            target_dir = demuxed_dir / folder_name
            target_dir.mkdir(exist_ok=True)

            for file_path in files:
                new_path = target_dir / file_path.name
                file_path.rename(new_path)
                self.logger.info(f"Moved {file_path.name} to {folder_name}/")

        self.logger.info(f"Organized {len(barcode_files)} folders including barcodes and unclassified")
        return True

    def run_command(self, command: str) -> bool:
        """Run a command and handle errors"""
        self.logger.info(f"Executing: {command}")
        
        try:
            result = subprocess.run(command, shell=True, check=True)
            self.logger.info("✓ Command completed successfully")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"✗ Command failed with exit code {e.returncode}")
            return False
    
    def run_basecalling(self, pod5_input: str) -> Optional[str]:
        """Run the basecalling step"""
        self.logger.info("=== STARTING BASECALLING ===")
        
        if not os.path.exists(pod5_input):
            self.logger.error(f"POD5 input not found: {pod5_input}")
            return None
        
        command, output_file = self.build_basecall_command(pod5_input)
        
        if self.run_command(command):
            # Find the actual output file (dorado might create it with a different timestamp)
            bam_files = list(self.dirs['rebasecalled'].glob("*.bam"))
            if bam_files:
                # Get the most recent BAM file
                latest_bam = max(bam_files, key=lambda x: x.stat().st_mtime)
                self.logger.info(f"✓ Basecalling completed. Output: {latest_bam}")
                return str(latest_bam)
            else:
                self.logger.error("No BAM output file found after basecalling")
                return None
        else:
            self.logger.error("Basecalling failed")
            return None

    def run_demuxing(self, basecalled_bam: str) -> bool:
        """Run the demultiplexing step"""
        self.logger.info("=== STARTING DEMULTIPLEXING ===")

        if not os.path.exists(basecalled_bam):
            self.logger.error(f"Basecalled BAM not found: {basecalled_bam}")
            return False

        command = self.build_demux_command(basecalled_bam)

        if self.run_command(command):
            # Organize files into barcode subdirectories
            if self.organize_demuxed_files():
                # Check if demuxed files were created
                demuxed_files = list(self.dirs['demuxed'].glob("**/*.bam"))
                self.logger.info(f"✓ Demultiplexing completed. Created {len(demuxed_files)} demuxed files")
                return True
            else:
                self.logger.error("Failed to organize demuxed files")
                return False
        else:
            self.logger.error("Demultiplexing failed")
            return False

    def run_post_processing(self, skip_nanotel: bool = False, reference_choice: str = "mouse") -> bool:
        """Run the post-processing steps using the post-processor"""
        self.logger.info("=== STARTING POST-PROCESSING ===")

        try:
            # Create proper LabConfig for post-processor
            post_config = LabConfig(
                references={
                    "mouse": self.config.mouse_reference,
                    "human": "/home/tzfati/Documents/Human/T2T-CHM13-150KchrHeadTail-Yincluded-telo-trimmed.fasta"
                    # Add human reference
                },
                nanotel_script=self.config.nanotel_script,
                default_kit=self.config.default_kit,
                telomere_pattern=self.config.telomere_pattern,
                min_density=self.config.min_density,
                fastq_subdir="fastqs",
                nanotel_subdir="nanotel_output",
                aligned_subdir="aligned"
            )

            # Initialize post-processor
            processor = DoradoPostProcessor(str(self.trial_dir), post_config)

            # Run the complete post-processing workflow
            processor.run_complete_workflow(
                str(self.dirs['demuxed']),
                reference_choice=reference_choice,
                skip_nanotel=skip_nanotel,
                skip_alignment=False,
                parallel=True
            )

            self.logger.info("✓ Post-processing completed successfully")
            return True

        except Exception as e:
            self.logger.error(f"Post-processing failed: {str(e)}")
            return False

    def run_r_analysis_pipeline(self) -> bool:
        """Run the complete R analysis pipeline"""
        self.logger.info("=== STARTING R ANALYSIS PIPELINE ===")

        try:
            # Generate R pipeline configuration
            r_config = self.generate_r_pipeline_config()
            r_config_file = self.dirs['logs'] / "r_pipeline_config.json"

            with open(r_config_file, 'w') as f:
                json.dump(r_config, f, indent=2)

            # Change to r_analysis directory and run R pipeline
            cmd = f"cd r_analysis && Rscript main_analysis_pipeline.R {r_config_file.absolute()}"

            if self.run_command(cmd):
                self.logger.info("R analysis pipeline completed successfully")
                return True
            else:
                self.logger.error("R analysis pipeline failed")
                return False

        except Exception as e:
            self.logger.error(f"R analysis pipeline failed: {str(e)}")
            return False

    def generate_r_pipeline_config(self) -> dict:
        """Generate R pipeline configuration from workflow config"""

        # Ensure R analysis directories exist
        (self.trial_dir / "nanotel_output").mkdir(exist_ok=True)
        (self.trial_dir / "mapping_output").mkdir(exist_ok=True)
        (self.trial_dir / "methylation_output").mkdir(exist_ok=True)

        return {
            "base_output_dir": str(self.trial_dir),
            "run_nanotel_analysis": self.config.run_nanotel_analysis,
            "run_mapping_analysis": self.config.run_mapping_analysis,
            "run_methylation_analysis": self.config.run_methylation_analysis,
            "stop_on_error": self.config.stop_on_error,

            "nanotel_analysis": {
                "input_dir": str(self.dirs['nanotel_output']),
                "output_dir": str(self.trial_dir / "nanotel_output"),
                "density_threshold": self.config.nanotel_density_threshold,
                "max_telomere_start": self.config.nanotel_max_telomere_start
            },

            "mapping_analysis": {
                "alignment_summary_path": str(self.dirs['aligned'] / "alignment_summary.txt"),
                "filtered_nanotel_dir": str(self.trial_dir / "nanotel_output"),
                "bam_dir": str(self.dirs['demuxed']),
                "output_dir": str(self.trial_dir / "mapping_output"),
                "min_mapq": self.config.mapping_min_mapq,
                "head_max_start": self.config.mapping_head_max_start,
                "tail_min_end": self.config.mapping_tail_min_end
            },

            "methylation_analysis": {
                "pileup_bed_dir": str(self.trial_dir / "mapping_output"),
                "output_dir": str(self.trial_dir / "methylation_output"),
                "head_max_pos": self.config.methylation_head_max_pos,
                "tail_min_pos": self.config.methylation_tail_min_pos,
                "create_plots": self.config.methylation_create_plots,
                "create_shiny_app": self.config.methylation_create_shiny_app
            }
        }

    def generate_command_script(self, pod5_input: str, script_name: str = None) -> str:
        """Generate a shell script with all the commands for manual execution"""
        if not script_name:
            script_name = f"{self.trial_name}_commands.sh"
        
        script_path = self.dirs['logs'] / script_name
        
        # Build all commands
        basecall_cmd, _ = self.build_basecall_command(pod5_input)
        
        # For script generation, we need to use a placeholder for the BAM file
        placeholder_bam = f"{self.dirs['rebasecalled']}/calls_TIMESTAMP.bam"
        demux_cmd = self.build_demux_command(placeholder_bam)

        script_content = f"""#!/bin/bash
# Generated Dorado workflow commands for {self.trial_name}
# Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

set -e  # Exit on any error

echo "Starting Dorado workflow for {self.trial_name}"

# Step 1: Basecalling
echo "=== BASECALLING ==="
{basecall_cmd}

# Step 2: Demultiplexing (update BAM path as needed)
echo "=== DEMULTIPLEXING ==="
# Note: Replace 'calls_TIMESTAMP.bam' with the actual output filename from basecalling
{demux_cmd.replace(placeholder_bam, '${self.dirs["rebasecalled"]}/calls_ACTUAL_TIMESTAMP.bam')}

# Step 3: Post-processing (BAM to FASTQ conversion)
echo "=== POST-PROCESSING ==="
python3 dorado_post_processor.py "{self.trial_dir}" "{self.dirs['demuxed']}"

echo "Workflow completed for {self.trial_name}"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make script executable
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"Command script generated: {script_path}")
        return str(script_path)

    def run_complete_workflow(self, pod5_input: str, reference_choice: str = "mouse",
                              skip_nanotel: bool = False, skip_r_analysis: bool = False) -> bool:
        """Run the complete workflow from basecalling to post-processing"""
        self.logger.info(f"=== STARTING COMPLETE WORKFLOW FOR {self.trial_name} ===")

        # Step 1: Basecalling
        existing_bams = list(self.dirs['rebasecalled'].glob("*.bam"))
        if existing_bams:
            latest_bam = max(existing_bams, key=lambda x: x.stat().st_mtime)
            self.logger.info(f"Using existing basecalled BAM: {latest_bam}")
            basecalled_bam = str(latest_bam)
        else:
            basecalled_bam = self.run_basecalling(pod5_input)
            if not basecalled_bam:
                return False

        # Step 2: Demultiplexing
        if not self.run_demuxing(basecalled_bam):
            return False

        # Step 3: Post-processing (pass reference choice)
        if not self.run_post_processing(skip_nanotel, reference_choice):
            return False

        # Step 4: R Analysis Pipeline (NEW)
        if not skip_r_analysis and not self.run_r_analysis_pipeline():
            return False

        self.logger.info(f"=== COMPLETE WORKFLOW FINISHED FOR {self.trial_name} ===")
        self.generate_workflow_summary()
        return True
    
    def generate_workflow_summary(self) -> None:
        """Generate a summary of the workflow results"""
        summary_lines = [
            "=" * 80,
            f"WORKFLOW SUMMARY FOR {self.trial_name}",
            "=" * 80,
            f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "DIRECTORY STRUCTURE:",
            f"  Base directory: {self.trial_dir}",
            f"  Rebasecalled: {self.dirs['rebasecalled']}",
            f"  Demuxed: {self.dirs['demuxed']}",
            f"  FASTQ files: {self.dirs['fastqs']}",
            f"  NanoTel output: {self.dirs['nanotel_output']}",
            f"  Aligned: {self.dirs['aligned']}",
            f"  Logs: {self.dirs['logs']}",
            "",
            "FILE COUNTS:",
        ]
        
        # Count files in each directory
        for name, path in self.dirs.items():
            if name != 'base':
                file_count = len(list(path.rglob("*"))) if path.exists() else 0
                summary_lines.append(f"  {name}: {file_count} files")
        
        # Find barcodes
        if self.dirs['demuxed'].exists():
            barcodes = [d.name for d in self.dirs['demuxed'].iterdir() if d.is_dir() and d.name.startswith('barcode')]
            summary_lines.extend([
                "",
                f"BARCODES PROCESSED: {len(barcodes)}",
                f"  {', '.join(sorted(barcodes))}"
            ])
        
        summary_lines.extend([
            "",
            "CONFIGURATION USED:",
            f"  Model: {self.config.default_model}",
            f"  Kit: {self.config.default_kit}",
            f"  Min Q-score: {self.config.min_qscore}",
            f"  Modified bases: {self.config.modified_bases}",
            f"  Reference: {os.path.basename(self.config.mouse_reference)}",
            "",
            "=" * 80
        ])
        
        # Save summary
        summary_file = self.dirs['logs'] / f"{self.trial_name}_summary.txt"
        with open(summary_file, 'w') as f:
            f.write('\n'.join(summary_lines))
        
        # Print to console
        print('\n'.join(summary_lines))
        self.logger.info(f"Workflow summary saved: {summary_file}")


def create_workflow_config(config_path: str = "workflow_config.json") -> None:
    """Create a workflow configuration template"""
    config = LabWorkflowConfig()

    # Convert to dictionary and update with human reference
    config_dict = asdict(config)

    with open(config_path, 'w') as f:
        json.dump(config_dict, f, indent=2)

    print(f"Workflow configuration template created: {config_path}")
    print("Edit this file to customize your lab's default settings.")
    print("Both mouse and human references are included.")

def main():
    """Main function for command line interface"""
    parser = argparse.ArgumentParser(
        description="Dorado Lab Workflow Manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete workflow for a new trial
  python dorado_lab_workflow.py Trial_75_New /media/tzfati/Storage/Nanopore_runs/T_75_New/combined_POD5
  
  # Only basecalling and demuxing
  python dorado_lab_workflow.py Trial_75_New /path/to/pod5 --only-demux
  
  # Only post-processing (if basecalling/demux already done)
  python dorado_lab_workflow.py Trial_75_New --only-post-process
  
  # Generate command script without running
  python dorado_lab_workflow.py Trial_75_New /path/to/pod5 --generate-script-only
  
  # Use custom configuration
  python dorado_lab_workflow.py Trial_75_New /path/to/pod5 --config workflow_config.json
  
  # Create configuration template
  python dorado_lab_workflow.py --create-config
        """
    )
    
    parser.add_argument('trial_name', nargs='?', help='Trial name (e.g., Trial_74_Aki)')
    parser.add_argument('pod5_input', nargs='?', help='POD5 input directory or file')
    
    parser.add_argument('--output-dir', type=str, default='/home/tzfati/Desktop/minknow_runs',
                       help='Base output directory (default: /home/tzfati/Desktop/minknow_runs)')
    parser.add_argument('--config', type=str, help='Workflow configuration JSON file')
    
    # Workflow control options
    parser.add_argument('--only-basecall', action='store_true', help='Only run basecalling')
    parser.add_argument('--only-demux', action='store_true', help='Run basecalling and demuxing only')
    parser.add_argument('--only-post-process', action='store_true', help='Only run post-processing')
    parser.add_argument('--only-align', action='store_true', help='Only run alignment')
    parser.add_argument('--skip-nanotel', action='store_true', help='Skip NanoTel analysis')
    parser.add_argument('--generate-script-only', action='store_true', help='Generate shell script without running')
    parser.add_argument('--reference', type=str, choices=['mouse', 'human'], default='mouse',
                        help='Reference genome to use for alignment and NanoTel (default: mouse)')

    # Configuration
    parser.add_argument('--create-config', action='store_true', help='Create workflow configuration template')

    # Rscript arguments
    parser.add_argument('--skip-r-analysis', action='store_true', help='Skip R analysis pipeline')
    parser.add_argument('--only-r-analysis', action='store_true', help='Only run R analysis')
    parser.add_argument('--r-config', type=str, help='Additional R analysis configuration')

    args = parser.parse_args()
    
    # Handle config creation
    if args.create_config:
        create_workflow_config()
        return
    
    # Validate required arguments
    if not args.trial_name:
        parser.error("trial_name is required (unless using --create-config)")

    if not args.only_post_process and not args.only_align and not args.only_r_analysis and not args.pod5_input:
        parser.error("pod5_input is required unless using --only-post-process, --only-align, or --only-r-analysis")
    
    # Load configuration
    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            config_data = json.load(f)
        config = LabWorkflowConfig(**config_data)
    else:
        config = LabWorkflowConfig()
    
    # Initialize workflow
    workflow = DoradoLabWorkflow(args.trial_name, args.output_dir, config)
    
    try:
        if args.generate_script_only:
            # Generate script without running
            if not args.pod5_input:
                parser.error("pod5_input required for script generation")
            script_path = workflow.generate_command_script(args.pod5_input)
            print(f"\n✓ Command script generated: {script_path}")
            print("Review and execute the script manually when ready.")
            
        elif args.only_basecall:
            # Only basecalling
            result = workflow.run_basecalling(args.pod5_input)
            if result:
                print(f"\n✓ Basecalling completed. Output: {result}")
            else:
                print("\n✗ Basecalling failed")
                sys.exit(1)
                
        elif args.only_demux:
            # Basecalling + demuxing
            basecalled_bam = workflow.run_basecalling(args.pod5_input)
            if basecalled_bam and workflow.run_demuxing(basecalled_bam):
                print(f"\n✓ Basecalling and demuxing completed")
            else:
                print("\n✗ Basecalling or demuxing failed")
                sys.exit(1)
                
        elif args.only_post_process:
            # Only post-processing
            if workflow.run_post_processing(skip_nanotel=args.skip_nanotel):
                print("\n✓ Post-processing completed")
            else:
                print("\n✗ Post-processing failed")
                sys.exit(1)

        elif args.only_r_analysis:
            # Only R analysis
            if workflow.run_r_analysis_pipeline():
                print("\nR analysis completed")
            else:
                print("\nR analysis failed")
                sys.exit(1)

        else:
            # Complete workflow
            if workflow.run_complete_workflow(args.pod5_input, reference_choice=args.reference,
                                              skip_nanotel=args.skip_nanotel):
                print(f"\n✓ Complete workflow finished for {args.trial_name}!")
            else:
                print(f"\n✗ Workflow failed for {args.trial_name}")
                sys.exit(1)
        
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Workflow failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()