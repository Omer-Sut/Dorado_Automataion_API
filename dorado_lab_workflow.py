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
            "--no-classify",
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
        #TODO: need to solve bam files not in barcode directories!!
        self.logger.info("=== STARTING DEMULTIPLEXING ===")
        
        if not os.path.exists(basecalled_bam):
            self.logger.error(f"Basecalled BAM not found: {basecalled_bam}")
            return False
        
        command = self.build_demux_command(basecalled_bam)
        
        if self.run_command(command):
            # Check if demuxed files were created
            demuxed_files = list(self.dirs['demuxed'].glob("**/*.bam"))
            self.logger.info(f"✓ Demultiplexing completed. Created {len(demuxed_files)} demuxed files")
            return True
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
                              skip_nanotel: bool = False) -> bool:
        """Run the complete workflow from basecalling to post-processing"""
        self.logger.info(f"=== STARTING COMPLETE WORKFLOW FOR {self.trial_name} ===")

        # Step 1: Basecalling
        basecalled_bam = self.run_basecalling(pod5_input)
        if not basecalled_bam:
            return False

        # Step 2: Demultiplexing
        if not self.run_demuxing(basecalled_bam):
            return False

        # Step 3: Post-processing (pass reference choice)
        if not self.run_post_processing(skip_nanotel, reference_choice):
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
    
    args = parser.parse_args()
    
    # Handle config creation
    if args.create_config:
        create_workflow_config()
        return
    
    # Validate required arguments
    if not args.trial_name:
        parser.error("trial_name is required (unless using --create-config)")
    
    if not args.only_post_process and not args.only_align and not args.pod5_input:
        parser.error("pod5_input is required unless using --only-post-process or --only-align")
    
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