#!/usr/bin/env python3
"""
Dorado Post-Processing Automation Script
=======================================

Automates the tedious post-processing steps after Dorado basecalling and demuxing:
1. Batch BAM to FASTQ conversion for all barcodes
2. Automated NanoTel telomere analysis
3. Directory organization and management
4. Batch alignment processing

This script handles the repetitive CLI commands you showed, making your workflow much more efficient.

Author: Omer Sutovsky
Date: 2025
"""

import os
import sys
import json
import logging
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List
from dataclasses import dataclass
from datetime import datetime
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

@dataclass
class LabConfig:
    """Lab-specific configuration settings"""
    references: Dict[str, str] = None
    default_reference: str = "mouse"
    nanotel_script: str = "/home/tzfati/Telomere-Analyzer/NanoTel.R"
    default_kit: str = "SQK-NBD114-24"
    telomere_pattern: str = "CCCTAA"
    min_density: float = 0.5
    fastq_subdir: str = "fastqs"
    nanotel_subdir: str = "nanotel_output"
    aligned_subdir: str = "aligned"

class DoradoPostProcessor:
    """Main class for post-processing automation"""
    
    def __init__(self, base_dir: str, config: LabConfig = None):
        self.logger = None
        self.base_dir = Path(base_dir)
        self.config = config or LabConfig()
        self.setup_logging()
        self.executed_commands = []
        self.command_lock = threading.Lock()
        # Validate tools
        self.validate_tools()
        
    def setup_logging(self):
        """Setup logging"""
        log_file = self.base_dir / f"post_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Post-processor initialized. Base dir: {self.base_dir}")
    
    def validate_tools(self):
        """Check if required tools are available"""
        required_tools = ['samtools', 'dorado', 'Rscript']
        missing_tools = []
        
        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            self.logger.error(f"Missing required tools: {missing_tools}")
            sys.exit(1)
        
        # Check if NanoTel script exists
        if not os.path.exists(self.config.nanotel_script):
            self.logger.warning(f"NanoTel script not found: {self.config.nanotel_script}")

    def run_command(self, command: str, capture_output: bool = False) -> subprocess.CompletedProcess:
        """Run shell command with logging"""
        self.logger.info(f"Running: {command}")

        # Thread-safe command tracking
        with self.command_lock:
            cmd_index = len(self.executed_commands)
            self.executed_commands.append({
                'timestamp': datetime.now().strftime('%H:%M:%S'),
                'command': command,
                'status': 'running'
            })

        try:
            if capture_output:
                result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(command, shell=True, check=True)

            # Update status to success using the specific index
            with self.command_lock:
                self.executed_commands[cmd_index]['status'] = 'success'
            return result

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {command}")
            self.logger.error(f"Error: {e}")

            # Update status to failed using the specific index
            with self.command_lock:
                self.executed_commands[cmd_index]['status'] = 'failed'
                self.executed_commands[cmd_index]['error'] = str(e)

            if capture_output and e.stderr:
                self.logger.error(f"stderr: {e.stderr}")
                with self.command_lock:
                    self.executed_commands[cmd_index]['stderr'] = e.stderr
            raise
    
    def find_demuxed_barcodes(self, demuxed_dir: str) -> List[str]:
        """Find all barcode directories in demuxed output"""
        demux_path = Path(demuxed_dir)
        
        if not demux_path.exists():
            raise FileNotFoundError(f"Demuxed directory not found: {demuxed_dir}")
        
        # Look for barcode directories
        barcode_dirs = []
        for item in demux_path.iterdir():
            if item.is_dir() and item.name.startswith('barcode'):
                barcode_dirs.append(item.name)
        
        barcode_dirs.sort(key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)
        
        self.logger.info(f"Found {len(barcode_dirs)} barcodes: {barcode_dirs}")
        return barcode_dirs
    
    def create_output_structure(self, barcodes: List[str]) -> Dict[str, Path]:
        """Create organized output directory structure"""
        output_dirs = {}
        
        # Create main output directories
        fastq_base = self.base_dir / self.config.fastq_subdir
        nanotel_base = self.base_dir / self.config.nanotel_subdir
        
        fastq_base.mkdir(exist_ok=True)
        nanotel_base.mkdir(exist_ok=True)
        
        # Create barcode-specific directories
        for barcode in barcodes:
            fastq_dir = fastq_base / barcode
            nanotel_dir = nanotel_base / barcode
            
            fastq_dir.mkdir(exist_ok=True)
            nanotel_dir.mkdir(exist_ok=True)
            
            output_dirs[barcode] = {
                'fastq': fastq_dir,
                'nanotel': nanotel_dir
            }
        
        self.logger.info(f"Created output structure for {len(barcodes)} barcodes")
        return output_dirs
    
    def convert_bam_to_fastq(self, demuxed_dir: str, output_dirs: Dict, parallel: bool = True) -> None:
        """Convert all BAM files to FASTQ for each barcode"""
        self.logger.info("Starting BAM to FASTQ conversion...")
        
        demux_path = Path(demuxed_dir)
        conversion_tasks = []
        
        # Collect all conversion tasks
        for barcode, dirs in output_dirs.items():
            barcode_dir = demux_path / barcode
            if not barcode_dir.exists():
                self.logger.warning(f"Barcode directory not found: {barcode_dir}")
                continue
            
            # Find all BAM files for this barcode
            bam_files = list(barcode_dir.glob("*.bam"))
            if not bam_files:
                self.logger.warning(f"No BAM files found in {barcode_dir}")
                continue
            
            # Create conversion tasks
            for i, bam_file in enumerate(bam_files, 1):
                output_fastq = dirs['fastq'] / f"{barcode}{chr(96+i)}.fastq"  # barcode2a.fastq, barcode2b.fastq, etc.
                
                conversion_tasks.append({
                    'barcode': barcode,
                    'bam_file': bam_file,
                    'output_fastq': output_fastq
                })
        
        self.logger.info(f"Found {len(conversion_tasks)} BAM files to convert")
        
        if parallel:
            self._convert_parallel(conversion_tasks)
        else:
            self._convert_sequential(conversion_tasks)
        
        self.logger.info("BAM to FASTQ conversion completed")
    
    def _convert_parallel(self, tasks: List[Dict]) -> None:
        """Convert BAM files in parallel"""
        def convert_single(task):
            command = f"samtools fastq {task['bam_file']} > {task['output_fastq']}"
            try:
                self.run_command(command)
                return f"✓ {task['barcode']}: {task['bam_file'].name} -> {task['output_fastq'].name}"
            except Exception as e:
                return f"✗ {task['barcode']}: Failed to convert {task['bam_file'].name} - {str(e)}"
        
        max_workers = min(4, len(tasks))  # Limit to 4 parallel processes
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(convert_single, task): task for task in tasks}
            
            for future in as_completed(future_to_task):
                result = future.result()
                self.logger.info(result)
    
    def _convert_sequential(self, tasks: List[Dict]) -> None:
        """Convert BAM files sequentially"""
        for task in tasks:
            command = f"samtools fastq {task['bam_file']} > {task['output_fastq']}"
            try:
                self.run_command(command)
                self.logger.info(f"✓ {task['barcode']}: {task['bam_file'].name} -> {task['output_fastq'].name}")
            except Exception as e:
                self.logger.error(f"✗ {task['barcode']}: Failed to convert {task['bam_file'].name} - {str(e)}")
    
    def run_nanotel_analysis(self, output_dirs: Dict, parallel: bool = True) -> None:
        """Run NanoTel telomere analysis on all barcodes"""
        if not os.path.exists(self.config.nanotel_script):
            self.logger.warning("NanoTel script not found, skipping telomere analysis")
            return
        
        self.logger.info("Starting NanoTel telomere analysis...")
        
        nanotel_tasks = []
        for barcode, dirs in output_dirs.items():
            fastq_dir = dirs['fastq']
            nanotel_output = dirs['nanotel']
            
            # Check if FASTQ files exist
            fastq_files = list(fastq_dir.glob("*.fastq"))
            if not fastq_files:
                self.logger.warning(f"No FASTQ files found for {barcode}, skipping NanoTel")
                continue
            
            nanotel_tasks.append({
                'barcode': barcode,
                'input_dir': fastq_dir,
                'output_dir': nanotel_output
            })
        
        if parallel:
            self._run_nanotel_parallel(nanotel_tasks)
        else:
            self._run_nanotel_sequential(nanotel_tasks)
        
        self.logger.info("NanoTel analysis completed")
    
    def _run_nanotel_parallel(self, tasks: List[Dict]) -> None:
        """Run NanoTel analysis in parallel"""
        def run_single_nanotel(task):
            command = (
                f"Rscript --vanilla {self.config.nanotel_script} "
                f"-i {task['input_dir']} "
                f"--save_path {task['output_dir']} "
                f"--patterns {self.config.telomere_pattern} "
                f"--min_density {self.config.min_density} "
                f"--use_filter"
            )
            
            try:
                self.run_command(command)
                return f"✓ NanoTel completed for {task['barcode']}"
            except Exception as e:
                return f"✗ NanoTel failed for {task['barcode']}: {str(e)}"
        
        max_workers = min(2, len(tasks))  # Limit R processes
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(run_single_nanotel, task): task for task in tasks}
            
            for future in as_completed(future_to_task):
                result = future.result()
                self.logger.info(result)
    
    def _run_nanotel_sequential(self, tasks: List[Dict]) -> None:
        """Run NanoTel analysis sequentially"""
        for task in tasks:
            command = (
                f"Rscript --vanilla {self.config.nanotel_script} "
                f"-i {task['input_dir']} "
                f"--save_path {task['output_dir']} "
                f"--patterns {self.config.telomere_pattern} "
                f"--min_density {self.config.min_density} "
                f"--use_filter"
            )
            
            try:
                self.run_command(command)
                self.logger.info(f"✓ NanoTel completed for {task['barcode']}")
            except Exception as e:
                self.logger.error(f"✗ NanoTel failed for {task['barcode']}: {str(e)}")
    
    def run_alignment(self, demuxed_dir: str, reference_choice: str = "mouse") -> None:
        """Run alignment on demuxed data"""
        self.logger.info("Starting alignment...")
        
        aligned_dir = self.base_dir / self.config.aligned_subdir
        aligned_dir.mkdir(exist_ok=True)

        if reference_choice == "human":
            reference_path = self.config.references["human"]
        else:
            reference_path = self.config.references["mouse"]
        command = (
            f"dorado aligner -r "
            f"--output-dir {aligned_dir} "
            f"--emit-summary "
            f"{reference_path} "
            f"{demuxed_dir}"
        )
        
        try:
            self.run_command(command)
            self.logger.info("✓ Alignment completed successfully")
        except Exception as e:
            self.logger.error(f"✗ Alignment failed: {str(e)}")
    
    def generate_summary_report(self, output_dirs: Dict) -> str:
        """Generate a summary report of the processing results"""
        report_lines = [
            "=" * 60,
            "DORADO POST-PROCESSING SUMMARY REPORT",
            "=" * 60,
            f"Processing date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Base directory: {self.base_dir}",
            "",
            "PROCESSED BARCODES:",
        ]
        
        for barcode, dirs in output_dirs.items():
            fastq_files = list(dirs['fastq'].glob("*.fastq"))
            nanotel_files = list(dirs['nanotel'].glob("*"))
            
            report_lines.extend([
                f"  {barcode}:",
                f"    FASTQ files: {len(fastq_files)}",
                f"    NanoTel outputs: {len(nanotel_files)}",
                f"    FASTQ directory: {dirs['fastq']}",
                f"    NanoTel directory: {dirs['nanotel']}",
                ""
            ])
        
        # Check alignment results
        aligned_dir = self.base_dir / self.config.aligned_subdir
        if aligned_dir.exists():
            aligned_files = list(aligned_dir.glob("*"))
            report_lines.extend([
                "ALIGNMENT:",
                f"  Output directory: {aligned_dir}",
                f"  Files generated: {len(aligned_files)}",
                ""
            ])

            # ADD COMMAND HISTORY SECTION
            if self.executed_commands:
                report_lines.extend([
                    "",
                    "EXECUTED COMMANDS:",
                    "=" * 50,
                ])

                for i, cmd_info in enumerate(self.executed_commands, 1):
                    status_icon = "✓" if cmd_info['status'] == 'success' else "✗"
                    report_lines.extend([
                        f"  {i}. [{cmd_info['timestamp']}] {status_icon} {cmd_info['status'].upper()}",
                        f"     Command: {cmd_info['command']}",
                    ])

                    if cmd_info['status'] == 'failed':
                        report_lines.append(f"     Error: {cmd_info.get('error', 'Unknown error')}")
                        if 'stderr' in cmd_info:
                            report_lines.append(f"     Stderr: {cmd_info['stderr']}")

                    report_lines.append("")  # Empty line between commands

        report_lines.extend([
            "NEXT STEPS:",
            "  1. Check FASTQ files for quality and completeness",
            "  2. Review NanoTel telomere analysis results",
            "  3. Examine alignment statistics and coverage",
            "  4. Proceed with downstream analysis",
            "=" * 60
        ])
        
        # Save report to file
        report_file = self.base_dir / f"processing_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        with open(report_file, 'w') as f:
            f.write('\n'.join(report_lines))

        # ALSO SAVE COMMANDS TO SEPARATE FILE
        self.save_command_history()

        self.logger.info(f"Summary report saved: {report_file}")
        return '\n'.join(report_lines)

    def save_command_history(self) -> str:
        """Save detailed command history to a separate file"""
        commands_file = self.base_dir / f"commands_executed_{datetime.now().strftime('%Y%m%d_%H%M%S')}.sh"

        script_lines = [
            "#!/bin/bash",
            f"# Commands executed by Dorado Post-Processor",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"# Base directory: {self.base_dir}",
            "",
            "# NOTE: This is a record of commands that were executed.",
            "# You can review, modify, and re-run these commands manually if needed.",
            "",
        ]

        for i, cmd_info in enumerate(self.executed_commands, 1):
            status_comment = "# SUCCESS" if cmd_info['status'] == 'success' else "# FAILED"
            script_lines.extend([
                f"# Command {i} - {cmd_info['timestamp']} - {status_comment}",
                cmd_info['command'],
                ""
            ])

        with open(commands_file, 'w') as f:
            f.write('\n'.join(script_lines))

        # Make it executable
        os.chmod(commands_file, 0o755)

        self.logger.info(f"Command history saved: {commands_file}")
        return str(commands_file)

    def run_complete_workflow(self, demuxed_dir: str, reference_choice: str = "mouse",
                              skip_nanotel: bool = False, skip_alignment: bool = False, parallel: bool = True) -> None:
        """Run the complete post-processing workflow"""
        self.logger.info("Starting complete post-processing workflow...")
        
        try:
            # Step 1: Find barcodes
            barcodes = self.find_demuxed_barcodes(demuxed_dir)
            if not barcodes:
                self.logger.error("No barcodes found in demuxed directory")
                return
            
            # Step 2: Create output structure
            output_dirs = self.create_output_structure(barcodes)
            
            # Step 3: Convert BAM to FASTQ
            self.convert_bam_to_fastq(demuxed_dir, output_dirs, parallel)
            
            # Step 4: Run NanoTel analysis (optional)
            if not skip_nanotel:
                self.run_nanotel_analysis(output_dirs, parallel)
            else:
                self.logger.info("Skipping NanoTel analysis")
            
            # Step 5: Run alignment (optional)
            if not skip_alignment:
                self.run_alignment(demuxed_dir, reference_choice)
            else:
                self.logger.info("Skipping alignment")
            
            # Step 6: Generate summary report
            report = self.generate_summary_report(output_dirs)
            print("\n" + report)
            
            self.logger.info("Complete workflow finished successfully!")
            
        except Exception as e:
            self.logger.error(f"Workflow failed: {str(e)}")
            raise


def create_lab_config_file(config_path: str = "lab_config.json") -> None:
    """Create a template lab configuration file"""
    default_config = {
        "references":{
            "human": "/home/tzfati/Documents/Human/T2T-CHM13-150KchrHeadTail-Yincluded-telo-trimmed.fasta",
            "mouse": "/home/tzfati/Documents/Mouse-assembly/mouse.241018.v1.1.0.-40k_edges-telo-trimmed-Y-included.fasta"},
        "default_reference": "mouse",
        "nanotel_script": "/home/tzfati/Telomere-Analyzer/NanoTel.R",
        "default_kit": "SQK-NBD114-24",
        "telomere_pattern": "CCCTAA",
        "min_density": 0.5,
        "fastq_subdir": "fastqs",
        "nanotel_subdir": "nanotel_output",
        "aligned_subdir": "aligned",

    }
    
    with open(config_path, 'w') as f:
        json.dump(default_config, f, indent=2)
    
    print(f"Lab configuration template created: {config_path}")
    print("Edit this file to customize paths for your lab setup.")
    print("\nTo use human reference:")
    print("1. Edit the 'human' reference path")
    print("2. Change 'default_reference' to 'human'")

def main():
    """Main function for command line interface"""
    parser = argparse.ArgumentParser(
        description="Dorado Post-Processing Automation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete workflow on demuxed data
  python dorado_post_processor.py /path/to/Trial_74_Aki /path/to/demuxed/
  
  # Skip NanoTel analysis
  python dorado_post_processor.py /path/to/Trial_74_Aki /path/to/demuxed/ --skip-nanotel
  
  # Only convert BAM to FASTQ
  python dorado_post_processor.py /path/to/Trial_74_Aki /path/to/demuxed/ --skip-nanotel --skip-alignment
  
  # Use custom config
  python dorado_post_processor.py /path/to/Trial_74_Aki /path/to/demuxed/ --config lab_config.json
  
  # Create config template
  python dorado_post_processor.py --create-config
        """
    )
    
    parser.add_argument('base_dir', nargs='?', help='Base output directory (e.g., /path/to/Trial_74_Aki)')
    parser.add_argument('demuxed_dir', nargs='?', help='Demuxed directory containing barcodes')
    
    parser.add_argument('--config', type=str, help='Lab configuration JSON file')
    parser.add_argument('--skip-nanotel', action='store_true', help='Skip NanoTel telomere analysis')
    parser.add_argument('--skip-alignment', action='store_true', help='Skip alignment step')
    parser.add_argument('--no-parallel', action='store_true', help='Disable parallel processing')
    parser.add_argument('--create-config', action='store_true', help='Create lab configuration template')
    
    # Individual step options
    parser.add_argument('--only-fastq', action='store_true', help='Only convert BAM to FASTQ')
    parser.add_argument('--only-nanotel', action='store_true', help='Only run NanoTel analysis')
    parser.add_argument('--only-alignment', action='store_true', help='Only run alignment')
    parser.add_argument('--reference', type=str, choices=['mouse', 'human'], default='mouse',
                        help='Reference genome to use for alignment (default: mouse)')
    
    args = parser.parse_args()
    
    # Handle config creation
    if args.create_config:
        create_lab_config_file()
        return
    
    # Validate required arguments
    if not args.base_dir or not args.demuxed_dir:
        parser.error("Both base_dir and demuxed_dir are required (unless using --create-config)")
    
    # Load configuration
    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            config_data = json.load(f)
        config = LabConfig(**config_data)
    else:
        config = LabConfig()
    
    # Initialize processor
    processor = DoradoPostProcessor(args.base_dir, config)
    
    try:
        if args.only_fastq:
            # Only BAM to FASTQ conversion
            barcodes = processor.find_demuxed_barcodes(args.demuxed_dir)
            output_dirs = processor.create_output_structure(barcodes)
            processor.convert_bam_to_fastq(args.demuxed_dir, output_dirs, not args.no_parallel)
            
        elif args.only_nanotel:
            # Only NanoTel analysis
            barcodes = processor.find_demuxed_barcodes(args.demuxed_dir)
            output_dirs = processor.create_output_structure(barcodes)
            processor.run_nanotel_analysis(output_dirs, not args.no_parallel)
            
        elif args.only_alignment:
            # Only alignment
            processor.run_alignment(args.demuxed_dir, args.reference)
            
        else:
            # Complete workflow
            processor.run_complete_workflow(
                args.demuxed_dir,
                reference_choice=args.reference,
                skip_nanotel=args.skip_nanotel,
                skip_alignment=args.skip_alignment,
                parallel=not args.no_parallel
            )
        
        print("\n✓ Post-processing completed successfully!")
        
    except Exception as e:
        print(f"\n✗ Post-processing failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
