
import traceback
from PySide6.QtCore import QObject, QThread, Signal
from dorado_api.pipline import RunPipeline


# ---- Worker that runs in background thread ----
class WorkerThread(QObject):
    log = Signal(str)
    done = Signal(bool, str)

    def __init__(self, trial_name: str, pod5_path: str = "", fastq_path: str = "",
                 bam_path: str = "", nanotel_path: str = "", output_dir: str = "",
                 organism: str = "mouse", do_pod5: bool = False, do_fastq: bool = False,
                 do_nanotel: bool = False, do_align: bool = False, do_r: bool = False,
                 run_filtration: bool = True, run_mapping: bool = True,
                 run_methylation: bool = True):
        super().__init__()
        self.trial_name = trial_name
        self.pod5_path = pod5_path
        self.fastq_path = fastq_path
        self.bam_path = bam_path
        self.nanotel_path = nanotel_path
        self.output_dir = output_dir
        self.organism = organism
        self.do_pod5 = do_pod5
        self.do_fastq = do_fastq
        self.do_nanotel = do_nanotel
        self.do_align = do_align
        self.do_r = do_r
        self.run_filtration = run_filtration
        self.run_mapping = run_mapping
        self.run_methylation = run_methylation
        self._stop_requested = False

    def run(self):
        try:
            self.log.emit("=" * 60)
            self.log.emit("WORKFLOW STARTING")
            self.log.emit("=" * 60)
            self.log.emit(f"Trial: {self.trial_name}")
            self.log.emit(f"Output: {self.output_dir}")
            self.log.emit(f"Organism: {self.organism}")
            self.log.emit("")
            self.log.emit("  Input Paths:")
            self.log.emit(f"  POD5:   {self.pod5_path if self.pod5_path else '(not provided)'}")
            self.log.emit(f"  FASTQ:  {self.fastq_path if self.fastq_path else '(not provided)'}")
            self.log.emit(f"  BAM:    {self.bam_path if self.bam_path else '(not provided)'}")
            self.log.emit(f"  NanoTel: {self.nanotel_path if self.nanotel_path else '(not provided)'}")
            self.log.emit("")
            self.log.emit("Steps to execute:")
            if self.do_pod5:
                self.log.emit("  • POD5 Workflow (basecall → demux → nanotel → align → R-analysis)")
            if self.do_fastq:
                self.log.emit("  • FASTQ Workflow (nanotel → align → R-analysis)")
            if self.do_nanotel:
                self.log.emit("  • NanoTel Analysis")
            if self.do_align:
                self.log.emit("  • Alignment")
            if self.do_r:
                self.log.emit("  • R-Analysis")
                self.log.emit(f"    - Filtration: {self.run_filtration}")
                self.log.emit(f"    - Mapping: {self.run_mapping}")
                self.log.emit(f"    - Methylation: {self.run_methylation}")
            self.log.emit("=" * 60)
            self.log.emit("")

            result, message  = RunPipeline(
                trial_name=self.trial_name,
                pod5_path=self.pod5_path,
                fastq_path=self.fastq_path,
                bam_path=self.bam_path,
                nanotel_path=self.nanotel_path,
                output_dir=self.output_dir,
                organism=self.organism,
                do_pod5=self.do_pod5,
                do_fastq=self.do_fastq,
                do_nanotel=self.do_nanotel,
                do_align=self.do_align,
                do_r=self.do_r,
                run_filtration=self.run_filtration,
                run_mapping=self.run_mapping,
                run_methylation=self.run_methylation,
                log_cb=self.log.emit,
                stop_cb = lambda: self._stop_requested
            )

            if result == 0:

                self.done.emit(True, message)
            else:  # result == 1 (error)
                self.done.emit(False, message)


        except Exception as e:
            self.log.emit(f"Error: {str(e)}")
            self.done.emit(False, f"Workflow failed: {str(e)}")

    def stop(self):
        self._stop_requested = True
        self.log.emit("Stopping worker thread...")

