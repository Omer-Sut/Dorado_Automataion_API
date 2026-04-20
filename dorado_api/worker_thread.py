import traceback
from PySide6.QtCore import QObject, Signal
from dorado_api.pipeline_runner import run_pipeline


"""Qt worker object that executes the pipeline outside the main UI thread."""
class WorkerThread(QObject):
    """Run the workflow asynchronously and report progress to the GUI.

    Signals:
        log(str): Emits textual log updates for display in the UI.
        done(bool, str): Emits completion status and a final message.
    """

    log = Signal(str)
    done = Signal(bool, str)

    def __init__(
        self,
        trial_name: str,
        pod5_path: str = "",
        fastq_path: str = "",
        bam_path: str = "",
        nanotel_path: str = "",
        output_dir: str = "",
        organism: str = "mouse",
        do_pod5: bool = False,
        do_fastq: bool = False,
        do_nanotel: bool = False,
        do_align: bool = False,
        do_r: bool = False,
        non_pod5_trim_status: str = "auto",
        run_filtration: bool = True,
        run_mapping: bool = True,
        run_methylation: bool = True,
    ):
        """Store workflow settings that will be passed to run_pipeline."""
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
        self.non_pod5_trim_status = non_pod5_trim_status

        self.run_filtration = run_filtration
        self.run_mapping = run_mapping
        self.run_methylation = run_methylation

        self._stop_requested = False

    def stop(self):
        """Request cancellation by setting a stop flag checked by the pipeline."""
        self._stop_requested = True
        self.log.emit("Stop requested by user")

    def run(self):
        """Execute the pipeline and emit final success/failure signals."""
        try:
            self.log.emit("=" * 60)
            self.log.emit("WORKFLOW STARTING")
            self.log.emit("=" * 60)

            # Pass all current settings and callbacks into the pipeline entrypoint.
            status_code, message = run_pipeline(
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
                non_pod5_trim_status=self.non_pod5_trim_status,
                run_filtration=self.run_filtration,
                run_mapping=self.run_mapping,
                run_methylation=self.run_methylation,
                log_cb=self.log.emit,
                stop_cb=lambda: self._stop_requested,
            )

            if self._stop_requested:
                self.done.emit(False, "Cancelled by user")
                return

            if status_code == 0:
                self.log.emit("WORKFLOW FINISHED SUCCESSFULLY")
                self.done.emit(True, message)
            else:
                self.done.emit(False, message)
        
        except Exception as e:
            if "Cancelled by user" in str(e):
                self.done.emit(False, "Cancelled by user")
                return

            err = traceback.format_exc()
            self.log.emit(err)
            self.done.emit(False, err)