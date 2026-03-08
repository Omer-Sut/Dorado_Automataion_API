import sys
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt, QThread
from PySide6.QtGui import QFont, QColor, QIcon
from PySide6.QtCore import QSize
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QFileDialog,
    QComboBox, QTextEdit, QCheckBox, QMessageBox,
    QFrame, QGraphicsDropShadowEffect, QToolButton
)

from dorado_api.stream_to_gui import StreamToGui
from dorado_api.worker_thread import WorkerThread


# =============================================================================
# THEME
# =============================================================================

class Theme:
    BG = "#F3F6F9"
    CARD = "#FFFFFF"
    PRIMARY = "#2F855A"
    BORDER = "#E5E7EB"
    TEXT = "#111827"
    MUTED = "#6B7280"
    LOG_BG = "#0F172A"

    @staticmethod
    def stylesheet():
        return f"""
        QWidget {{
            background-color: {Theme.BG};
            color: {Theme.TEXT};
            font-family: 'Segoe UI';
        }}

        QLineEdit, QComboBox {{
            background-color: white;
            border: 1px solid {Theme.BORDER};
            border-radius: 8px;
            padding: 10px 12px;
            min-height: 38px;
            font-size: 11pt;
        }}

        QPushButton {{
            background-color: {Theme.PRIMARY};
            color: white;
            border-radius: 8px;
            padding: 8px 16px;
            font-size: 11pt;
            font-weight: bold;
        }}

        QPushButton:hover {{
            background-color: #276749;
        }}

        QPushButton:disabled {{
            background-color: #D1D5DB;
            color: #6B7280;
        }}

        QTextEdit {{
            background-color: {Theme.LOG_BG};
            color: #10B981;
            border-radius: 8px;
            font-family: Consolas;
            font-size: 10pt;
        }}
        """


# =============================================================================
# HELPERS
# =============================================================================

def shadow(widget):
    """ Apply a drop shadow effect to a widget"""
    s = QGraphicsDropShadowEffect()
    s.setBlurRadius(24)
    s.setYOffset(5)
    s.setColor(QColor(0, 0, 0, 50))
    widget.setGraphicsEffect(s)


class Card(QFrame):
    """ A card-like widget with a title and a body"""
    def __init__(self, title):
        super().__init__()
        self.setStyleSheet(f"background:{Theme.CARD}; border-radius:12px;")
        shadow(self)

        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.setContentsMargins(14, 14, 14, 14)

        lbl = QLabel(title)
        lbl.setFont(QFont("Segoe UI", 12, QFont.Bold))
        layout.addWidget(lbl)

        self.body = QVBoxLayout()
        self.body.setSpacing(10)
        layout.addLayout(self.body)


# =============================================================================
# VALIDATION
# =============================================================================

def validate_advanced_selection(
        *,
        do_pod5,
        do_fastq,
        do_nanotel,
        do_align,
        do_r,
        do_methylation,
        pod5_path,
        fastq_path,
        bam_path,
        nanotel_path,
):
    """ Validate the advanced workflow configuration. Returns a list of human-readable error messages. """
    errors = []

    has_pod5 = bool(pod5_path)
    has_fastq = bool(fastq_path)
    has_bam = bool(bam_path)
    has_nanotel = bool(nanotel_path)

    if do_pod5 and not has_pod5:
        errors.append("POD5 workflow selected, but no POD5 directory was provided.")

    if do_fastq and not has_fastq:
        errors.append("FASTQ workflow selected, but no FASTQ directory was provided.")

    if do_nanotel and not (has_fastq or has_nanotel):
        errors.append(
            "NanoTel selected, but no FASTQ input or existing NanoTel summaries were provided."
        )

    if do_align and not (has_fastq or has_bam):
        errors.append(
            "Alignment selected, but no FASTQ or BAM files were provided."
        )

    if do_r and not (has_nanotel or has_bam):
        errors.append(
            "R analysis selected, but no NanoTel summaries or BAM files were provided."
        )

    if do_methylation and not has_bam:
        errors.append(
            "Methylation analysis selected, but no BAM file was provided."
        )

    return errors


def validate_input_directories(pod5, fastq, bam, nanotel):
    """
    Validate that provided directories actually contain the expected files.
    Returns a list of human-readable error messages.
    """
    errors = []

    def has_files(path, patterns):
        if not path:
            return False
        p = Path(path)
        if not p.exists():
            return False
        for pat in patterns:
            if list(p.rglob(pat)):
                return True
        return False

    if pod5 and not has_files(pod5, ["*.pod5"]):
        errors.append("POD5 directory does not contain any .pod5 files.")

    if fastq and not has_files(fastq, ["*.fastq", "*.fastq.gz"]):
        errors.append("FASTQ directory does not contain any FASTQ files (.fastq / .fastq.gz).")

    if bam and not has_files(bam, ["*.bam"]):
        errors.append("BAM directory does not contain any .bam files.")

    if nanotel and not has_files(nanotel, ["*summary*.csv"]):
        errors.append("NanoTel directory does not contain summary CSV files.")

    return errors


# =============================================================================
# MAIN WINDOW
# =============================================================================

class AppWindow(QWidget):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Telomere Analyzer")
        self.resize(1000, 760)
        self.setStyleSheet(Theme.stylesheet())

        # Thread management
        self.thread = None
        self.worker = None

        main = QVBoxLayout(self)
        main.setSpacing(14)
        main.setContentsMargins(16, 16, 16, 16)

        # ================= HEADER =================
        header = QFrame()
        header.setStyleSheet("""
        QFrame {
            background: qlineargradient(
                x1:0, y1:0, x2:1, y2:0,
                stop:0 #2F855A,
                stop:1 #2563EB
            );
            border-radius: 16px;
        }
        """)
        shadow(header)

        hl = QVBoxLayout(header)
        hl.setContentsMargins(18, 18, 18, 18)

        title = QLabel("Telomere Analyzer")
        title.setStyleSheet("color:white; font-size:22pt; font-weight:bold")

        hl.addWidget(title)
        main.addWidget(header)

        # ================= INPUT PATHS =================
        inputs = Card("Input Paths")

        def path_row(label):
            r = QHBoxLayout()
            lbl = QLabel(label)
            lbl.setFixedWidth(80)
            lbl.setFont(QFont("Segoe UI", 11, QFont.Bold))

            edit = QLineEdit()
            edit.setPlaceholderText(f"Select {label} directory…")
            edit.setMinimumWidth(520)

            btn = QToolButton()
            btn.setIcon(QIcon(QIcon.fromTheme("folder")))
            btn.setIconSize(QSize(28, 28))
            btn.setFixedSize(44, 44)
            btn.setToolTip(f"Browse for {label} directory")
            btn.clicked.connect(lambda: self.browse(edit))

            r.addWidget(lbl)
            r.addWidget(edit, 1)
            r.addWidget(btn)
            inputs.body.addLayout(r)
            return edit

        self.pod5 = path_row("POD5")
        self.fastq = path_row("FASTQ")
        self.bam = path_row("BAM")
        self.nanotel = path_row("NanoTel")

        main.addWidget(inputs)

        # ================= OUTPUT PATH =================
        output_card = Card("Output Directory")

        out_row = QHBoxLayout()
        out_label = QLabel("Output")
        out_label.setFixedWidth(80)
        out_label.setFont(QFont("Segoe UI", 11, QFont.Bold))

        self.output_dir = QLineEdit()
        self.output_dir.setPlaceholderText("Select output directory…")

        out_btn = QToolButton()
        out_btn.setIcon(QIcon(QIcon.fromTheme("folder")))
        out_btn.setIconSize(QSize(28, 28))
        out_btn.setFixedSize(44, 44)
        out_btn.setToolTip("Browse for output directory")
        out_btn.clicked.connect(lambda: self.browse(self.output_dir))

        out_row.addWidget(out_label)
        out_row.addWidget(self.output_dir, 1)
        out_row.addWidget(out_btn)

        output_card.body.addLayout(out_row)
        main.addWidget(output_card)

        # ================= ORGANISM SELECTION =================
        organism_card = Card("Configuration")

        org_row = QHBoxLayout()
        org_label = QLabel("Organism")
        org_label.setFixedWidth(80)
        org_label.setFont(QFont("Segoe UI", 11, QFont.Bold))

        self.organism_combo = QComboBox()
        self.organism_combo.addItems(["mouse", "human", "fish"])

        org_row.addWidget(org_label)
        org_row.addWidget(self.organism_combo, 1)
        organism_card.body.addLayout(org_row)
        main.addWidget(organism_card)

        # ================= QUICK SETUP =================
        quick = Card("Quick Setup")

        self.quick_mode = QComboBox()
        self.quick_mode.addItems([
            "Complete POD5 Workflow",
            "FASTQ Workflow",
            "NanoTel Analysis",
            "R Analysis"
        ])

        hint = QLabel("Choose a high-level mode. All internal steps are auto-configured.")
        hint.setStyleSheet("color:#6B7280; font-size:9pt")

        quick.body.addWidget(QLabel("Workflow mode:"))
        quick.body.addWidget(self.quick_mode)
        quick.body.addWidget(hint)

        main.addWidget(quick)

        # ================= ADVANCED TOGGLE =================
        toggle_row = QHBoxLayout()
        toggle_row.addStretch()

        self.advanced_btn = QToolButton()
        self.advanced_btn.setText("⚙ Advanced options")
        self.advanced_btn.setCheckable(True)
        self.advanced_btn.setStyleSheet("font-size:10pt;")
        toggle_row.addWidget(self.advanced_btn)

        main.addLayout(toggle_row)

        # ================= ADVANCED PANEL =================
        self.advanced = Card("Advanced Manual Selection")
        self.advanced.setVisible(False)

        self.cb_pod5 = QCheckBox("POD5")
        self.cb_fastq = QCheckBox("FASTQ")
        self.cb_nanotel = QCheckBox("NanoTel")
        self.cb_align = QCheckBox("Alignment")
        self.cb_r = QCheckBox("R analysis")

        for cb in [self.cb_pod5, self.cb_fastq, self.cb_nanotel, self.cb_align, self.cb_r]:
            self.advanced.body.addWidget(cb)

        self.advanced.body.addSpacing(10)

        self.cb_filtration = QCheckBox("✓ NanoTel filtration")
        self.cb_mapping = QCheckBox("✓ Mapping analysis")
        self.cb_methylation = QCheckBox("✓ Methylation analysis")

        for cb in [self.cb_filtration, self.cb_mapping, self.cb_methylation]:
            cb.setEnabled(False)
            self.advanced.body.addWidget(cb)

        main.addWidget(self.advanced)

        # ================= RUN BUTTONS =================
        btn_row = QHBoxLayout()
        btn_row.setSpacing(10)

        self.run_btn = QPushButton("▶ Run Workflow")
        self.run_btn.setMinimumHeight(44)

        self.cancel_btn = QPushButton("✕ Cancel")
        self.cancel_btn.setMinimumHeight(44)
        self.cancel_btn.setEnabled(False)

        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.cancel_btn)
        main.addLayout(btn_row)

        # ================= LOG =================
        log = Card("Execution Log")
        self.log = QTextEdit()
        self.log.setFixedHeight(200)
        log.body.addWidget(self.log)
        main.addWidget(log)

        # ================= SIGNALS =================
        self.advanced_btn.toggled.connect(self.advanced.setVisible)
        self.cb_r.toggled.connect(self.toggle_r_opts)
        self.run_btn.clicked.connect(self.run)
        self.cancel_btn.clicked.connect(self.cancel)

        sys.stdout = StreamToGui(sys.stdout)
        sys.stdout.text_written.connect(self.append_log)

    # =============================================================================
    # METHODS
    # =============================================================================

    def browse(self, edit):
        path = QFileDialog.getExistingDirectory(self, "Select directory")
        if path:
            edit.setText(path)

    def toggle_r_opts(self, checked):
        self.cb_filtration.setEnabled(checked)
        self.cb_mapping.setEnabled(checked)
        self.cb_methylation.setEnabled(checked)

    def append_log(self, text):
        ts = datetime.now().strftime("%H:%M:%S")
        self.log.append(f"[{ts}] {text}")

    def apply_quick_mode(self):
        """Configure checkboxes based on quick mode selection"""
        mode = self.quick_mode.currentText()

        # Reset all
        for cb in [self.cb_pod5, self.cb_fastq, self.cb_nanotel, self.cb_align, self.cb_r]:
            cb.setChecked(False)

        # Set based on mode
        if mode == "Complete POD5 Workflow":
            self.cb_pod5.setChecked(True)
            self.cb_fastq.setChecked(True)
            self.cb_nanotel.setChecked(True)
            self.cb_align.setChecked(True)
        elif mode == "FASTQ Workflow":
            self.cb_fastq.setChecked(True)
            self.cb_nanotel.setChecked(True)
            self.cb_align.setChecked(True)
        elif mode == "NanoTel Analysis":
            self.cb_nanotel.setChecked(True)
        elif mode == "R Analysis":
            self.cb_r.setChecked(True)

    def run(self):
        """Validate and run the workflow"""
        # ===== VALIDATION: Output directory =====
        if not self.output_dir.text():
            QMessageBox.warning(
                self,
                " Missing output directory",
                "Please select an output directory."
            )
            return

        output_path = Path(self.output_dir.text())
        if not output_path.exists():
            QMessageBox.warning(
                self,
                " Invalid output directory",
                f"Output directory does not exist: {self.output_dir.text()}"
            )
            return

        # ===== VALIDATION: Input directories content =====
        input_errors = validate_input_directories(
            pod5=self.pod5.text(),
            fastq=self.fastq.text(),
            bam=self.bam.text(),
            nanotel=self.nanotel.text(),
        )

        if input_errors:
            QMessageBox.critical(
                self,
                "⚠ Invalid input directories",
                "Some input directories are invalid:\n\n• "
                + "\n• ".join(input_errors),
            )
            return

        # ===== DETERMINE WORKFLOW CONFIGURATION =====
        if self.advanced.isVisible():
            # Advanced mode: use manual selection
            do_pod5 = self.cb_pod5.isChecked()
            do_fastq = self.cb_fastq.isChecked()
            do_nanotel = self.cb_nanotel.isChecked()
            do_align = self.cb_align.isChecked()
            do_r = self.cb_r.isChecked()

            # Validate advanced selection
            errors = validate_advanced_selection(
                do_pod5=do_pod5,
                do_fastq=do_fastq,
                do_nanotel=do_nanotel,
                do_align=do_align,
                do_r=do_r,
                do_methylation=self.cb_methylation.isChecked(),
                pod5_path=self.pod5.text(),
                fastq_path=self.fastq.text(),
                bam_path=self.bam.text(),
                nanotel_path=self.nanotel.text(),
            )

            if errors:
                QMessageBox.critical(
                    self,
                    "⚠ Invalid configuration",
                    "The selected workflow cannot be run:\n\n• "
                    + "\n• ".join(errors),
                )
                return

            run_filtration = self.cb_filtration.isChecked()
            run_mapping = self.cb_mapping.isChecked()
            run_methylation = self.cb_methylation.isChecked()
        else:
            # Quick mode: auto-configure based on selection
            self.apply_quick_mode()
            do_pod5 = self.cb_pod5.isChecked()
            do_fastq = self.cb_fastq.isChecked()
            do_nanotel = self.cb_nanotel.isChecked()
            do_align = self.cb_align.isChecked()
            do_r = self.cb_r.isChecked()

            # Default R-analysis options
            run_filtration = True
            run_mapping = True
            run_methylation = True

        # ===== VALIDATION: At least one workflow step selected =====
        if not any([do_pod5, do_fastq, do_nanotel, do_align, do_r]):
            QMessageBox.warning(
                self,
                "⚠ No workflow selected",
                "Please select at least one workflow step."
            )
            return

        # ===== ALL VALIDATION PASSED - START WORKFLOW =====
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.log.clear()

        self.append_log("=" * 80)
        self.append_log("WORKFLOW STARTING")
        self.append_log("=" * 80)

        # Create worker thread
        self.thread = QThread()
        self.worker = WorkerThread(
            trial_name=output_path.name,
            pod5_path=self.pod5.text(),
            fastq_path=self.fastq.text(),
            bam_path=self.bam.text(),
            nanotel_path=self.nanotel.text(),
            output_dir=self.output_dir.text(),
            organism=self.organism_combo.currentText(),
            do_pod5=do_pod5,
            do_fastq=do_fastq,
            do_nanotel=do_nanotel,
            do_align=do_align,
            do_r=do_r,
            run_filtration=run_filtration,
            run_mapping=run_mapping,
            run_methylation=run_methylation,
        )

        # Move worker to thread
        self.worker.moveToThread(self.thread)

        # Connect signals
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.append_log)
        self.worker.done.connect(self.on_done)
        self.worker.done.connect(self.thread.quit)
        self.thread.finished.connect(self.thread.deleteLater)

        # Start thread
        self.thread.start()

    def cancel(self):
        """Cancel running workflow"""
        if not self.worker:
            return

        self.append_log("⚠ Cancelling workflow...")

        self.worker.stop()
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(False)

    def on_done(self, ok, msg):
        """Handle workflow completion or cancellation"""

        # Re-enable UI
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)

        # Cleanup worker & thread
        if self.worker:
            self.worker.deleteLater()
            self.worker = None

        if self.thread:
            self.thread.quit()
            self.thread.wait()
            self.thread.deleteLater()
            self.thread = None

        if ok:
            self.append_log("Workflow finished successfully")
            QMessageBox.information(self, "Success", msg)
        else:
            self.append_log("⚠ Workflow stopped")
            QMessageBox.information(self, "Stopped", msg)


# =============================================================================
# ENTRY
# =============================================================================

def main():
    app = QApplication(sys.argv)
    w = AppWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

