"""
YR-Trans: Transcriptome Analysis Tools Plugin
Includes HISAT2, FeatureCounts, PyDESeq2, and integrated workflow
"""
from PyQt5.QtWidgets import (QWidget, QDialog, QVBoxLayout, QTabWidget, QHBoxLayout, 
                             QPushButton, QCheckBox, QPlainTextEdit, QProgressBar, 
                             QMessageBox, QFileDialog, QFrame, QLabel, QTextEdit, 
                             QComboBox, QSpinBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSizePolicy, QScrollArea, QTableWidget, QTableWidgetItem,
                             QDoubleSpinBox, QSplitter, QButtonGroup, QRadioButton)
from PyQt5.QtCore import QThread, pyqtSignal, QUrl, Qt
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtGui import QFont, QTextCursor
import os
import json
import tempfile
import platform
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available. Install it to view histogram plots.")

# Try to import PyDESeq2 if available
try:
    import pydeseq2
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False
    print("Warning: PyDESeq2 not available. Install it to use differential expression analysis.")


# ============================================================================
# Template and Base Classes (adapted from YR-MPE)
# ============================================================================

class CommandThread(QThread):
    """Base command thread for executing tools"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, list)  # output_files, html_files
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)  # message, msg_type
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__()
        self.tool_path = tool_path
        self.input_files = input_files if isinstance(input_files, list) else [input_files]
        self.parameters = parameters
        self.imported_files = imported_files or []
    
    def get_tool_name(self):
        return "Tool"
        
    def run(self):
        # Override in subclasses
        pass


# ============================================================================
# HISAT2 Wrapper
# ============================================================================

class HISAT2_wrapper(QWidget):
    """HISAT2 wrapper for RNA-seq alignment"""
    
    def __init__(self, config=None, plugin_path=None):
        super().__init__()
        self.config = config
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.tool_path = None
        self.build_path = None
        self.is_running = False
        self.load_config()

        print(f"Plugin path: {self.plugin_path}")   

        print(f"HISAT2 tool path: {self.tool_path}")
        print(f"HISAT2 build path: {self.build_path}")

        self.init_ui()
        
    def load_config(self):
        """Load HISAT2 and samtools paths from config.json"""
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                for tool in config_data:
                    if tool.get("name") == "HISAT2":
                        self.tool_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                        if "build_path" in tool:
                            self.build_path = os.path.join(self.plugin_path, tool["build_path"].lstrip("/"))
                    elif tool.get("name") == "samtools":
                        self.samtools_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                        # break
        except Exception as e:
            print(f"Config load failed: {e}")
    
    def init_ui(self):
        self.setWindowTitle("HISAT2 Alignment Tool")
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Mode selection
        mode_group = QGroupBox("Operation Mode")
        mode_layout = QVBoxLayout()
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Build Index", "Align Reads"])
        self.mode_combo.currentIndexChanged.connect(self.on_mode_changed)
        mode_layout.addWidget(self.mode_combo)
        mode_group.setLayout(mode_layout)
        main_layout.addWidget(mode_group)
        
        # Tab widget for different modes
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Build Index Tab
        self.build_tab = QWidget()
        self.tab_widget.addTab(self.build_tab, "Build Index")
        self.setup_build_tab()
        
        # Align Tab
        self.align_tab = QWidget()
        self.tab_widget.addTab(self.align_tab, "Align Reads")
        self.setup_align_tab()
        
        # Console Tab
        self.console_tab = QWidget()
        self.tab_widget.addTab(self.console_tab, "Console")
        self.setup_console_tab()

        # TODO: As Transcriptome alignment is a big task, add a button to support 'run in dedicated terminal', 
        # & another button to support 'copy command to clipboard', 
        # & another button to support 'save command to file'.
        
        # Control panel
        self.setup_control_panel(main_layout)
        
    def on_mode_changed(self, index):
        """Switch tabs based on mode"""
        if index == 0:
            self.tab_widget.setCurrentIndex(0)
        else:
            self.tab_widget.setCurrentIndex(1)
    
    def setup_build_tab(self):
        layout = QVBoxLayout()
        self.build_tab.setLayout(layout)
        
        # Reference genome input
        ref_group = QGroupBox("Reference Genome")
        ref_layout = QFormLayout()
        ref_group.setLayout(ref_layout)
        
        ref_file_layout = QHBoxLayout()
        self.ref_genome_edit = QLineEdit()
        self.ref_genome_btn = QPushButton("Browse")
        self.ref_genome_btn.clicked.connect(lambda: self.browse_file(self.ref_genome_edit, "FASTA files (*.fa *.fasta *.fas)"))
        ref_file_layout.addWidget(self.ref_genome_edit)
        ref_file_layout.addWidget(self.ref_genome_btn)
        ref_layout.addRow("Genome file:", ref_file_layout)
        
        # Index prefix
        self.index_prefix_edit = QLineEdit()
        self.index_prefix_edit.setPlaceholderText("Output index prefix")
        ref_layout.addRow("Index prefix:", self.index_prefix_edit)
        
        layout.addWidget(ref_group)
        layout.addStretch()
        
    def setup_align_tab(self):
        layout = QVBoxLayout()
        self.align_tab.setLayout(layout)
        
        # Input files
        input_group = QGroupBox("Input Files")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        
        # Index file
        index_layout = QHBoxLayout()
        self.index_path_edit = QLineEdit()
        self.index_path_edit.setPlaceholderText("Select any HISAT2 index file or enter index prefix")
        self.index_path_btn = QPushButton("Browse")
        self.index_path_btn.clicked.connect(self.browse_index_file)
        index_layout.addWidget(self.index_path_edit)
        index_layout.addWidget(self.index_path_btn)
        input_layout.addRow("Index prefix:", index_layout)
        
        # Multiple reads files
        reads_group = QGroupBox("Reads Files (Support multiple samples)")
        reads_layout = QVBoxLayout()
        reads_group.setLayout(reads_layout)
        
        # Add file tags container for multiple samples
        self.reads_tags_container = QFrame()
        self.reads_tags_layout = QVBoxLayout()
        self.reads_tags_container.setLayout(self.reads_tags_layout)
        self.reads_tags_container.setVisible(False)
        reads_layout.addWidget(self.reads_tags_container)
        
        # Button to add reads files
        add_reads_btn = QPushButton("Add Reads Files")
        add_reads_btn.clicked.connect(self.add_reads_files)
        reads_layout.addWidget(add_reads_btn)
        
        input_layout.addRow(reads_group)
        
        layout.addWidget(input_group)
        
        # Parameters
        params_group = QGroupBox("Alignment Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        
        # Threads
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 32)
        self.threads_spin.setValue(4)
        params_layout.addRow("Threads (-p):", self.threads_spin)
        
        # Output directory
        output_layout = QHBoxLayout()
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory for BAM files")
        self.output_dir_btn = QPushButton("Browse")
        self.output_dir_btn.clicked.connect(self.browse_output_dir)
        output_layout.addWidget(self.output_dir_edit)
        output_layout.addWidget(self.output_dir_btn)
        params_layout.addRow("Output Directory:", output_layout)
        
        layout.addWidget(params_group)
        layout.addStretch()
        
    def setup_console_tab(self):
        layout = QVBoxLayout()
        self.console_tab.setLayout(layout)
        
        controls = QHBoxLayout()
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_console)
        controls.addWidget(clear_btn)
        controls.addStretch()
        layout.addLayout(controls)
        
        self.console_text = QPlainTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Courier New", 9))
        if platform.system() == "Windows":
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    font-family: 'Consolas', monospace;
                }
            """)
        layout.addWidget(self.console_text)
        
    def setup_control_panel(self, main_layout):
        control_layout = QHBoxLayout()
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_command)
        control_layout.addWidget(self.run_button)
        
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_command)
        self.stop_button.setEnabled(False)
        control_layout.addWidget(self.stop_button)
        
        main_layout.addLayout(control_layout)
        
    def browse_index_file(self):
        """Browse HISAT2 index files"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select HISAT2 Index File", "", 
                                                   "HISAT2 Index Files (*.ht2 *.ht1);;All files (*)")
        if file_path:
            # Extract the index prefix by removing the .N.ht2 or .N.ht1 suffix
            index_prefix = self.extract_index_prefix(file_path)
            self.index_path_edit.setText(index_prefix)
            
    def extract_index_prefix(self, file_path):
        """
        Extract the index prefix from a HISAT2 index file path.
        HISAT2 index files are named like: prefix.1.ht2, prefix.2.ht2, etc.
        We need to remove the .N.ht2 or .N.ht1 suffix to get the prefix.
        """
        # Get the file name without directory
        file_name = os.path.basename(file_path)
        
        # Check if it's a HISAT2 index file (.N.ht2 or .N.ht1)
        import re
        # Pattern to match suffixes like .1.ht2, .2.ht2, .1.ht1, .2.ht1, etc.
        pattern = r'\.\d+\.ht[12]$'
        match = re.search(pattern, file_name)
        
        if match:
            # Extract the prefix by removing the suffix
            prefix = file_name[:match.start()]
            # Return the full path with the prefix
            dir_path = os.path.dirname(file_path)
            return os.path.join(dir_path, prefix)
        else:
            # If it doesn't match the pattern, return the file path without extension
            return os.path.splitext(file_path)[0]
            
    def browse_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
            
    def browse_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir_edit.setText(dir_path)
            
    def clear_console(self):
        self.console_text.clear()
        
    def add_console_message(self, message, msg_type="info"):
        timestamp = __import__('datetime').datetime.now().strftime("%H:%M:%S")
        if msg_type == "command":
            formatted = f"[{timestamp}] $ {message}"
        elif msg_type == "error":
            formatted = f"[{timestamp}] ERROR: {message}"
        elif msg_type == "warning":
            formatted = f"[{timestamp}] WARNING: {message}"
        else:
            formatted = f"[{timestamp}] {message}"
        self.console_text.appendPlainText(formatted)
        self.console_text.moveCursor(QTextCursor.End)
        
    def add_reads_files(self):
        """Add multiple reads files (paired-end or single-end)"""
        # Create dialog for adding reads files
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Reads Files")
        dialog.setMinimumSize(500, 300)
        layout = QVBoxLayout()
        dialog.setLayout(layout)
        
        # R1 file
        r1_layout = QHBoxLayout()
        r1_layout.addWidget(QLabel("Reads 1 (R1):"))
        r1_edit = QLineEdit()
        r1_btn = QPushButton("Browse")
        r1_btn.clicked.connect(lambda: self.browse_file(r1_edit, "FastQ files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        r1_layout.addWidget(r1_edit)
        r1_layout.addWidget(r1_btn)
        layout.addLayout(r1_layout)
        
        # R2 file (optional)
        r2_layout = QHBoxLayout()
        r2_layout.addWidget(QLabel("Reads 2 (R2):"))
        r2_edit = QLineEdit()
        r2_btn = QPushButton("Browse")
        r2_btn.clicked.connect(lambda: self.browse_file(r2_edit, "FastQ files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        r2_layout.addWidget(r2_edit)
        r2_layout.addWidget(r2_btn)
        layout.addLayout(r2_layout)
        
        # Sample name
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("Sample Name:"))
        name_edit = QLineEdit()
        name_edit.setPlaceholderText("e.g., sample1")
        name_layout.addWidget(name_edit)
        layout.addLayout(name_layout)
        
        # Dialog buttons
        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("Add")
        cancel_btn = QPushButton("Cancel")
        btn_layout.addStretch()
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)
        
        ok_btn.clicked.connect(dialog.accept)
        cancel_btn.clicked.connect(dialog.reject)
        
        if dialog.exec_() == QDialog.Accepted:
            r1_file = r1_edit.text()
            r2_file = r2_edit.text()
            sample_name = name_edit.text()
            
            if not r1_file:
                QMessageBox.warning(self, "Warning", "Please select R1 reads file!")
                return
                
            if not sample_name:
                QMessageBox.warning(self, "Warning", "Please specify a sample name!")
                return
                
            # Add file tag for this sample
            self.add_reads_tag(r1_file, r2_file, sample_name)
            
    def add_reads_tag(self, r1_file, r2_file, sample_name):
        """Add a tag widget for a sample's reads files"""
        if not hasattr(self, 'reads_samples'):
            self.reads_samples = []
            
        # Store sample info
        sample_info = {
            'r1': r1_file,
            'r2': r2_file if r2_file else None,
            'name': sample_name
        }
        self.reads_samples.append(sample_info)
        
        # Create tag widget
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("""
            QFrame {
                background-color: #e9ecef;
                border-radius: 15px;
                margin: 2px;
            }
        """)
        
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        
        # Sample info label
        if sample_info['r2']:
            info_text = f"{sample_name}: {os.path.basename(r1_file)} + {os.path.basename(r2_file)}"
        else:
            info_text = f"{sample_name}: {os.path.basename(r1_file)}"
            
        info_label = QLabel(info_text)
        info_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(info_label)
        
        # Close button
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #dc3545;
                color: white;
                border: none;
                border-radius: 10px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_reads_tag(sample_info, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.reads_tags_layout.addWidget(tag_widget)
        if not hasattr(self, 'reads_tags'):
            self.reads_tags = []
        self.reads_tags.append((sample_info, tag_widget))
        
        # Show container
        self.reads_tags_container.setVisible(True)
        
    def remove_reads_tag(self, sample_info, tag_widget):
        """Remove a reads tag"""
        if hasattr(self, 'reads_samples') and sample_info in self.reads_samples:
            self.reads_samples.remove(sample_info)
        
        # Remove from tags list
        if hasattr(self, 'reads_tags'):
            self.reads_tags = [(si, tw) for si, tw in self.reads_tags if si != sample_info]
        
        # Remove widget
        self.reads_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Hide container if no samples left
        if not hasattr(self, 'reads_samples') or not self.reads_samples:
            self.reads_tags_container.setVisible(False)
            
    def run_command(self):
        if self.is_running:
            return
            
        mode = self.mode_combo.currentIndex()
        if mode == 0:  # Build index
            self.run_build_index()
        else:  # Align
            self.run_align()
            
    def run_build_index(self):
        if not self.ref_genome_edit.text() or not os.path.exists(self.ref_genome_edit.text()):
            QMessageBox.warning(self, "Error", "Please select a valid reference genome file!")
            return
            
        if not self.index_prefix_edit.text():
            QMessageBox.warning(self, "Error", "Please specify an index prefix!")
            return
            
        if not self.build_path or not os.path.exists(self.build_path):
            QMessageBox.critical(self, "Error", f"HISAT2-build executable not found!\nIn directory:{self.build_path}")
            return
            
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        cmd = [
            self.build_path,
            "-p", str(self.threads_spin.value()),
            self.ref_genome_edit.text(),
            self.index_prefix_edit.text()
        ]
        
        self.add_console_message(" ".join(cmd), "command")
        self.tab_widget.setCurrentIndex(2)  # Switch to console
        
        self.build_thread = BuildIndexThread(cmd)
        self.build_thread.progress.connect(self.progress_bar.setFormat)
        self.build_thread.finished.connect(self.command_finished)
        self.build_thread.error.connect(self.command_error)
        self.build_thread.console_output.connect(self.add_console_message)
        self.build_thread.start()
        
    def run_align(self):
        if not self.index_path_edit.text():
            QMessageBox.warning(self, "Error", "Please select an index file!")
            return
            
        if not hasattr(self, 'reads_samples') or not self.reads_samples:
            QMessageBox.warning(self, "Error", "Please add at least one sample with reads files!")
            return
            
        if not self.output_dir_edit.text():
            QMessageBox.warning(self, "Error", "Please specify output directory!")
            return
            
        if not os.path.exists(self.output_dir_edit.text()):
            QMessageBox.warning(self, "Error", "Output directory does not exist!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", f"HISAT2-align executable not found!\nIn directory: {self.tool_path}")
            return
            
        if not self.samtools_path or not os.path.exists(self.samtools_path):
            QMessageBox.critical(self, "Error", f"samtools executable not found!\nIn directory: {self.samtools_path}")
            return
            
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        # Run alignment for each sample
        self.sample_index = 0
        self.total_samples = len(self.reads_samples)
        self.align_next_sample()
        
    def align_next_sample(self):
        """Align next sample in the list"""
        if self.sample_index >= self.total_samples:
            # All samples processed
            self.command_finished("All samples aligned successfully")
            return
            
        sample = self.reads_samples[self.sample_index]
        self.progress_bar.setFormat(f"Processing sample {self.sample_index+1}/{self.total_samples}: {sample['name']}")
        
        # Create output BAM file path
        output_bam = os.path.join(self.output_dir_edit.text(), f"{sample['name']}.bam")
        
        cmd = [
            self.tool_path,
            "-x", self.index_path_edit.text().replace("\\", "/"),
            "-p", str(self.threads_spin.value()),
        ]
        
        # Add reads
        if sample['r2'] and os.path.exists(sample['r2']):
            cmd.extend(["-1", sample['r1'].replace("\\", "/"), "-2", sample['r2'].replace("\\", "/")])
        else:
            cmd.extend(["-U", sample['r1'].replace("\\", "/")])
        
        # Generate SAM file
        temp_sam = os.path.join(self.output_dir_edit.text(), f"{sample['name']}.sam")
        cmd.extend([
            "-S",
            temp_sam
        ])

        self.add_console_message(" ".join(cmd), "command")
        self.tab_widget.setCurrentIndex(2)
        
        self.align_thread = AlignThread(cmd, temp_sam, self.samtools_path, output_bam)
        self.align_thread.progress.connect(self.progress_bar.setFormat)
        self.align_thread.error.connect(self.command_error)
        self.align_thread.console_output.connect(self.add_console_message)
        self.align_thread.finished.connect(lambda: self.samtools_finished(output_bam))
        self.align_thread.start()
        
    def samtools_finished(self, output_bam):
        """Called when samtools conversion is finished for a sample"""
        self.sample_index += 1
        self.align_next_sample()
        
    def command_finished(self, output_file):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Command completed successfully!", "info")
        QMessageBox.information(self, "Success", "Operation completed successfully!")
        
    def command_error(self, error_msg):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"Operation failed.\nSee console for details.")
        self.tab_widget.setCurrentIndex(2)
        
    def stop_command(self):
        if hasattr(self, 'build_thread') and self.build_thread.isRunning():
            self.build_thread.terminate()
            self.build_thread.wait()
        if hasattr(self, 'align_thread') and self.align_thread.isRunning():
            self.align_thread.terminate()
            self.align_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)


class BuildIndexThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd
        
    def run(self):
        try:
            self.progress.emit("Building index...")
            result = subprocess.run(self.cmd, capture_output=True, text=True, timeout=3600)
            
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(result.stderr)
            else:
                self.finished.emit("Index built successfully")
        except subprocess.TimeoutExpired:
            self.error.emit("Operation timed out")
        except Exception as e:
            self.error.emit(str(e))


class AlignThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal()
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd, temp_sam, samtools_path, output_bam):
        super().__init__()
        self.cmd = cmd
        self.temp_sam = temp_sam
        self.samtools_path = samtools_path
        self.output_bam = output_bam
        
    def run(self):
        try:
            self.progress.emit("Aligning reads...")
            result = subprocess.run(self.cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
            
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(result.stderr)
                return
                
            # Run samtools to convert SAM to BAM
            self.progress.emit("Converting SAM to BAM...")
            samtools_cmd = [
                self.samtools_path,
                "view",
                "-bS",
                self.temp_sam,
                "-o",
                self.output_bam
            ]
            
            samtools_result = subprocess.run(samtools_cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
            
            if samtools_result.stdout:
                for line in samtools_result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if samtools_result.stderr:
                for line in samtools_result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if samtools_result.returncode != 0 else "output")
            
            # Clean up temp SAM file
            try:
                os.remove(self.temp_sam)
            except:
                pass
                
            if samtools_result.returncode != 0:
                self.error.emit(samtools_result.stderr or "Samtools conversion failed")
            else:
                self.finished.emit()
        except subprocess.TimeoutExpired:
            self.error.emit("Operation timed out")
        except Exception as e:
            self.error.emit(str(e))


class HISAT2_wrapper_entry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return HISAT2_wrapper(config=self.config, plugin_path=self.plugin_path)

# ============================================================================
# FeatureCounts Wrapper
# ============================================================================

class FeatureCounts_wrapper(QWidget):
    """FeatureCounts wrapper for read counting"""
    
    def __init__(self, config=None, plugin_path=None):
        super().__init__()
        self.config = config
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.tool_path = None
        self.is_running = False
        self.count_data = None  # To store count data for TPM/FPKM calculation
        self.expr_buttons_added = False  # To track if expression buttons are added
        self.load_config()
        self.init_ui()

        self.init_ui()
        
    def load_config(self):
        """Load HISAT2 and samtools paths from config.json"""
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                for tool in config_data:
                    if tool.get("name") == "HISAT2":
                        self.tool_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                        if "build_path" in tool:
                            self.build_path = os.path.join(self.plugin_path, tool["build_path"].lstrip("/"))
                    elif tool.get("name") == "samtools":
                        self.samtools_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                        # break
        except Exception as e:
            print(f"Config load failed: {e}")
    
    def init_ui(self):
        self.setWindowTitle("HISAT2 Alignment Tool")
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Mode selection
        mode_group = QGroupBox("Operation Mode")
        mode_layout = QVBoxLayout()
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Build Index", "Align Reads"])
        self.mode_combo.currentIndexChanged.connect(self.on_mode_changed)
        mode_layout.addWidget(self.mode_combo)
        mode_group.setLayout(mode_layout)
        main_layout.addWidget(mode_group)
        
        # Tab widget for different modes
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Build Index Tab
        self.build_tab = QWidget()
        self.tab_widget.addTab(self.build_tab, "Build Index")
        self.setup_build_tab()
        
        # Align Tab
        self.align_tab = QWidget()
        self.tab_widget.addTab(self.align_tab, "Align Reads")
        self.setup_align_tab()
        
        # Console Tab
        self.console_tab = QWidget()
        self.tab_widget.addTab(self.console_tab, "Console")
        self.setup_console_tab()

        # TODO: As Transcriptome alignment is a big task, add a button to support 'run in dedicated terminal', 
        # & another button to support 'copy command to clipboard', 
        # & another button to support 'save command to file'.
        
        # Control panel
        self.setup_control_panel(main_layout)
        
    def on_mode_changed(self, index):
        """Switch tabs based on mode"""
        if index == 0:
            self.tab_widget.setCurrentIndex(0)
        else:
            self.tab_widget.setCurrentIndex(1)
    
    def setup_build_tab(self):
        layout = QVBoxLayout()
        self.build_tab.setLayout(layout)
        
        # Reference genome input
        ref_group = QGroupBox("Reference Genome")
        ref_layout = QFormLayout()
        ref_group.setLayout(ref_layout)
        
        ref_file_layout = QHBoxLayout()
        self.ref_genome_edit = QLineEdit()
        self.ref_genome_btn = QPushButton("Browse")
        self.ref_genome_btn.clicked.connect(lambda: self.browse_file(self.ref_genome_edit, "FASTA files (*.fa *.fasta *.fas)"))
        ref_file_layout.addWidget(self.ref_genome_edit)
        ref_file_layout.addWidget(self.ref_genome_btn)
        ref_layout.addRow("Genome file:", ref_file_layout)
        
        # Index prefix
        self.index_prefix_edit = QLineEdit()
        self.index_prefix_edit.setPlaceholderText("Output index prefix")
        ref_layout.addRow("Index prefix:", self.index_prefix_edit)
        
        layout.addWidget(ref_group)
        layout.addStretch()
        
    def setup_align_tab(self):
        layout = QVBoxLayout()
        self.align_tab.setLayout(layout)
        
        # Input files
        input_group = QGroupBox("Input Files")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        
        # Index file
        index_layout = QHBoxLayout()
        self.index_path_edit = QLineEdit()
        self.index_path_edit.setPlaceholderText("Select any HISAT2 index file or enter index prefix")
        self.index_path_btn = QPushButton("Browse")
        self.index_path_btn.clicked.connect(self.browse_index_file)
        index_layout.addWidget(self.index_path_edit)
        index_layout.addWidget(self.index_path_btn)
        input_layout.addRow("Index prefix:", index_layout)
        
        # Multiple reads files
        reads_group = QGroupBox("Reads Files (Support multiple samples)")
        reads_layout = QVBoxLayout()
        reads_group.setLayout(reads_layout)
        
        # Add file tags container for multiple samples
        self.reads_tags_container = QFrame()
        self.reads_tags_layout = QVBoxLayout()
        self.reads_tags_container.setLayout(self.reads_tags_layout)
        self.reads_tags_container.setVisible(False)
        reads_layout.addWidget(self.reads_tags_container)
        
        # Button to add reads files
        add_reads_btn = QPushButton("Add Reads Files")
        add_reads_btn.clicked.connect(self.add_reads_files)
        reads_layout.addWidget(add_reads_btn)
        
        input_layout.addRow(reads_group)
        
        layout.addWidget(input_group)
        
        # Parameters
        params_group = QGroupBox("Alignment Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        
        # Threads
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 32)
        self.threads_spin.setValue(4)
        params_layout.addRow("Threads (-p):", self.threads_spin)
        
        # Output directory
        output_layout = QHBoxLayout()
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory for BAM files")
        self.output_dir_btn = QPushButton("Browse")
        self.output_dir_btn.clicked.connect(self.browse_output_dir)
        output_layout.addWidget(self.output_dir_edit)
        output_layout.addWidget(self.output_dir_btn)
        params_layout.addRow("Output Directory:", output_layout)
        
        layout.addWidget(params_group)
        layout.addStretch()
        
    def setup_console_tab(self):
        layout = QVBoxLayout()
        self.console_tab.setLayout(layout)
        
        controls = QHBoxLayout()
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_console)
        controls.addWidget(clear_btn)
        controls.addStretch()
        layout.addLayout(controls)
        
        self.console_text = QPlainTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Courier New", 9))
        if platform.system() == "Windows":
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    font-family: 'Consolas', monospace;
                }
            """)
        layout.addWidget(self.console_text)
        
    def setup_control_panel(self, main_layout):
        control_layout = QHBoxLayout()
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_command)
        control_layout.addWidget(self.run_button)
        
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_command)
        self.stop_button.setEnabled(False)
        control_layout.addWidget(self.stop_button)
        
        main_layout.addLayout(control_layout)
        
    def browse_index_file(self):
        """Browse HISAT2 index files"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select HISAT2 Index File", "", 
                                                   "HISAT2 Index Files (*.ht2 *.ht1);;All files (*)")
        if file_path:
            # Extract the index prefix by removing the .N.ht2 or .N.ht1 suffix
            index_prefix = self.extract_index_prefix(file_path)
            self.index_path_edit.setText(index_prefix)
            
    def extract_index_prefix(self, file_path):
        """
        Extract the index prefix from a HISAT2 index file path.
        HISAT2 index files are named like: prefix.1.ht2, prefix.2.ht2, etc.
        We need to remove the .N.ht2 or .N.ht1 suffix to get the prefix.
        """
        # Get the file name without directory
        file_name = os.path.basename(file_path)
        
        # Check if it's a HISAT2 index file (.N.ht2 or .N.ht1)
        import re
        # Pattern to match suffixes like .1.ht2, .2.ht2, .1.ht1, .2.ht1, etc.
        pattern = r'\.\d+\.ht[12]$'
        match = re.search(pattern, file_name)
        
        if match:
            # Extract the prefix by removing the suffix
            prefix = file_name[:match.start()]
            # Return the full path with the prefix
            dir_path = os.path.dirname(file_path)
            return os.path.join(dir_path, prefix)
        else:
            # If it doesn't match the pattern, return the file path without extension
            return os.path.splitext(file_path)[0]
            
    def browse_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
            
    def browse_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir_edit.setText(dir_path)
            
    def clear_console(self):
        self.console_text.clear()
        
    def add_console_message(self, message, msg_type="info"):
        timestamp = __import__('datetime').datetime.now().strftime("%H:%M:%S")
        if msg_type == "command":
            formatted = f"[{timestamp}] $ {message}"
        elif msg_type == "error":
            formatted = f"[{timestamp}] ERROR: {message}"
        elif msg_type == "warning":
            formatted = f"[{timestamp}] WARNING: {message}"
        else:
            formatted = f"[{timestamp}] {message}"
        self.console_text.appendPlainText(formatted)
        self.console_text.moveCursor(QTextCursor.End)
        
    def add_reads_files(self):
        """Add multiple reads files (paired-end or single-end)"""
        # Create dialog for adding reads files
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Reads Files")
        dialog.setMinimumSize(500, 300)
        layout = QVBoxLayout()
        dialog.setLayout(layout)
        
        # R1 file
        r1_layout = QHBoxLayout()
        r1_layout.addWidget(QLabel("Reads 1 (R1):"))
        r1_edit = QLineEdit()
        r1_btn = QPushButton("Browse")
        r1_btn.clicked.connect(lambda: self.browse_file(r1_edit, "FastQ files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        r1_layout.addWidget(r1_edit)
        r1_layout.addWidget(r1_btn)
        layout.addLayout(r1_layout)
        
        # R2 file (optional)
        r2_layout = QHBoxLayout()
        r2_layout.addWidget(QLabel("Reads 2 (R2):"))
        r2_edit = QLineEdit()
        r2_btn = QPushButton("Browse")
        r2_btn.clicked.connect(lambda: self.browse_file(r2_edit, "FastQ files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        r2_layout.addWidget(r2_edit)
        r2_layout.addWidget(r2_btn)
        layout.addLayout(r2_layout)
        
        # Sample name
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("Sample Name:"))
        name_edit = QLineEdit()
        name_edit.setPlaceholderText("e.g., sample1")
        name_layout.addWidget(name_edit)
        layout.addLayout(name_layout)
        
        # Dialog buttons
        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("Add")
        cancel_btn = QPushButton("Cancel")
        btn_layout.addStretch()
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)
        
        ok_btn.clicked.connect(dialog.accept)
        cancel_btn.clicked.connect(dialog.reject)
        
        if dialog.exec_() == QDialog.Accepted:
            r1_file = r1_edit.text()
            r2_file = r2_edit.text()
            sample_name = name_edit.text()
            
            if not r1_file:
                QMessageBox.warning(self, "Warning", "Please select R1 reads file!")
                return
                
            if not sample_name:
                QMessageBox.warning(self, "Warning", "Please specify a sample name!")
                return
                
            # Add file tag for this sample
            self.add_reads_tag(r1_file, r2_file, sample_name)
            
    def add_reads_tag(self, r1_file, r2_file, sample_name):
        """Add a tag widget for a sample's reads files"""
        if not hasattr(self, 'reads_samples'):
            self.reads_samples = []
            
        # Store sample info
        sample_info = {
            'r1': r1_file,
            'r2': r2_file if r2_file else None,
            'name': sample_name
        }
        self.reads_samples.append(sample_info)
        
        # Create tag widget
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("""
            QFrame {
                background-color: #e9ecef;
                border-radius: 15px;
                margin: 2px;
            }
        """)
        
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        
        # Sample info label
        if sample_info['r2']:
            info_text = f"{sample_name}: {os.path.basename(r1_file)} + {os.path.basename(r2_file)}"
        else:
            info_text = f"{sample_name}: {os.path.basename(r1_file)}"
            
        info_label = QLabel(info_text)
        info_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(info_label)
        
        # Close button
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #dc3545;
                color: white;
                border: none;
                border-radius: 10px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_reads_tag(sample_info, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.reads_tags_layout.addWidget(tag_widget)
        if not hasattr(self, 'reads_tags'):
            self.reads_tags = []
        self.reads_tags.append((sample_info, tag_widget))
        
        # Show container
        self.reads_tags_container.setVisible(True)
        
    def remove_reads_tag(self, sample_info, tag_widget):
        """Remove a reads tag"""
        if hasattr(self, 'reads_samples') and sample_info in self.reads_samples:
            self.reads_samples.remove(sample_info)
        
        # Remove from tags list
        if hasattr(self, 'reads_tags'):
            self.reads_tags = [(si, tw) for si, tw in self.reads_tags if si != sample_info]
        
        # Remove widget
        self.reads_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Hide container if no samples left
        if not hasattr(self, 'reads_samples') or not self.reads_samples:
            self.reads_tags_container.setVisible(False)
            
    def run_command(self):
        if self.is_running:
            return
            
        mode = self.mode_combo.currentIndex()
        if mode == 0:  # Build index
            self.run_build_index()
        else:  # Align
            self.run_align()
            
    def run_build_index(self):
        if not self.ref_genome_edit.text() or not os.path.exists(self.ref_genome_edit.text()):
            QMessageBox.warning(self, "Error", "Please select a valid reference genome file!")
            return
            
        if not self.index_prefix_edit.text():
            QMessageBox.warning(self, "Error", "Please specify an index prefix!")
            return
            
        if not self.build_path or not os.path.exists(self.build_path):
            QMessageBox.critical(self, "Error", f"HISAT2-build executable not found!\nIn directory:{self.build_path}")
            return
            
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        cmd = [
            self.build_path,
            "-p", str(self.threads_spin.value()),
            self.ref_genome_edit.text(),
            self.index_prefix_edit.text()
        ]
        
        self.add_console_message(" ".join(cmd), "command")
        self.tab_widget.setCurrentIndex(2)  # Switch to console
        
        self.build_thread = BuildIndexThread(cmd)
        self.build_thread.progress.connect(self.progress_bar.setFormat)
        self.build_thread.finished.connect(self.command_finished)
        self.build_thread.error.connect(self.command_error)
        self.build_thread.console_output.connect(self.add_console_message)
        self.build_thread.start()
        
    def run_align(self):
        if not self.index_path_edit.text():
            QMessageBox.warning(self, "Error", "Please select an index file!")
            return
            
        if not hasattr(self, 'reads_samples') or not self.reads_samples:
            QMessageBox.warning(self, "Error", "Please add at least one sample with reads files!")
            return
            
        if not self.output_dir_edit.text():
            QMessageBox.warning(self, "Error", "Please specify output directory!")
            return
            
        if not os.path.exists(self.output_dir_edit.text()):
            QMessageBox.warning(self, "Error", "Output directory does not exist!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", f"HISAT2-align executable not found!\nIn directory: {self.tool_path}")
            return
            
        if not self.samtools_path or not os.path.exists(self.samtools_path):
            QMessageBox.critical(self, "Error", f"samtools executable not found!\nIn directory: {self.samtools_path}")
            return
            
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        # Run alignment for each sample
        self.sample_index = 0
        self.total_samples = len(self.reads_samples)
        self.align_next_sample()
        
    def align_next_sample(self):
        """Align next sample in the list"""
        if self.sample_index >= self.total_samples:
            # All samples processed
            self.command_finished("All samples aligned successfully")
            return
            
        sample = self.reads_samples[self.sample_index]
        self.progress_bar.setFormat(f"Processing sample {self.sample_index+1}/{self.total_samples}: {sample['name']}")
        
        # Create output BAM file path
        output_bam = os.path.join(self.output_dir_edit.text(), f"{sample['name']}.bam")
        
        cmd = [
            self.tool_path,
            "-x", self.index_path_edit.text().replace("\\", "/"),
            "-p", str(self.threads_spin.value()),
        ]
        
        # Add reads
        if sample['r2'] and os.path.exists(sample['r2']):
            cmd.extend(["-1", sample['r1'].replace("\\", "/"), "-2", sample['r2'].replace("\\", "/")])
        else:
            cmd.extend(["-U", sample['r1'].replace("\\", "/")])
        
        # Generate SAM file
        temp_sam = os.path.join(self.output_dir_edit.text(), f"{sample['name']}.sam")
        cmd.extend([
            "-S",
            temp_sam
        ])

        self.add_console_message(" ".join(cmd), "command")
        self.tab_widget.setCurrentIndex(2)
        
        self.align_thread = AlignThread(cmd, temp_sam, self.samtools_path, output_bam)
        self.align_thread.progress.connect(self.progress_bar.setFormat)
        self.align_thread.error.connect(self.command_error)
        self.align_thread.console_output.connect(self.add_console_message)
        self.align_thread.finished.connect(lambda: self.samtools_finished(output_bam))
        self.align_thread.start()
        
    def samtools_finished(self, output_bam):
        """Called when samtools conversion is finished for a sample"""
        self.sample_index += 1
        self.align_next_sample()
        
    def command_finished(self, output_file):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Command completed successfully!", "info")
        QMessageBox.information(self, "Success", "Operation completed successfully!")
        
    def command_error(self, error_msg):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"Operation failed.\nSee console for details.")
        self.tab_widget.setCurrentIndex(2)
        
    def stop_command(self):
        if hasattr(self, 'build_thread') and self.build_thread.isRunning():
            self.build_thread.terminate()
            self.build_thread.wait()
        if hasattr(self, 'align_thread') and self.align_thread.isRunning():
            self.align_thread.terminate()
            self.align_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)


class BuildIndexThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd
        
    def run(self):
        try:
            self.progress.emit("Building index...")
            result = subprocess.run(self.cmd, capture_output=True, text=True, timeout=3600)
            
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(result.stderr)
            else:
                self.finished.emit("Index built successfully")
        except subprocess.TimeoutExpired:
            self.error.emit("Operation timed out")
        except Exception as e:
            self.error.emit(str(e))


class AlignThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal()
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd, temp_sam, samtools_path, output_bam):
        super().__init__()
        self.cmd = cmd
        self.temp_sam = temp_sam
        self.samtools_path = samtools_path
        self.output_bam = output_bam
        
    def run(self):
        try:
            self.progress.emit("Aligning reads...")
            result = subprocess.run(self.cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
            
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(result.stderr)
                return
                
            # Run samtools to convert SAM to BAM
            self.progress.emit("Converting SAM to BAM...")
            samtools_cmd = [
                self.samtools_path,
                "view",
                "-bS",
                self.temp_sam,
                "-o",
                self.output_bam
            ]
            
            samtools_result = subprocess.run(samtools_cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
            
            if samtools_result.stdout:
                for line in samtools_result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if samtools_result.stderr:
                for line in samtools_result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if samtools_result.returncode != 0 else "output")
            
            # Clean up temp SAM file
            try:
                os.remove(self.temp_sam)
            except:
                pass
                
            if samtools_result.returncode != 0:
                self.error.emit(samtools_result.stderr or "Samtools conversion failed")
            else:
                self.finished.emit()
        except subprocess.TimeoutExpired:
            self.error.emit("Operation timed out")
        except Exception as e:
            self.error.emit(str(e))


class HISAT2_wrapper_entry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return HISAT2_wrapper(config=self.config, plugin_path=self.plugin_path)


# ============================================================================
# FeatureCounts Wrapper
# ============================================================================

class FeatureCounts_wrapper(QWidget):
    """FeatureCounts wrapper for read counting"""
    
    def __init__(self, config=None, plugin_path=None):
        super().__init__()
        self.config = config
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.tool_path = None
        self.is_running = False
        self.load_config()
        self.init_ui()
        
        
    def load_config(self):
        """Load FeatureCounts path from config.json"""
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                for tool in config_data:
                    if tool.get("name").lower() == "featurecounts":
                        self.tool_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                        # QMessageBox.information(self, "Info", f"FeatureCounts executable found in: {self.tool_path}")
                        break
        except Exception as e:
            print(f"Config load failed: {e}")
    
    def init_ui(self):
        self.setWindowTitle("FeatureCounts Read Counting Tool")
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Tabs
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Input Tab
        self.input_tab = QWidget()
        self.tab_widget.addTab(self.input_tab, "Input & Parameters")
        self.setup_input_tab()
        
        # Results Tab
        self.results_tab = QWidget()
        self.tab_widget.addTab(self.results_tab, "Results")
        self.setup_results_tab()
        
        # Console Tab
        self.console_tab = QWidget()
        self.tab_widget.addTab(self.console_tab, "Console")
        self.setup_console_tab()
        
        # Control panel
        self.setup_control_panel(main_layout)
        
    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # Input files
        input_group = QGroupBox("Input Files")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        
        # Annotation file
        anno_layout = QHBoxLayout()
        self.anno_file_edit = QLineEdit()
        self.anno_file_btn = QPushButton("Browse")
        self.anno_file_btn.clicked.connect(lambda: self.browse_file(self.anno_file_edit, "GTF files (*.gtf);;SAF files (*.saf)"))
        anno_layout.addWidget(self.anno_file_edit)
        anno_layout.addWidget(self.anno_file_btn)
        input_layout.addRow("Annotation file:", anno_layout)
        
        # BAM/SAM files with multi-file support
        bam_group = QGroupBox("BAM/SAM Files (Support multiple samples)")
        bam_layout = QVBoxLayout()
        bam_group.setLayout(bam_layout)
        
        # Add file tags container for multiple BAM/SAM files
        self.bam_tags_container = QFrame()
        self.bam_tags_layout = QVBoxLayout()
        self.bam_tags_container.setLayout(self.bam_tags_layout)
        self.bam_tags_container.setVisible(False)
        bam_layout.addWidget(self.bam_tags_container)
        
        # Button to add BAM/SAM files
        add_bam_btn = QPushButton("Add BAM/SAM Files")
        add_bam_btn.clicked.connect(self.browse_bam_files)
        bam_layout.addWidget(add_bam_btn)
        
        input_layout.addRow(bam_group)
        
        layout.addWidget(input_group)
        
        # Parameters
        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        
        # Library type
        self.lib_type_choice = QButtonGroup()
        lib_type_layout = QHBoxLayout()
        self.single_end_radio = QRadioButton("Single-end")
        self.paired_end_radio = QRadioButton("Paired-end")
        self.lib_type_choice.addButton(self.single_end_radio)
        self.lib_type_choice.addButton(self.paired_end_radio)
        self.lib_type_choice.setExclusive(True)
        self.paired_end_radio.setChecked(True)  # Set default selection
        
        lib_type_layout.addWidget(self.single_end_radio)
        lib_type_layout.addWidget(self.paired_end_radio)
        params_layout.addRow("Library type:", lib_type_layout)
        
        # Feature type
        self.feature_type_combo = QComboBox()
        self.feature_type_combo.addItems(["exon", "gene", "CDS", "transcript"])
        self.feature_type_combo.setCurrentText("exon")
        params_layout.addRow("Feature type (-t):", self.feature_type_combo)
        
        # Attribute type
        self.attr_type_edit = QLineEdit()
        self.attr_type_edit.setText("gene_id")
        self.attr_type_edit.setPlaceholderText("e.g., gene_id")
        params_layout.addRow("Attribute type (-g):", self.attr_type_edit)
        
        # Threads
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 32)
        self.threads_spin.setValue(4)
        params_layout.addRow("Threads (-T):", self.threads_spin)
        
        # Output file
        output_layout = QHBoxLayout()
        self.output_file_edit = QLineEdit()
        self.output_file_edit.setPlaceholderText("Output count file")
        self.output_file_btn = QPushButton("Browse")
        self.output_file_btn.clicked.connect(lambda: self.save_file(self.output_file_edit, "Text files (*.txt)"))
        output_layout.addWidget(self.output_file_edit)
        output_layout.addWidget(self.output_file_btn)
        params_layout.addRow("Output file:", output_layout)
        
        layout.addWidget(params_group)
        layout.addStretch()
        
    def setup_results_tab(self):
        layout = QVBoxLayout()
        self.results_tab.setLayout(layout)
        
        # Create matplotlib figure for histogram
        if MATPLOTLIB_AVAILABLE:
            self.figure = Figure(figsize=(10, 6))
            self.canvas = FigureCanvas(self.figure)
            layout.addWidget(self.canvas)
        else:
            # Fallback to table if matplotlib is not available
            self.results_table = QTableWidget()
            self.results_table.setAlternatingRowColors(True)
            layout.addWidget(self.results_table)
        
    def setup_console_tab(self):
        layout = QVBoxLayout()
        self.console_tab.setLayout(layout)
        
        controls = QHBoxLayout()
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_console)
        controls.addWidget(clear_btn)
        controls.addStretch()
        layout.addLayout(controls)
        
        self.console_text = QPlainTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Courier New", 9))
        if platform.system() == "Windows":
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    font-family: 'Consolas', monospace;
                }
            """)
        layout.addWidget(self.console_text)
        
    def setup_control_panel(self, main_layout):
        control_layout = QHBoxLayout()
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        self.run_button = QPushButton("Run FeatureCounts")
        self.run_button.clicked.connect(self.run_command)
        control_layout.addWidget(self.run_button)
        
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_command)
        self.stop_button.setEnabled(False)
        control_layout.addWidget(self.stop_button)
        
        main_layout.addLayout(control_layout)
        
    def browse_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
            
    def browse_bam_files(self):
        """Browse and add multiple BAM/SAM files"""
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select BAM/SAM Files", "", 
                                                      "BAM files (*.bam);;SAM files (*.sam);;All files (*)")
        if file_paths:
            for file_path in file_paths:
                self.add_bam_tag(file_path)
            
    def save_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
    
    def clear_console(self):
        self.console_text.clear()
        
    def add_console_message(self, message, msg_type="info"):
        timestamp = __import__('datetime').datetime.now().strftime("%H:%M:%S")
        if msg_type == "command":
            formatted = f"[{timestamp}] $ {message}"
        elif msg_type == "error":
            formatted = f"[{timestamp}] ERROR: {message}"
        else:
            formatted = f"[{timestamp}] {message}"
        self.console_text.appendPlainText(formatted)
        self.console_text.moveCursor(QTextCursor.End)
        
    def add_bam_tag(self, file_path):
        """Add a BAM/SAM file tag widget"""
        if not hasattr(self, 'bam_files'):
            self.bam_files = []
            
        if file_path in self.bam_files:
            return  # Avoid duplicates
            
        self.bam_files.append(file_path)
        
        # Create tag widget
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("""
            QFrame {
                background-color: #e9ecef;
                border-radius: 15px;
                margin: 2px;
            }
        """)
        
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        
        # Get display name
        display_name = os.path.basename(file_path)
        
        # File name label
        name_label = QLabel(display_name)
        name_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(name_label)
        
        # Close button
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #dc3545;
                color: white;
                border: none;
                border-radius: 10px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_bam_tag(file_path, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.bam_tags_layout.addWidget(tag_widget)
        if not hasattr(self, 'bam_tags'):
            self.bam_tags = []
        self.bam_tags.append((file_path, tag_widget))
        
        # Show container
        self.bam_tags_container.setVisible(True)
        
    def remove_bam_tag(self, file_path, tag_widget):
        """Remove a BAM/SAM file tag"""
        if hasattr(self, 'bam_files') and file_path in self.bam_files:
            self.bam_files.remove(file_path)
        
        # Remove from tags list
        if hasattr(self, 'bam_tags'):
            self.bam_tags = [(fp, tw) for fp, tw in self.bam_tags if fp != file_path]
        
        # Remove widget
        self.bam_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Hide container if no files left
        if not hasattr(self, 'bam_files') or not self.bam_files:
            self.bam_tags_container.setVisible(False)
            
    def run_command(self):
        if self.is_running:
            return
            
        if not self.anno_file_edit.text() or not os.path.exists(self.anno_file_edit.text()):
            QMessageBox.warning(self, "Error", "Please select an annotation file!")
            return
            
        if not hasattr(self, 'bam_files') or not self.bam_files:
            QMessageBox.warning(self, "Error", "Please add at least one BAM/SAM file!")
            return
            
        if not self.output_file_edit.text():
            QMessageBox.warning(self, "Error", "Please specify output file!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "FeatureCounts executable not found!")
            return
            
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        cmd = [
            self.tool_path,
            "-a", self.anno_file_edit.text(),
            "-t", self.feature_type_combo.currentText(),
            "-g", self.attr_type_edit.text(),
            "-T", str(self.threads_spin.value()),
            "-o", self.output_file_edit.text(),
            "-p" if self.paired_end_radio.isChecked() else ""
        ]
        
        # Add all BAM/SAM files
        cmd.extend(self.bam_files)
        
        self.add_console_message(" ".join(cmd), "command")
        self.tab_widget.setCurrentIndex(2)
        
        self.count_thread = CountThread(cmd, self.output_file_edit.text())
        self.count_thread.progress.connect(self.progress_bar.setFormat)
        self.count_thread.finished.connect(self.command_finished)
        self.count_thread.error.connect(self.command_error)
        self.count_thread.console_output.connect(self.add_console_message)
        self.count_thread.start()
        
    def command_finished(self, output_file):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Counting completed successfully!", "info")
        
        # Load and display results
        try:
            self.load_results(output_file)
        except Exception as e:
            self.add_console_message(f"Failed to load results: {e}", "error")
            
        QMessageBox.information(self, "Success", "Counting completed successfully!")
        self.tab_widget.setCurrentIndex(1)  # Switch to results tab
        
    def load_results(self, output_file):
        """Load count results and display histogram"""
        try:
            # FeatureCounts output has comments at the beginning
            with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
            
            # Find the header line (starts with Geneid)
            header_idx = 0
            for i, line in enumerate(lines):
                if line.startswith("Geneid") or line.startswith("GeneID"):
                    header_idx = i
                    break
            
            # Read data starting from header
            df = pd.read_csv(output_file, sep='\t', skiprows=header_idx)
            
            # Store the dataframe for later use
            self.count_data = df
            
            # Check if matplotlib is available
            if MATPLOTLIB_AVAILABLE:
                # Clear the previous plot
                self.figure.clear()
                
                # Calculate total counts per gene (sum across all samples)
                # Make sure we're working with numeric data
                count_columns = df.columns[1:-4]  # Exclude first column (Geneid) and last 4 columns
                numeric_data = df[count_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
                total_counts_per_gene = numeric_data.sum(axis=1)
                
                # Create histogram
                ax = self.figure.add_subplot(111)
                ax.hist(total_counts_per_gene, bins=50, edgecolor='black')
                ax.set_xlabel('Count Value')
                ax.set_ylabel('Number of Genes')
                ax.set_title('Distribution of Gene Counts')
                # ax.grid(True, alpha=0.3)

                # disable xticks
                ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                
                # Refresh canvas
                self.canvas.draw()
            else:
                # Fallback to table display if matplotlib is not available
                self.results_table.setRowCount(len(df))
                self.results_table.setColumnCount(len(df.columns))
                self.results_table.setHorizontalHeaderLabels(df.columns.tolist())
                
                for i, row in enumerate(df.itertuples()):
                    for j, val in enumerate(row[1:], 0):
                        item = QTableWidgetItem(str(val))
                        self.results_table.setItem(i, j, item)
                        
                self.results_table.resizeColumnsToContents()
            
            # Add buttons for TPM/FPKM calculation if annotation file is available
            self.add_expression_buttons()
            
            # Add filter button
            self.add_filter_button()
            
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not load results: {str(e)}")
    
    def add_filter_button(self):
        """Add button for filtering zero-count rows"""
        # Check if we already have the filter button
        if hasattr(self, 'filter_button_added') and self.filter_button_added:
            return
            
        layout = self.results_tab.layout()
        
        # Create a horizontal layout for filter button
        filter_layout = QHBoxLayout()
        
        # Add a separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        layout.addWidget(separator)
        
        # Add label
        label = QLabel("Filter data:")
        filter_layout.addWidget(label)
        
        # Add filter button
        filter_btn = QPushButton("Filter Zero-Count Rows")
        filter_btn.clicked.connect(self.filter_zero_count_rows)
        filter_layout.addWidget(filter_btn)
        
        # Add layout to the results tab
        layout.addLayout(filter_layout)
        
        # Mark that we've added the filter button
        self.filter_button_added = True
    
    def filter_zero_count_rows(self):
        """Filter out rows where all count values are zero"""
        if not hasattr(self, 'count_data') or self.count_data is None:
            QMessageBox.warning(self, "Warning", "No count data available!")
            return
        
        try:
            # Get count columns (exclude Geneid column and last 4 columns)
            count_columns = self.count_data.columns[1:-4]
            
            # Convert count data to numeric
            count_data = self.count_data[count_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
            
            # Identify rows where all counts are zero
            zero_rows = (count_data == 0).all(axis=1)
            
            # Filter out zero-count rows
            filtered_df = self.count_data[~zero_rows].copy()
            
            # Save filtered data
            base_name = os.path.splitext(self.output_file_edit.text())[0]
            filtered_file = f"{base_name}_filtered.txt"
            
            filtered_df.to_csv(filtered_file, sep='\t', index=False)
            
            # Update stored count data
            self.count_data = filtered_df
            
            QMessageBox.information(self, "Success", f"Filtered data saved to:\n{filtered_file}\n\n"
                                                  f"Removed {sum(zero_rows)} zero-count rows.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to filter data:\n{str(e)}")

    def add_expression_buttons(self):
        """Add buttons for TPM/FPKM calculation"""
        # Check if we already have the buttons
        if hasattr(self, 'expr_buttons_added') and self.expr_buttons_added:
            return
            
        layout = self.results_tab.layout()
        
        # Create a horizontal layout for expression calculation buttons
        expr_layout = QHBoxLayout()
        
        # Add a separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        layout.addWidget(separator)
        
        # Add label
        label = QLabel("Calculate expression values:")
        expr_layout.addWidget(label)
        
        # Add TPM button
        tpm_btn = QPushButton("Calculate TPM")
        tpm_btn.clicked.connect(self.calculate_tpm)
        expr_layout.addWidget(tpm_btn)
        
        # Add FPKM button
        fpkm_btn = QPushButton("Calculate FPKM")
        fpkm_btn.clicked.connect(self.calculate_fpkm)
        expr_layout.addWidget(fpkm_btn)
        
        # Add layout to the results tab
        layout.addLayout(expr_layout)
        
        # Mark that we've added the buttons
        self.expr_buttons_added = True
    
    def calculate_tpm(self):
        """Calculate TPM values based on count data and annotation file"""
        if not hasattr(self, 'count_data') or self.count_data is None:
            QMessageBox.warning(self, "Warning", "No count data available!")
            return
            
        if not self.anno_file_edit.text() or not os.path.exists(self.anno_file_edit.text()):
            QMessageBox.warning(self, "Warning", "Annotation file is required for TPM calculation!")
            return
            
        try:
            # Calculate TPM
            tpm_df = self._calculate_expression_values("TPM")
            
            # Save to file
            base_name = os.path.splitext(self.output_file_edit.text())[0]
            tpm_file = f"{base_name}_tpm.txt"
            
            tpm_df.to_csv(tpm_file, sep='\t', index=False)
            
            # Also save filtered version
            filtered_tpm_df = self._filter_zero_count_rows(tpm_df)
            filtered_tpm_file = f"{base_name}_tpm_filtered.txt"
            filtered_tpm_df.to_csv(filtered_tpm_file, sep='\t', index=False)
            
            QMessageBox.information(self, "Success", f"TPM values saved to:\n{tpm_file}\n\n"
                                                   f"Filtered TPM values saved to:\n{filtered_tpm_file}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to calculate TPM values:\n{str(e)}")
    
    def calculate_fpkm(self):
        """Calculate FPKM values based on count data and annotation file"""
        if not hasattr(self, 'count_data') or self.count_data is None:
            QMessageBox.warning(self, "Warning", "No count data available!")
            return
            
        if not self.anno_file_edit.text() or not os.path.exists(self.anno_file_edit.text()):
            QMessageBox.warning(self, "Warning", "Annotation file is required for FPKM calculation!")
            return
            
        try:
            # Calculate FPKM
            fpkm_df = self._calculate_expression_values("FPKM")
            
            # Save to file
            base_name = os.path.splitext(self.output_file_edit.text())[0]
            fpkm_file = f"{base_name}_fpkm.txt"
            
            fpkm_df.to_csv(fpkm_file, sep='\t', index=False)
            
            # Also save filtered version
            filtered_fpkm_df = self._filter_zero_count_rows(fpkm_df)
            filtered_fpkm_file = f"{base_name}_fpkm_filtered.txt"
            filtered_fpkm_df.to_csv(filtered_fpkm_file, sep='\t', index=False)
            
            QMessageBox.information(self, "Success", f"FPKM values saved to:\n{fpkm_file}\n\n"
                                                   f"Filtered FPKM values saved to:\n{filtered_fpkm_file}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to calculate FPKM values:\n{str(e)}")
    
    def _filter_zero_count_rows(self, df):
        """Filter out rows where all count values are zero"""
        try:
            # Get count columns (exclude Geneid column and any metadata columns)
            count_columns = [col for col in df.columns if col != 'Geneid']
            
            # Convert count data to numeric
            count_data = df[count_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
            
            # Identify rows where all counts are zero
            zero_rows = (count_data == 0).all(axis=1)
            
            # Filter out zero-count rows
            filtered_df = df[~zero_rows].copy()
            
            return filtered_df
            
        except Exception as e:
            raise Exception(f"Failed to filter data:\n{str(e)}")

    def _calculate_expression_values(self, expr_type):
        """Calculate TPM or FPKM values"""
        # Get gene lengths from annotation file (GTF format)
        gene_lengths = self._get_gene_lengths_from_gtf()
        
        # Get count data (exclude Geneid column and last 4 columns)
        count_columns = self.count_data.columns[1:-4]  # Sample columns
        count_data = self.count_data[count_columns]
        gene_ids = self.count_data.iloc[:, 0]  # Gene IDs
        
        # Ensure count_data contains only numeric values
        count_data = count_data.apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Match gene lengths with count data
        lengths = []
        valid_genes = []
        valid_count_data = []
        
        for i, gene_id in enumerate(gene_ids):
            # Convert gene_id to string to match with gene_lengths keys
            gene_id_str = str(gene_id)
            if gene_id_str in gene_lengths:
                lengths.append(gene_lengths[gene_id_str])
                valid_genes.append(gene_id_str)
                valid_count_data.append(count_data.iloc[i])
        
        if not valid_genes:
            raise Exception("No matching genes found between count data and annotation file")
        
        lengths = np.array(lengths)
        valid_count_data = pd.DataFrame(valid_count_data, columns=count_columns)
        
        # Calculate expression values
        if expr_type == "TPM":
            # TPM calculation:
            # 1. Normalize counts by gene length (in kilobases) - RPK
            # 2. Normalize RPK values by the sum of all RPKs (per million) - TPM
            rpk = valid_count_data.div(lengths / 1000, axis=0)
            tpm = rpk.div(rpk.sum(axis=0) / 1e6, axis=1)
            result_df = tpm.copy()
            result_df.insert(0, "Geneid", valid_genes)
            
        elif expr_type == "FPKM":
            # FPKM calculation:
            # 1. Normalize counts by gene length (in kilobases)
            # 2. Normalize by total count (per million)
            total_counts = count_data.sum(axis=0)  # Total counts per sample
            rpk = valid_count_data.div(lengths / 1000, axis=0)
            fpkm = rpk.div(total_counts / 1e6, axis=1)
            result_df = fpkm.copy()
            result_df.insert(0, "Geneid", valid_genes)
            
        return result_df
    
    def _get_gene_lengths_from_gtf(self):
        """Extract gene lengths from GTF annotation file"""
        gene_lengths = {}
        
        try:
            with open(self.anno_file_edit.text(), 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('#'):
                        continue  # Skip comment lines
                        
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue  # Not a valid GTF line
                        
                    feature_type = parts[2]
                    if feature_type != self.feature_type_combo.currentText():
                        continue  # Only process selected feature type
                        
                    # Parse attributes column (column 9)
                    attributes = parts[8]
                    gene_id = None
                    
                    # Extract gene_id from attributes
                    if 'gene_id' in attributes:
                        # Handle different attribute formats
                        if 'gene_id "' in attributes:
                            # Format: gene_id "ENSG00000187634";
                            gene_id_match = attributes.split('gene_id "')[1].split('";')[0]
                            gene_id = gene_id_match
                        else:
                            # Format: gene_id ENSG00000187634;
                            gene_id_match = attributes.split('gene_id ')[1].split(';')[0]
                            gene_id = gene_id_match
                    
                    if not gene_id:
                        continue
                        
                    # Calculate length of feature
                    try:
                        start = int(parts[3])
                        end = int(parts[4])
                        length = end - start + 1
                        
                        # For gene_id, sum up lengths of all exons (or other features)
                        if gene_id in gene_lengths:
                            gene_lengths[gene_id] += length
                        else:
                            gene_lengths[gene_id] = length
                    except ValueError:
                        # Skip lines with invalid coordinates
                        continue
                        
        except Exception as e:
            raise Exception(f"Failed to parse annotation file: {str(e)}")
            
        return gene_lengths

    def command_error(self, error_msg):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"Counting failed: {error_msg}")
        
    def stop_command(self):
        if hasattr(self, 'count_thread') and self.count_thread.isRunning():
            self.count_thread.terminate()
            self.count_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)


class CountThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd, output_file):
        super().__init__()
        self.cmd = cmd
        self.output_file = output_file
        
    def run(self):
        try:
            self.progress.emit("Counting reads...")
            result = subprocess.run(self.cmd, capture_output=True, text=True, timeout=3600)
            
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(result.stderr)
            else:
                self.finished.emit(self.output_file)
        except subprocess.TimeoutExpired:
            self.error.emit("Operation timed out")
        except Exception as e:
            self.error.emit(str(e))


class FeatureCounts_wrapper_entry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return FeatureCounts_wrapper(config=self.config, plugin_path=self.plugin_path)


# ============================================================================
# PyDESeq2 Wrapper
# ============================================================================

class PyDESeq2_wrapper(QWidget):
    """PyDESeq2 wrapper for differential expression analysis"""
    
    def __init__(self, config=None, plugin_path=None):
        super().__init__()
        self.config = config
        self.plugin_path = plugin_path or os.path.dirname(os.path.abspath(__file__))
        self.is_running = False
        self.count_data = None
        self.current_design_df = None
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("PyDESeq2 Differential Expression Analysis")
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Tabs
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Input Tab
        self.input_tab = QWidget()
        self.tab_widget.addTab(self.input_tab, "Input Data")
        self.setup_input_tab()
        
        # Design Tab
        self.design_tab = QWidget()
        self.tab_widget.addTab(self.design_tab, "Experimental Design")
        self.setup_design_tab()
        
        # Results Tab
        self.results_tab = QWidget()
        self.tab_widget.addTab(self.results_tab, "Results")
        self.setup_results_tab()
        
        # Control panel
        self.setup_control_panel(main_layout)
        
    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # Count file input
        input_group = QGroupBox("Count Data")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        
        count_layout = QHBoxLayout()
        self.count_file_edit = QLineEdit()
        self.count_file_btn = QPushButton("Browse")
        self.count_file_btn.clicked.connect(lambda: self.browse_file(self.count_file_edit, "Count files (*.txt *.tsv *.csv)"))
        count_layout.addWidget(self.count_file_edit)
        count_layout.addWidget(self.count_file_btn)
        input_layout.addRow("Count file:", count_layout)
        
        layout.addWidget(input_group)
        
        # Info label
        info_label = QLabel("Note: Count file should have genes in rows and samples in columns. "
                           "The first column should be gene IDs, and the header row should contain sample names.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; padding: 10px;")
        layout.addWidget(info_label)
        layout.addStretch()
        
    def setup_design_tab(self):
        layout = QVBoxLayout()
        self.design_tab.setLayout(layout)
        
        # Design table
        design_group = QGroupBox("Sample Groups")
        design_layout = QVBoxLayout()
        design_group.setLayout(design_layout)
        
        info_label = QLabel("Specify which samples belong to which condition/group:")
        design_layout.addWidget(info_label)
        
        self.design_table = QTableWidget()
        self.design_table.setColumnCount(2)
        self.design_table.setHorizontalHeaderLabels(["Sample", "Condition"])
        self.design_table.horizontalHeader().setStretchLastSection(True)
        design_layout.addWidget(self.design_table)
        
        table_controls = QHBoxLayout()
        add_row_btn = QPushButton("Add Row")
        add_row_btn.clicked.connect(lambda: self.design_table.insertRow(self.design_table.rowCount()))
        remove_row_btn = QPushButton("Remove Row")
        remove_row_btn.clicked.connect(self.remove_design_row)
        table_controls.addWidget(add_row_btn)
        table_controls.addWidget(remove_row_btn)
        table_controls.addStretch()
        design_layout.addLayout(table_controls)
        
        layout.addWidget(design_group)
        
        # Parameters
        params_group = QGroupBox("Analysis Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.001, 0.1)
        self.alpha_spin.setSingleStep(0.01)
        self.alpha_spin.setValue(0.05)
        self.alpha_spin.setDecimals(3)
        params_layout.addRow("FDR threshold (alpha):", self.alpha_spin)
        
        self.fc_spin = QDoubleSpinBox()
        self.fc_spin.setRange(1.0, 10.0)
        self.fc_spin.setSingleStep(0.1)
        self.fc_spin.setValue(1.5)
        self.fc_spin.setDecimals(1)
        params_layout.addRow("Fold change threshold:", self.fc_spin)
        
        # Output file
        output_layout = QHBoxLayout()
        self.output_file_edit = QLineEdit()
        self.output_file_edit.setPlaceholderText("Output results file")
        self.output_file_btn = QPushButton("Browse")
        self.output_file_btn.clicked.connect(lambda: self.save_file(self.output_file_edit, "CSV files (*.csv)"))
        output_layout.addWidget(self.output_file_edit)
        output_layout.addWidget(self.output_file_btn)
        params_layout.addRow("Output file:", output_layout)
        
        layout.addWidget(params_group)
        layout.addStretch()
        
    def setup_results_tab(self):
        layout = QVBoxLayout()
        self.results_tab.setLayout(layout)
        
        # Summary label
        self.summary_label = QLabel("Results will be displayed here after analysis.")
        self.summary_label.setWordWrap(True)
        layout.addWidget(self.summary_label)
        
    def setup_control_panel(self, main_layout):
        control_layout = QHBoxLayout()
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        if not PYDESEQ2_AVAILABLE:
            self.run_button = QPushButton("Run PyDESeq2 (Not Available)")
            self.run_button.setEnabled(False)
            self.run_button.setToolTip("PyDESeq2 library is not installed. Please install it to use this feature.")
        else:
            self.run_button = QPushButton("Run PyDESeq2 Analysis")
            self.run_button.clicked.connect(self.run_analysis)
        control_layout.addWidget(self.run_button)
        
        main_layout.addLayout(control_layout)
        
    def browse_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
            self.load_count_file(file_path)
            
    def save_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
            
    def remove_design_row(self):
        current_row = self.design_table.currentRow()
        if current_row >= 0:
            self.design_table.removeRow(current_row)
            
    def load_count_file(self, file_path):
        """Load count file and populate design table"""
        try:
            # Try to read the count file
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            
            # Get sample names (columns)
            samples = df.columns.tolist()
            
            # Populate design table
            self.design_table.setRowCount(len(samples))
            for i, sample in enumerate(samples):
                self.design_table.setItem(i, 0, QTableWidgetItem(sample))
                self.design_table.setItem(i, 1, QTableWidgetItem("condition1"))
                
            QMessageBox.information(self, "Success", f"Loaded {len(samples)} samples. Please assign conditions in the Design tab.")
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not auto-load design: {e}")
            
    def run_analysis(self):
        if not PYDESEQ2_AVAILABLE:
            QMessageBox.critical(self, "Error", "PyDESeq2 library is not installed!")
            return
            
        if not self.count_file_edit.text() or not os.path.exists(self.count_file_edit.text()):
            QMessageBox.warning(self, "Error", "Please select a count file!")
            return
            
        # Get design from table
        design = {}
        for i in range(self.design_table.rowCount()):
            sample_item = self.design_table.item(i, 0)
            condition_item = self.design_table.item(i, 1)
            if sample_item and condition_item:
                design[sample_item.text()] = condition_item.text()
                
        if not design:
            QMessageBox.warning(self, "Error", "Please specify sample conditions in the Design tab!")
            return
            
        # Check that we have at least 2 conditions
        conditions = set(design.values())
        if len(conditions) < 2:
            QMessageBox.warning(self, "Error", "Please specify at least 2 different conditions!")
            return
            
        self.is_running = True
        self.run_button.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        # Run analysis in thread
        self.analysis_thread = DESeq2AnalysisThread(
            self.count_file_edit.text(),
            design,
            self.alpha_spin.value(),
            self.fc_spin.value(),
            self.output_file_edit.text() if self.output_file_edit.text() else None
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.start()
        
    def analysis_finished(self, results_df, design_df):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        # Store design_df for later use
        self.current_design_df = design_df
        
        # Display results
        self.display_results(results_df)
        self.tab_widget.setCurrentIndex(2)
        QMessageBox.information(self, "Success", "Analysis completed successfully!")
        
    def display_results(self, results_df):
        """Display results in table or DEGs distribution plot"""
        try:
            # Check if matplotlib is available
            if not MATPLOTLIB_AVAILABLE:
                self.display_results_table(results_df)
                return
                
            # Check if we have design information
            if not hasattr(self, 'current_design_df'):
                self.display_results_table(results_df)
                return
                
            # Extract conditions from design matrix
            conditions = list(set(self.current_design_df["condition"]))
            
            # Generate all pairwise comparisons
            comparisons = []
            for i in range(len(conditions)):
                for j in range(i + 1, len(conditions)):
                    comparisons.append((conditions[i], conditions[j]))
            
            # If we don't have comparison column, just show table
            if 'comparison' not in results_df.columns:
                self.display_results_table(results_df)
                return
            
            # Prepare data for plotting
            up_regulated = []
            down_regulated = []
            comparison_labels = []
            
            for cond1, cond2 in comparisons:
                comparison_name = f"{cond1}_vs_{cond2}"
                # Filter results for this comparison
                comp_results = results_df[results_df['comparison'] == comparison_name]
                
                # Filter significant genes with log2FoldChange > 0 (up-regulated) and < 0 (down-regulated)
                sig_genes = comp_results[comp_results['padj'] < self.alpha_spin.value()]
                up_count = len(sig_genes[sig_genes['log2FoldChange'] > 0])
                down_count = len(sig_genes[sig_genes['log2FoldChange'] < 0])
                
                up_regulated.append(up_count)
                down_regulated.append(down_count)
                comparison_labels.append(f"{cond1} vs {cond2}")
            
            # Create the plot
            self.figure = Figure(figsize=(12, 8))
            ax = self.figure.add_subplot(111)
            
            # Define colors
            up_color = '#ba3e45'  # Red for up-regulated
            down_color = '#4e6691'  # Blue for down-regulated
            
            # Plot bars
            x_pos = np.arange(len(comparison_labels))
            bar_width = 0.35
            
            bars_up = ax.bar(x_pos - bar_width/2, up_regulated, bar_width, color=up_color, label='Up-regulated')
            bars_down = ax.bar(x_pos + bar_width/2, down_regulated, bar_width, color=down_color, label='Down-regulated')
            
            # Add value labels on top of bars
            for i, (up_val, down_val) in enumerate(zip(up_regulated, down_regulated)):
                ax.text(i - bar_width/2, up_val + max(up_regulated + [1])*0.01, str(up_val), ha='center', va='bottom', fontsize=10)
                ax.text(i + bar_width/2, down_val + max(down_regulated + [1])*0.01, str(down_val), ha='center', va='bottom', fontsize=10)
            
            # Customize the plot
            ax.set_xlabel('Comparison Groups')
            ax.set_ylabel('Number of Differentially Expressed Genes')
            ax.set_title('Distribution of DEGs across All Pairwise Comparisons')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(comparison_labels, rotation=45, ha='right')
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
            
            # Adjust layout to prevent label cutoff
            self.figure.tight_layout()
            
            # Display the plot
            canvas = FigureCanvas(self.figure)
            layout = self.results_tab.layout()
            # Remove previous widgets
            for i in reversed(range(layout.count())): 
                layout.itemAt(i).widget().setParent(None)
            layout.addWidget(canvas)
            
            # Also show summary info
            total_sig = sum(up_regulated) + sum(down_regulated)
            summary = f"Total significant genes across all comparisons: {total_sig}\n"
            summary += f"Significant genes threshold (padj < {self.alpha_spin.value()}): {len(results_df[results_df['padj'] < self.alpha_spin.value()])}"
            self.summary_label.setText(summary)
            
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not display results: {e}")
            self.display_results_table(results_df)
    
    def display_results_table(self, results_df):
        """Fallback method to display results in table format"""
        try:
            # Clear existing layout
            layout = self.results_tab.layout()
            for i in reversed(range(layout.count())): 
                layout.itemAt(i).widget().setParent(None)
                
            # Add table widget
            self.results_table = QTableWidget()
            self.results_table.setAlternatingRowColors(True)
            self.results_table.setSortingEnabled(True)
            layout.addWidget(self.results_table)
            
            self.results_table.setRowCount(len(results_df))
            self.results_table.setColumnCount(len(results_df.columns))
            self.results_table.setHorizontalHeaderLabels(results_df.columns.tolist())
            
            for i, row in results_df.iterrows():
                # All columns including gene names
                for j, (col, val) in enumerate(row.items()):
                    item = QTableWidgetItem(str(val))
                    if isinstance(val, (int, float)):
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    self.results_table.setItem(i, j, item)
                    
            self.results_table.resizeColumnsToContents()
            
            # Summary
            if 'padj' in results_df.columns:
                sig_genes = results_df[results_df['padj'] < self.alpha_spin.value()].shape[0]
                up_genes = sig_genes if 'log2FoldChange' in results_df.columns else 0
                summary = f"Total genes: {len(results_df)}\n"
                summary += f"Significant genes (padj < {self.alpha_spin.value()}): {sig_genes}"
                self.summary_label.setText(summary)
                
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not display results: {e}")
            
    def analysis_error(self, error_msg):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Error", f"Analysis failed: {error_msg}")


class DESeq2AnalysisThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(object, object)  # results DataFrame, design DataFrame
    error = pyqtSignal(str)
    
    def __init__(self, count_file, design, alpha, fc_threshold, output_file):
        super().__init__()
        self.count_file = count_file
        self.design = design
        self.alpha = alpha
        self.fc_threshold = fc_threshold
        self.output_file = output_file
        
    def run(self):
        try:
            if not PYDESEQ2_AVAILABLE:
                self.error.emit("PyDESeq2 library not available")
                return
                
            self.progress.emit("Loading count data...")
            # Load count data
            count_df = pd.read_csv(self.count_file, sep='\t', index_col=0)
            
            # Prepare design DataFrame
            design_df = pd.DataFrame(list(self.design.items()), columns=['sample', 'condition'])
            design_df.set_index('sample', inplace=True)
            
            # Filter samples that exist in both count and design
            common_samples = [s for s in count_df.columns if s in design_df.index]
            if not common_samples:
                self.error.emit("No common samples found between count file and design!")
                return
                
            count_df = count_df[common_samples]
            design_df = design_df.loc[common_samples]
            
            self.progress.emit("Running DESeq2 analysis...")
            
            # Create DESeq2 object and run analysis
            # Use the correct PyDESeq2 API
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
            
            # Transpose count data to have samples as rows and genes as columns for PyDESeq2
            counts_transposed = count_df.T
            
            # Create DeseqDataSet object
            dds = DeseqDataSet(
                counts=counts_transposed,
                metadata=design_df,
                design_factors="condition",
                quiet=False
            )
            
            # Run the differential expression analysis
            dds.deseq2()
            
            # Get all unique conditions
            conditions = list(set(design_df["condition"]))
            
            # If we have at least 2 conditions, perform pairwise comparisons
            if len(conditions) >= 2:
                # Create a list to store all comparison results
                all_results = []
                
                # Perform all pairwise comparisons
                for i in range(len(conditions)):
                    for j in range(i + 1, len(conditions)):
                        cond1 = conditions[i]
                        cond2 = conditions[j]
                        
                        # Create contrast and get results
                        stat_res = DeseqStats(dds, contrast=("condition", cond1, cond2))
                        stat_res.summary()
                        results = stat_res.results_df
                        
                        # Reset index to make gene names a column
                        results = results.reset_index()
                        # Rename the index column to 'geneid'
                        if 'index' in results.columns:
                            results = results.rename(columns={'index': 'geneid'})
                        
                        # Add comparison info to results
                        results['comparison'] = f"{cond1}_vs_{cond2}"
                        
                        # Filter by fold change if specified
                        if 'log2FoldChange' in results.columns and self.fc_threshold > 1.0:
                            results = results[
                                (results['log2FoldChange'].abs() >= np.log2(self.fc_threshold))
                            ]
                        
                        all_results.append(results)
                
                # Combine all results
                if all_results:
                    combined_results = pd.concat(all_results, ignore_index=True)
                else:
                    # Fallback to default results if no comparisons worked
                    results = dds.summary()
                    combined_results = results.reset_index()
                    if 'index' in combined_results.columns:
                        combined_results = combined_results.rename(columns={'index': 'geneid'})
            else:
                # If we can't make a proper contrast, just get the results
                results = dds.summary()
                combined_results = results.reset_index()
                if 'index' in combined_results.columns:
                    combined_results = combined_results.rename(columns={'index': 'geneid'})
            
            # Ensure geneid column is present
            if combined_results.index.name is not None and 'geneid' not in combined_results.columns:
                combined_results['geneid'] = combined_results.index
                combined_results = combined_results.reset_index(drop=True)
            
            # Reorder columns to put geneid first
            cols = combined_results.columns.tolist()
            if 'geneid' in cols:
                cols.insert(0, cols.pop(cols.index('geneid')))
                combined_results = combined_results[cols]
            
            # Save if output file specified
            if self.output_file:
                combined_results.to_csv(self.output_file, sep='\t', index=False)
                
            self.finished.emit(combined_results, design_df)
            
        except Exception as e:
            import traceback
            self.error.emit(f"{str(e)}\n{traceback.format_exc()}")


class PyDESeq2_wrapper_entry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return PyDESeq2_wrapper(config=self.config, plugin_path=self.plugin_path)


# ============================================================================
# YR-Trans Workflow (Integrated Pipeline)
# ============================================================================

class YR_Trans_workflow(QWidget):
    """Integrated workflow for transcriptome analysis"""
    
    def __init__(self, config=None, plugin_path=None):
        super().__init__()
        self.config = config
        self.plugin_path = plugin_path or os.path.dirname(os.path.abspath(__file__))
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("YR-Trans: Transcriptome Analysis Workflow")
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Workflow steps
        steps_group = QGroupBox("Analysis Workflow")
        steps_layout = QVBoxLayout()
        steps_group.setLayout(steps_layout)
        
        # Step 1: HISAT2
        step1_layout = QHBoxLayout()
        step1_btn = QPushButton("Step 1: Build Index & Align (HISAT2)")
        step1_btn.clicked.connect(self.open_hisat2)
        step1_layout.addWidget(step1_btn)
        step1_status = QLabel("⏳ Pending")
        step1_status.setStyleSheet("color: #666; font-weight: bold;")
        step1_layout.addWidget(step1_status)
        step1_layout.addStretch()
        steps_layout.addLayout(step1_layout)
        self.step1_status = step1_status
        
        # Step 2: FeatureCounts
        step2_layout = QHBoxLayout()
        step2_btn = QPushButton("Step 2: Count Reads (FeatureCounts)")
        step2_btn.clicked.connect(self.open_featurecounts)
        step2_layout.addWidget(step2_btn)
        step2_status = QLabel("⏳ Pending")
        step2_status.setStyleSheet("color: #666; font-weight: bold;")
        step2_layout.addWidget(step2_status)
        step2_layout.addStretch()
        steps_layout.addLayout(step2_layout)
        self.step2_status = step2_status
        
        # Step 3: PyDESeq2
        step3_layout = QHBoxLayout()
        step3_btn = QPushButton("Step 3: Differential Expression (PyDESeq2)")
        step3_btn.clicked.connect(self.open_pydeseq2)
        step3_layout.addWidget(step3_btn)
        step3_status = QLabel("⏳ Pending")
        step3_status.setStyleSheet("color: #666; font-weight: bold;")
        step3_layout.addWidget(step3_status)
        step3_layout.addStretch()
        steps_layout.addLayout(step3_layout)
        self.step3_status = step3_status
        
        main_layout.addWidget(steps_group)
        
        # Quick workflow button
        quick_workflow_btn = QPushButton("Run Complete Workflow")
        quick_workflow_btn.setStyleSheet("font-size: 14px; font-weight: bold; padding: 10px;")
        quick_workflow_btn.clicked.connect(self.run_complete_workflow)
        main_layout.addWidget(quick_workflow_btn)
        
        # Info label
        info_label = QLabel(
            "YR-Trans Workflow integrates HISAT2, FeatureCounts, and PyDESeq2 for complete "
            "transcriptome analysis. You can run each step individually or use the complete workflow. "
            "For the complete workflow, ensure you have all input files ready."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("padding: 10px; background-color: #f0f0f0; border-radius: 5px;")
        main_layout.addWidget(info_label)
        
        main_layout.addStretch()
        
    def open_hisat2(self):
        widget = HISAT2_wrapper(config=self.config, plugin_path=self.plugin_path)
        self._open_in_new_tab("HISAT2", widget)
        
    def open_featurecounts(self):
        widget = FeatureCounts_wrapper(config=self.config, plugin_path=self.plugin_path)
        self._open_in_new_tab("FeatureCounts", widget)
        
    def open_pydeseq2(self):
        widget = PyDESeq2_wrapper(config=self.config, plugin_path=self.plugin_path)
        self._open_in_new_tab("PyDESeq2", widget)
        
    def _open_in_new_tab(self, name, widget):
        # This would open in a new tab in the main application
        # For now, just show as a dialog
        dialog = QDialog(self)
        dialog.setWindowTitle(name)
        dialog.setMinimumSize(1200, 800)
        layout = QVBoxLayout()
        layout.addWidget(widget)
        dialog.setLayout(layout)
        dialog.exec_()
        
    def run_complete_workflow(self):
        QMessageBox.information(
            self, 
            "Complete Workflow", 
            "The complete workflow will guide you through each step sequentially. "
            "Please use the individual step buttons to proceed through the analysis."
        )


class YR_Trans_entry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return YR_Trans_workflow(config=self.config, plugin_path=self.plugin_path)

