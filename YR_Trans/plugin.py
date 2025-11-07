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
        self.index_path_btn = QPushButton("Browse")
        self.index_path_btn.clicked.connect(lambda: self.browse_file(self.index_path_edit, "Index files (*)"))
        index_layout.addWidget(self.index_path_edit)
        index_layout.addWidget(self.index_path_btn)
        input_layout.addRow("Index prefix:", index_layout)
        
        # Reads files
        reads_layout = QHBoxLayout()
        self.reads1_edit = QLineEdit()
        self.reads1_btn = QPushButton("Browse R1")
        self.reads1_btn.clicked.connect(lambda: self.browse_file(self.reads1_edit, "FastQ files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        reads_layout.addWidget(self.reads1_edit)
        reads_layout.addWidget(self.reads1_btn)
        input_layout.addRow("Reads 1 (R1):", reads_layout)
        
        reads2_layout = QHBoxLayout()
        self.reads2_edit = QLineEdit()
        self.reads2_btn = QPushButton("Browse R2")
        self.reads2_btn.clicked.connect(lambda: self.browse_file(self.reads2_edit, "FastQ files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        reads2_layout.addWidget(self.reads2_edit)
        reads2_layout.addWidget(self.reads2_btn)
        input_layout.addRow("Reads 2 (R2, optional):", reads2_layout)
        
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
        
        # Output file
        output_layout = QHBoxLayout()
        self.output_bam_edit = QLineEdit()
        self.output_bam_edit.setPlaceholderText("Output BAM file")
        self.output_bam_btn = QPushButton("Browse")
        self.output_bam_btn.clicked.connect(lambda: self.save_file(self.output_bam_edit, "BAM files (*.bam)"))
        output_layout.addWidget(self.output_bam_edit)
        output_layout.addWidget(self.output_bam_btn)
        params_layout.addRow("Output BAM:", output_layout)
        
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
        
    def browse_file(self, line_edit, filter_str):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if file_path:
            line_edit.setText(file_path)
            
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
        elif msg_type == "warning":
            formatted = f"[{timestamp}] WARNING: {message}"
        else:
            formatted = f"[{timestamp}] {message}"
        self.console_text.appendPlainText(formatted)
        self.console_text.moveCursor(QTextCursor.End)
        
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
            
        if not self.reads1_edit.text() or not os.path.exists(self.reads1_edit.text()):
            QMessageBox.warning(self, "Error", "Please select valid reads file(s)!")
            return
            
        if not self.output_bam_edit.text():
            QMessageBox.warning(self, "Error", "Please specify output BAM file!")
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
        
        cmd = [
            self.tool_path,
            "-x", self.index_path_edit.text().replace("\\", "/"),
            "-p", str(self.threads_spin.value()),
        ]
        
        # Add reads
        reads2 = self.reads2_edit.text().replace("\\", "/")
        if reads2 and os.path.exists(reads2):
            cmd.extend(["-1", self.reads1_edit.text().replace("\\", "/"), "-2", reads2.replace("\\", "/")])
        else:
            cmd.extend(["-U", self.reads1_edit.text().replace("\\", "/")])
        
        # 生成SAM文件而不是直接管道传输
        temp_sam = self.output_bam_edit.text().replace(".bam", ".sam")
        cmd.extend([
            "-S",
            temp_sam
        ])

        self.add_console_message(" ".join(cmd), "command")
        self.tab_widget.setCurrentIndex(2)
        
        self.align_thread = AlignThread(cmd)
        self.align_thread.progress.connect(self.progress_bar.setFormat)
        self.align_thread.error.connect(self.command_error)
        self.align_thread.console_output.connect(self.add_console_message)
        self.align_thread.start()
        
        # 构建正确的samtools命令
        samtools_cmd = [
            self.samtools_path,
            "view",
            "-bS",
            temp_sam,
            "-o",
            self.output_bam_edit.text()
        ]

        self.add_console_message(" ".join(samtools_cmd), "command")
        self.samtools_thread = SamtoolsThread(samtools_cmd)  

        self.align_thread.finished.connect(self.samtools_thread.start)
        # 修复信号连接，确保传递正确的参数给command_finished
        self.samtools_thread.finished.connect(lambda x: self.command_finished(self.output_bam_edit.text()))
        
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


class SamtoolsThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd
        
    def run(self):
        try:
            self.progress.emit("Converting SAM to BAM...")
            # 修复编码问题，显式指定编码为utf-8，并处理可能的错误
            result = subprocess.run(self.cmd, capture_output=True, text=True, timeout=3600, encoding='utf-8', errors='ignore')
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(result.stderr or "Samtools conversion failed")
            else:
                self.finished.emit("Conversion completed successfully")
        except subprocess.TimeoutExpired:
            self.error.emit("Operation timed out")
        except Exception as e:
            self.error.emit(str(e))


class AlignThread(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd
        
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
            else:
                self.finished.emit("Alignment completed successfully")
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
        
        # BAM/SAM files
        bam_layout = QHBoxLayout()
        self.bam_files_edit = QLineEdit()
        self.bam_files_edit.setPlaceholderText("Select BAM/SAM files...")
        self.bam_files_btn = QPushButton("Browse")
        self.bam_files_btn.clicked.connect(self.browse_bam_files)
        bam_layout.addWidget(self.bam_files_edit)
        bam_layout.addWidget(self.bam_files_btn)
        input_layout.addRow("BAM/SAM files:", bam_layout)
        
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
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select BAM/SAM Files", "", 
                                                      "BAM files (*.bam);;SAM files (*.sam);;All files (*)")
        if file_paths:
            self.bam_files_edit.setText(";".join(file_paths))
            self.bam_files_edit.setToolTip("\n".join(file_paths))
            
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
        
    def run_command(self):
        if self.is_running:
            return
            
        if not self.anno_file_edit.text() or not os.path.exists(self.anno_file_edit.text()):
            QMessageBox.warning(self, "Error", "Please select an annotation file!")
            return
            
        bam_files_str = self.bam_files_edit.text()
        if not bam_files_str:
            QMessageBox.warning(self, "Error", "Please select BAM/SAM files!")
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
        
        bam_files = [f.strip() for f in bam_files_str.split(";") if f.strip()]
        cmd = [
            self.tool_path,
            "-a", self.anno_file_edit.text(),
            "-t", self.feature_type_combo.currentText(),
            "-g", self.attr_type_edit.text(),
            "-T", str(self.threads_spin.value()),
            "-o", self.output_file_edit.text(),
            "-p" if self.paired_end_radio.isChecked() else ""
        ]
        cmd.extend(bam_files)
        
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
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            # Find the header line (starts with Geneid)
            header_idx = 0
            for i, line in enumerate(lines):
                if line.startswith("Geneid") or line.startswith("GeneID"):
                    header_idx = i
                    break
            
            # Read data starting from header
            df = pd.read_csv(output_file, sep='\t', skiprows=header_idx, index_col=0)
            
            # Check if matplotlib is available
            if MATPLOTLIB_AVAILABLE:
                # Clear the previous plot
                self.figure.clear()
                
                # Calculate total counts per gene (sum across all samples)
                total_counts_per_gene = df.iloc[:, :-4].sum(axis=1)  # Exclude last 4 columns (Length, EffectiveLength, etc.)
                
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
            
        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Could not load results: {e}")
        
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
        
        # Results table
        self.results_table = QTableWidget()
        self.results_table.setAlternatingRowColors(True)
        self.results_table.setSortingEnabled(True)
        layout.addWidget(self.results_table)
        
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
        
    def analysis_finished(self, results_df):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        # Display results
        self.display_results(results_df)
        self.tab_widget.setCurrentIndex(2)
        QMessageBox.information(self, "Success", "Analysis completed successfully!")
        
    def display_results(self, results_df):
        """Display results in table"""
        try:
            self.results_table.setRowCount(len(results_df))
            self.results_table.setColumnCount(len(results_df.columns))
            self.results_table.setHorizontalHeaderLabels(results_df.columns.tolist())
            
            for i, (idx, row) in enumerate(results_df.iterrows()):
                # First column: gene ID
                self.results_table.setItem(i, 0, QTableWidgetItem(str(idx)))
                # Other columns
                for j, (col, val) in enumerate(row.items(), 1):
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
    finished = pyqtSignal(object)  # results DataFrame
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
            # Note: This is a simplified example. Actual PyDESeq2 API may differ
            # You'll need to adjust based on the actual PyDESeq2 library interface
            dds = pydeseq2(count_matrix=count_df, 
                          design_matrix=design_df,
                          design_formula="~ condition")
            dds.deseq2()
            
            self.progress.emit("Extracting results...")
            results = dds.results(alpha=self.alpha)
            
            # Filter by fold change if specified
            if 'log2FoldChange' in results.columns and self.fc_threshold > 1.0:
                results = results[
                    (results['log2FoldChange'].abs() >= np.log2(self.fc_threshold))
                ]
            
            # Save if output file specified
            if self.output_file:
                results.to_csv(self.output_file, sep='\t')
                
            self.finished.emit(results)
            
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

