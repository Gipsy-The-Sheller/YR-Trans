"""
TransHub: Transcriptome Data Management and Analysis Tool
Manages transcriptome projects and datasets with visualization capabilities
"""

# TODO: add input portal for GXF annotation

import os
import json
import pandas as pd
import re
import subprocess
from pathlib import Path
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidget,
                             QTableWidgetItem, QTabWidget, QDialog, QFormLayout, QLineEdit,
                             QFileDialog, QMessageBox, QLabel, QComboBox, QGroupBox,
                             QHeaderView, QApplication, QCheckBox, QDoubleSpinBox, QFrame,
                             QTextEdit, QSpinBox, QProgressBar, QSplitter)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QTextCursor


class AddFastqDialog(QDialog):
    """Dialog for adding FASTQ files for a sample"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add FASTQ Files")
        self.setModal(True)
        self.resize(500, 200)
        
        self.sample_name = ""
        self.r1_file = ""
        self.r2_file = ""
        
        self.setup_ui()
        
    def setup_ui(self):
        layout = QFormLayout()
        
        # Sample name
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("e.g., sample1")
        layout.addRow("Sample Name:", self.name_edit)
        
        # R1 file
        r1_layout = QHBoxLayout()
        self.r1_edit = QLineEdit()
        self.r1_btn = QPushButton("Browse")
        self.r1_btn.clicked.connect(lambda: self.browse_file(self.r1_edit, "FASTQ Files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        r1_layout.addWidget(self.r1_edit)
        r1_layout.addWidget(self.r1_btn)
        layout.addRow("Reads 1 (R1):", r1_layout)
        
        # R2 file
        r2_layout = QHBoxLayout()
        self.r2_edit = QLineEdit()
        self.r2_btn = QPushButton("Browse")
        self.r2_btn.clicked.connect(lambda: self.browse_file(self.r2_edit, "FASTQ Files (*.fq *.fastq *.fq.gz *.fastq.gz)"))
        r2_layout.addWidget(self.r2_edit)
        r2_layout.addWidget(self.r2_btn)
        layout.addRow("Reads 2 (R2):", r2_layout)
        
        # Buttons
        button_layout = QHBoxLayout()
        ok_btn = QPushButton("Add")
        cancel_btn = QPushButton("Cancel")
        button_layout.addStretch()
        button_layout.addWidget(ok_btn)
        button_layout.addWidget(cancel_btn)
        layout.addRow(button_layout)
        
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def browse_file(self, line_edit, filter_str):
        file, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if file:
            line_edit.setText(file)
            
    def accept(self):
        self.sample_name = self.name_edit.text().strip()
        self.r1_file = self.r1_edit.text().strip()
        self.r2_file = self.r2_edit.text().strip()
        
        if not self.sample_name:
            QMessageBox.warning(self, "Warning", "Please enter sample name")
            return
            
        if not self.r1_file:
            QMessageBox.warning(self, "Warning", "Please select R1 file")
            return
            
        super().accept()


class SampleDesignDialog(QDialog):
    """Dialog for specifying sample grouping info"""
    
    def __init__(self, sample_names, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Specify Sample Grouping")
        self.setModal(True)
        self.resize(400, 300)
        
        self.sample_names = sample_names
        self.design_data = []
        
        self.setup_ui()
        self.populate_table()
        
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Info label
        info_label = QLabel("Specify which group each sample belongs to:")
        layout.addWidget(info_label)
        
        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Sample Name", "Group"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        layout.addWidget(self.table)
        
        # Buttons
        button_layout = QHBoxLayout()
        add_btn = QPushButton("Add Row")
        add_btn.clicked.connect(lambda: self.table.insertRow(self.table.rowCount()))
        remove_btn = QPushButton("Remove Row")
        remove_btn.clicked.connect(self.remove_row)
        button_layout.addWidget(add_btn)
        button_layout.addWidget(remove_btn)
        button_layout.addStretch()
        layout.addLayout(button_layout)
        
        # OK/Cancel buttons
        ok_cancel_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        ok_cancel_layout.addStretch()
        ok_cancel_layout.addWidget(ok_btn)
        ok_cancel_layout.addWidget(cancel_btn)
        layout.addLayout(ok_cancel_layout)
        
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def populate_table(self):
        self.table.setRowCount(len(self.sample_names))
        for i, sample_name in enumerate(self.sample_names):
            sample_item = QTableWidgetItem(sample_name)
            group_item = QTableWidgetItem("group1")
            
            # Make sample name non-editable
            sample_item.setFlags(sample_item.flags() & ~Qt.ItemIsEditable)
            
            self.table.setItem(i, 0, sample_item)
            self.table.setItem(i, 1, group_item)
            
    def remove_row(self):
        current_row = self.table.currentRow()
        if current_row >= 0:
            self.table.removeRow(current_row)
            
    def accept(self):
        self.design_data = []
        for i in range(self.table.rowCount()):
            sample_item = self.table.item(i, 0)
            group_item = self.table.item(i, 1)
            
            if sample_item and group_item:
                self.design_data.append({
                    'sample': sample_item.text(),
                    'group': group_item.text()
                })
                
        super().accept()


class NewProjectDialog(QDialog):
    """Dialog for creating a new transcriptome project"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("New Project")
        self.setModal(True)
        self.resize(600, 500)
        
        self.workspace_path = ""
        self.project_name = ""
        self.samples = []  # List of sample dictionaries
        self.index_file = ""
        self.annotation_file = ""  # GTF/GFF annotation file
        self.design_data = []  # List of sample-group mappings
        
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Workspace selection
        workspace_layout = QHBoxLayout()
        self.workspace_edit = QLineEdit()
        self.workspace_edit.setReadOnly(True)
        self.workspace_btn = QPushButton("Select Workspace")
        self.workspace_btn.clicked.connect(self.select_workspace)
        workspace_layout.addWidget(QLabel("Workspace Directory:"))
        workspace_layout.addWidget(self.workspace_edit)
        workspace_layout.addWidget(self.workspace_btn)
        layout.addLayout(workspace_layout)
        
        # Project name
        project_layout = QHBoxLayout()
        self.project_name_edit = QLineEdit()
        project_layout.addWidget(QLabel("Project Name:"))
        project_layout.addWidget(self.project_name_edit)
        layout.addLayout(project_layout)
        
        # Samples section
        samples_group = QGroupBox("Sample Files")
        samples_layout = QVBoxLayout()
        samples_group.setLayout(samples_layout)
        
        # Add sample button
        add_sample_btn = QPushButton("Add FASTQ Files")
        add_sample_btn.clicked.connect(self.add_sample)
        samples_layout.addWidget(add_sample_btn)
        
        # Samples table
        self.samples_table = QTableWidget()
        self.samples_table.setColumnCount(3)
        self.samples_table.setHorizontalHeaderLabels(["Sample Name", "R1 File", "R2 File"])
        self.samples_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        samples_layout.addWidget(self.samples_table)
        
        layout.addWidget(samples_group)
        
        # Index file
        index_layout = QHBoxLayout()
        self.index_edit = QLineEdit()
        self.index_edit.setReadOnly(True)
        self.index_btn = QPushButton("Select Index File")
        self.index_btn.clicked.connect(self.select_index_file)
        index_layout.addWidget(QLabel("HISAT2 Index:"))
        index_layout.addWidget(self.index_edit)
        index_layout.addWidget(self.index_btn)
        layout.addLayout(index_layout)
        
        # Annotation file
        annotation_layout = QHBoxLayout()
        self.annotation_edit = QLineEdit()
        self.annotation_edit.setReadOnly(True)
        self.annotation_btn = QPushButton("Select Annotation File")
        self.annotation_btn.clicked.connect(self.select_annotation_file)
        annotation_layout.addWidget(QLabel("Annotation File:"))
        annotation_layout.addWidget(self.annotation_edit)
        annotation_layout.addWidget(self.annotation_btn)
        layout.addLayout(annotation_layout)
        
        # Design info
        design_layout = QHBoxLayout()
        self.design_btn = QPushButton("Specify Grouping Info")
        self.design_btn.clicked.connect(self.specify_design)
        self.design_label = QLabel("Grouping info not specified")
        design_layout.addWidget(QLabel("Experimental Design:"))
        design_layout.addWidget(self.design_btn)
        design_layout.addWidget(self.design_label)
        layout.addLayout(design_layout)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.ok_btn = QPushButton("OK")
        self.ok_btn.clicked.connect(self.accept)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        button_layout.addStretch()
        button_layout.addWidget(self.ok_btn)
        button_layout.addWidget(self.cancel_btn)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
        
    def select_workspace(self):
        path = QFileDialog.getExistingDirectory(self, "Select Workspace Directory")
        if path:
            self.workspace_path = path
            self.workspace_edit.setText(path)
            
    def add_sample(self):
        dialog = AddFastqDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            sample_info = {
                'name': dialog.sample_name,
                'r1': dialog.r1_file,
                'r2': dialog.r2_file if dialog.r2_file else None
            }
            
            self.samples.append(sample_info)
            self.update_samples_table()
            
    def update_samples_table(self):
        self.samples_table.setRowCount(len(self.samples))
        for i, sample in enumerate(self.samples):
            name_item = QTableWidgetItem(sample['name'])
            r1_item = QTableWidgetItem(os.path.basename(sample['r1']))
            r2_item = QTableWidgetItem(os.path.basename(sample['r2']) if sample['r2'] else "None")
            
            # Make items non-editable
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            r1_item.setFlags(r1_item.flags() & ~Qt.ItemIsEditable)
            r2_item.setFlags(r2_item.flags() & ~Qt.ItemIsEditable)
            
            self.samples_table.setItem(i, 0, name_item)
            self.samples_table.setItem(i, 1, r1_item)
            self.samples_table.setItem(i, 2, r2_item)
            
    def select_index_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select HISAT2 Index File", "", 
                                              "HISAT2 Index Files (*.ht2 *.ht1);;All files (*)")
        if file:
            # Extract the index prefix by removing the .N.ht2 or .N.ht1 suffix
            index_prefix = self.extract_index_prefix(file)
            self.index_file = index_prefix
            self.index_edit.setText(index_prefix)
            
    def select_annotation_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Annotation File", "", 
                                              "Annotation Files (*.gtf *.gff *.gff3);;All files (*)")
        if file:
            self.annotation_file = file
            self.annotation_edit.setText(file)
            
    def extract_index_prefix(self, file_path):
        """
        Extract the index prefix from a HISAT2 index file path.
        HISAT2 index files are named like: prefix.1.ht2, prefix.2.ht2, etc.
        We need to remove the .N.ht2 or .N.ht1 suffix to get the prefix.
        """
        # Get the file name without directory
        file_name = os.path.basename(file_path)
        
        # Check if it's a HISAT2 index file (.N.ht2 or .N.ht1)
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
            
    def specify_design(self):
        if not self.samples:
            QMessageBox.warning(self, "Warning", "Please add sample files first")
            return
            
        sample_names = [s['name'] for s in self.samples]
        dialog = SampleDesignDialog(sample_names, self)
        if dialog.exec_() == QDialog.Accepted:
            self.design_data = dialog.design_data
            self.design_label.setText(f"Grouping info specified for {len(self.design_data)} samples")
            
    def accept(self):
        self.project_name = self.project_name_edit.text().strip()
        
        if not self.workspace_path:
            QMessageBox.warning(self, "Warning", "Please select workspace directory")
            return
            
        if not self.project_name:
            QMessageBox.warning(self, "Warning", "Please enter project name")
            return
            
        if not self.samples:
            QMessageBox.warning(self, "Warning", "Please add sample files")
            return
            
        if not self.index_file:
            QMessageBox.warning(self, "Warning", "Please select HISAT2 index file")
            return
            
        if not self.annotation_file:
            QMessageBox.warning(self, "Warning", "Please select annotation file")
            return
            
        if len(self.design_data) != len(self.samples):
            QMessageBox.warning(self, "Warning", "Please specify grouping info for all samples")
            return
            
        # Create project structure
        self.create_project_structure()
        super().accept()
        
    def create_project_structure(self):
        """Create project directory structure"""
        project_path = Path(self.workspace_path) / self.project_name
        project_path.mkdir(parents=True, exist_ok=True)
        
        # Create results directory
        results_path = project_path / f"{self.project_name}_results"
        results_path.mkdir(exist_ok=True)
        
        # Save project info to JSON
        project_info = {
            "name": self.project_name,
            "workspace": self.workspace_path,
            "samples": self.samples,
            "index_file": self.index_file,
            "annotation_file": self.annotation_file,
            "design_data": self.design_data,
            "created": pd.Timestamp.now().isoformat(),
            "status": "unprocessed"  # Initial status
        }
        
        with open(project_path / "project.json", "w", encoding="utf-8") as f:
            json.dump(project_info, f, indent=2, ensure_ascii=False)


class FilterWidget(QWidget):
    """Widget for filtering gene data"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_ui()
        
    def setup_ui(self):
        layout = QHBoxLayout()
        
        self.column_combo = QComboBox()
        self.column_combo.setMinimumWidth(150)
        self.operator_combo = QComboBox()
        self.operator_combo.addItems(["Greater than or equal to", "Greater than", "Equal to", "Less than or equal to", "Less than"])
        self.value_spinbox = QDoubleSpinBox()
        self.value_spinbox.setRange(-999999, 999999)
        self.value_spinbox.setDecimals(6)
        self.value_spinbox.setValue(0.0)
        
        layout.addWidget(QLabel("Column:"))
        layout.addWidget(self.column_combo)
        layout.addWidget(QLabel("Condition:"))
        layout.addWidget(self.operator_combo)
        layout.addWidget(QLabel("Value:"))
        layout.addWidget(self.value_spinbox)
        
        self.setLayout(layout)


class ProcessThread(QThread):
    """Thread for running the transcriptome analysis workflow"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)
    
    def __init__(self, project_data, project_path, plugin_path):
        super().__init__()
        self.project_data = project_data
        self.project_path = project_path
        self.plugin_path = plugin_path
        self.is_running = True
        self.checkpoint_file = os.path.join(project_path, "checkpoint.json")
        
    def run(self):
        try:
            # Load checkpoint if exists
            checkpoint = self.load_checkpoint()
            
            # Step 1: HISAT2 alignment
            if not checkpoint.get("hisat2_align", False):
                self.progress.emit("Starting HISAT2 alignment...")
                if not self.hisat2_align():
                    return
                # Update checkpoint
                checkpoint["hisat2_align"] = True
                self.save_checkpoint(checkpoint)
            
            # Step 2: FeatureCounts
            if not checkpoint.get("feature_counts", False):
                self.progress.emit("Starting FeatureCounts counting...")
                if not self.feature_counts():
                    return
                # Update checkpoint
                checkpoint["feature_counts"] = True
                self.save_checkpoint(checkpoint)
                
            # Step 3: PyDESeq2 analysis
            if not checkpoint.get("pydeseq2_analysis", False):
                self.progress.emit("Starting PyDESeq2 differential analysis...")
                if not self.pydeseq2_analysis():
                    return
                # Update checkpoint
                checkpoint["pydeseq2_analysis"] = True
                self.save_checkpoint(checkpoint)
                
            self.finished.emit("Analysis workflow completed")
        except Exception as e:
            self.error.emit(str(e))
            
    def load_checkpoint(self):
        """Load checkpoint from file"""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                self.console_output.emit(f"Failed to load checkpoint file: {str(e)}", "error")
        return {}
        
    def save_checkpoint(self, checkpoint):
        """Save checkpoint to file"""
        try:
            with open(self.checkpoint_file, 'w', encoding='utf-8') as f:
                json.dump(checkpoint, f, indent=2, ensure_ascii=False)
        except Exception as e:
            self.console_output.emit(f"Failed to save checkpoint file: {str(e)}", "error")
            
    def hisat2_align(self):
        """Run HISAT2 alignment for all samples"""
        try:
            self.console_output.emit("Starting HISAT2 alignment workflow...", "info")
            
            # Get HISAT2 path
            config_path = os.path.join(self.plugin_path, "config.json")
            hisat2_path = None
            samtools_path = None
            
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                for tool in config_data:
                    if tool.get("name").lower() == "hisat2":
                        hisat2_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                    elif tool.get("name").lower() == "samtools":
                        samtools_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
            
            if not hisat2_path or not os.path.exists(hisat2_path):
                self.error.emit("HISAT2 executable not found")
                return False
                
            if not samtools_path or not os.path.exists(samtools_path):
                self.error.emit("samtools executable not found")
                return False
                
            # Create output directory
            results_dir = os.path.join(self.project_path, f"{self.project_data['name']}_results")
            bam_dir = os.path.join(results_dir, "bam_files")
            os.makedirs(bam_dir, exist_ok=True)
            
            # Align each sample
            samples = self.project_data.get('samples', [])
            for i, sample in enumerate(samples):
                if not self.is_running:
                    return False
                    
                sample_name = sample['name']
                self.console_output.emit(f"Aligning sample {i+1}/{len(samples)}: {sample_name}", "info")
                
                # Build HISAT2 command
                cmd = [
                    hisat2_path,
                    "-x", self.project_data['index_file'],
                    "-p", "4"  # threads
                ]
                
                # Add reads files
                if sample.get('r2'):  # paired-end
                    cmd.extend(["-1", sample['r1'], "-2", sample['r2']])
                else:  # single-end
                    cmd.extend(["-U", sample['r1']])
                
                # Output to SAM file first
                sam_file = os.path.join(bam_dir, f"{sample_name}.sam")
                cmd.extend(["-S", sam_file])
                
                self.console_output.emit(" ".join(cmd), "command")
                
                # Run HISAT2
                result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
                
                if result.stdout:
                    for line in result.stdout.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "output")
                if result.stderr:
                    for line in result.stderr.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                            
                if result.returncode != 0:
                    self.error.emit(f"HISAT2 alignment failed for {sample_name}: {result.stderr}")
                    return False
                    
                # Convert SAM to BAM using samtools
                bam_file = os.path.join(bam_dir, f"{sample_name}.bam")
                sort_cmd = [
                    samtools_path,
                    "sort",
                    "-o", bam_file,
                    sam_file
                ]
                
                self.console_output.emit(" ".join(sort_cmd), "command")
                
                sort_result = subprocess.run(sort_cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
                
                if sort_result.stdout:
                    for line in sort_result.stdout.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "output")
                if sort_result.stderr:
                    for line in sort_result.stderr.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "error" if sort_result.returncode != 0 else "output")
                
                if sort_result.returncode != 0:
                    self.error.emit(f"samtools sort failed for {sample_name}: {sort_result.stderr}")
                    return False
                    
                # Remove temporary SAM file
                try:
                    os.remove(sam_file)
                except Exception as e:
                    self.console_output.emit(f"Warning: Failed to remove temporary SAM file {sam_file}: {str(e)}", "warning")
                    
            return True
        except Exception as e:
            self.error.emit(f"Error during HISAT2 alignment: {str(e)}")
            return False
            
    def feature_counts(self):
        """Run FeatureCounts on all BAM files"""
        try:
            # Get FeatureCounts path
            config_path = os.path.join(self.plugin_path, "config.json")
            featurecounts_path = None
            
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                for tool in config_data:
                    if tool.get("name").lower() == "featurecounts":
                        featurecounts_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
                        break
            
            if not featurecounts_path or not os.path.exists(featurecounts_path):
                self.error.emit("FeatureCounts executable not found")
                return False
                
            # Find annotation file (GTF)
            # For now, we'll assume it's in the same directory as the index or needs to be specified
            # In a real implementation, this would be part of the project configuration
            gtf_file = self.project_data["annotation_file"]  # This would need to be specified in the project configuration
            
            if not gtf_file or not os.path.exists(gtf_file):
                self.error.emit("Annotation file not found or invalid")
                return False
            
            # Create output directory
            results_dir = os.path.join(self.project_path, f"{self.project_data['name']}_results")
            bam_dir = os.path.join(results_dir, "bam_files")
            count_file = os.path.join(results_dir, "counts.txt")
            
            # Find all BAM files
            bam_files = []
            if os.path.exists(bam_dir):
                for file in os.listdir(bam_dir):
                    if file.endswith(".bam"):
                        bam_files.append(os.path.join(bam_dir, file))
            
            if not bam_files:
                self.error.emit("No BAM files found for counting")
                return False
                
            # Build FeatureCounts command
            cmd = [
                featurecounts_path,
                "-a", gtf_file,
                "-o", count_file,
                "-T", "4",  # threads
                "-p"  # paired-end
            ]
            cmd.extend(bam_files)
            
            self.console_output.emit(" ".join(cmd), "command")
            
            # Run FeatureCounts
            result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
            
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "output")
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line.strip():
                        self.console_output.emit(line.strip(), "error" if result.returncode != 0 else "output")
                        
            if result.returncode != 0:
                self.error.emit(f"FeatureCounts counting failed: {result.stderr}")
                return False
                
            return True
        except Exception as e:
            self.error.emit(f"Error during FeatureCounts counting: {str(e)}")
            return False
            
    def pydeseq2_analysis(self):
        """Run PyDESeq2 differential expression analysis"""
        try:
            self.progress.emit("Running differential expression analysis...")
            
            # Import PyDESeq2
            try:
                from pydeseq2.dds import DeseqDataSet
                from pydeseq2.ds import DeseqStats
            except ImportError:
                self.error.emit("PyDESeq2 library not installed, please install pydeseq2 library")
                return False
            
            # Load count data
            self.console_output.emit("Loading count data...", "info")
            results_dir = os.path.join(self.project_path, f"{self.project_data['name']}_results")
            count_file = os.path.join(results_dir, "counts.txt")
            
            if not os.path.exists(count_file):
                self.error.emit(f"Count file not found: {count_file}")
                return False
                
            # Load count data - handle FeatureCounts output format correctly
            self.console_output.emit("Parsing count file...", "info")
            try:
                # First, find the header line in the FeatureCounts output
                with open(count_file, 'r', encoding='utf-8', errors='ignore') as f:
                    lines = f.readlines()
                
                # Find the header line (starts with Geneid)
                header_idx = 0
                for i, line in enumerate(lines):
                    if line.startswith("Geneid"):
                        header_idx = i
                        break
                
                # Read data starting from header
                count_df = pd.read_csv(count_file, sep='\t', skiprows=header_idx)
                
                # Set Geneid as index
                count_df.set_index(count_df.columns[0], inplace=True)
                
                # Remove any extra metadata columns that FeatureCounts adds 
                # For FeatureCounts, the last N columns are the sample columns
                # where N is the number of samples
                sample_names = [s['name'] for s in self.project_data['samples']]
                # Select only the sample columns (last N columns)
                count_df = count_df.iloc[:, -len(sample_names):]
                # Rename columns to match sample names
                count_df.columns = sample_names
                
            except Exception as e:
                self.error.emit(f"Failed to parse count file: {str(e)}")
                return False
            
            # Prepare design DataFrame from project data
            self.console_output.emit("Building experimental design...", "info")
            design_data = self.project_data.get('design_data', [])
            if not design_data:
                self.error.emit("Experimental design data missing")
                return False
                
            # Create design DataFrame
            sample_names = [s['name'] for s in self.project_data['samples']]
            design_dict = {item['sample']: item['group'] for item in design_data}
            
            design_df = pd.DataFrame(list(design_dict.items()), columns=['sample', 'condition'])
            design_df.set_index('sample', inplace=True)
            
            # Filter samples that exist in both count and design
            common_samples = [s for s in count_df.columns if s in design_df.index]
            if not common_samples:
                self.error.emit("No common samples found between count file and experimental design")
                return False
                
            count_df = count_df[common_samples]
            design_df = design_df.loc[common_samples]
            
            self.console_output.emit("Running DESeq2 analysis...", "info")
            
            # Create DESeq2 object and run analysis
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
            deseq2_results_dir = os.path.join(results_dir, "deseq2_results")
            os.makedirs(deseq2_results_dir, exist_ok=True)
            
            if len(conditions) >= 2:
                # Perform all pairwise comparisons
                for i in range(len(conditions)):
                    for j in range(i + 1, len(conditions)):
                        cond1 = conditions[i]
                        cond2 = conditions[j]
                        
                        self.console_output.emit(f"Comparing {cond1} vs {cond2}...", "info")
                        
                        # Create contrast and get results
                        stat_res = DeseqStats(dds, contrast=("condition", cond1, cond2))
                        stat_res.summary()
                        results = stat_res.results_df
                        
                        # Reset index to make gene names a column
                        results = results.reset_index()
                        # Rename the index column to 'geneid'
                        if 'index' in results.columns:
                            results = results.rename(columns={'index': 'geneid'})
                        
                        # Save individual comparison results
                        comparison_name = f"{cond1}_vs_{cond2}"
                        comparison_file = os.path.join(deseq2_results_dir, f"{comparison_name}_results.txt")
                        results.to_csv(comparison_file, sep='\t', index=False)
                        
                # Create combined results for all comparisons
                self.console_output.emit("Extracting results...", "info")
                
                # Create a list to store all comparison results
                all_results = []
                
                # Perform all pairwise comparisons again for combined results
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
            
            # Save combined results
            deseq2_results = os.path.join(results_dir, "deseq2_results.txt")
            combined_results.to_csv(deseq2_results, sep='\t', index=False)
            
            # Filter out rows where baseMean = 0 and save to filtered file
            if 'baseMean' in combined_results.columns:
                filtered_results = combined_results[combined_results['baseMean'] != 0]
                deseq2_filtered_results = os.path.join(results_dir, "deseq2_results_filtered.txt")
                filtered_results.to_csv(deseq2_filtered_results, sep='\t', index=False)
                self.console_output.emit(f"Filtered out rows with baseMean=0, filtered results saved to: {deseq2_filtered_results}", "info")
            
            return True
        except Exception as e:
            import traceback
            self.error.emit(f"Error during PyDESeq2 analysis: {str(e)}\n{traceback.format_exc()}")
            return False


class TransHub(QWidget):
    """Main transcriptome data hub widget"""
    
    def __init__(self, config=None, plugin_path=None):
        super().__init__()
        self.config = config
        self.plugin_path = plugin_path or os.path.dirname(os.path.abspath(__file__))
        self.projects = []
        self.analyses = []
        self.current_data = None
        self.selected_project_row = -1
        self.process_thread = None
        
        self.init_ui()
        self.load_projects()
        
    def init_ui(self):
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Title
        title_label = QLabel("Transcriptome Data Hub")
        title_font = QFont()
        title_font.setPointSize(16)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)
        
        # Tab widget
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Projects tab
        self.projects_tab = QWidget()
        self.setup_projects_tab()
        self.tab_widget.addTab(self.projects_tab, "Projects")
        
        # Analyses tab
        self.analyses_tab = QWidget()
        self.setup_analyses_tab()
        self.tab_widget.addTab(self.analyses_tab, "Analysis")
        
        # Console tab
        self.console_tab = QWidget()
        self.setup_console_tab()
        self.tab_widget.addTab(self.console_tab, "Processing Log")
        
    def setup_projects_tab(self):
        layout = QVBoxLayout()
        self.projects_tab.setLayout(layout)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.new_project_btn = QPushButton("New Project")
        self.new_project_btn.clicked.connect(self.new_project)
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.clicked.connect(self.load_projects)
        
        # Workflow buttons
        self.process_btn = QPushButton("Start Processing")
        self.process_btn.clicked.connect(self.process_project)
        self.process_btn.setEnabled(False)
        
        self.stop_btn = QPushButton("Stop Processing")
        self.stop_btn.clicked.connect(self.stop_process)
        self.stop_btn.setEnabled(False)
        
        # Import to analysis button
        self.import_analysis_btn = QPushButton("Import to Analysis")
        self.import_analysis_btn.clicked.connect(self.import_to_analysis)
        self.import_analysis_btn.setEnabled(False)
        
        button_layout.addWidget(self.new_project_btn)
        button_layout.addWidget(self.refresh_btn)
        button_layout.addWidget(self.process_btn)
        button_layout.addWidget(self.stop_btn)
        button_layout.addWidget(self.import_analysis_btn)
        button_layout.addStretch()
        layout.addLayout(button_layout)
        
        # Projects table
        self.projects_table = QTableWidget()
        self.projects_table.setColumnCount(5)
        self.projects_table.setHorizontalHeaderLabels(["Project Name", "Workspace", "Created", "Status", "Progress"])
        self.projects_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.projects_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.projects_table.itemSelectionChanged.connect(self.on_project_selected)
        layout.addWidget(self.projects_table)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
    def setup_analyses_tab(self):
        layout = QVBoxLayout()
        self.analyses_tab.setLayout(layout)
        
        # Create splitter
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)
        
        # Left side - data table and filter controls
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_widget.setLayout(left_layout)
        
        # Data table
        self.data_table = QTableWidget()
        left_layout.addWidget(self.data_table)
        
        # Filter section
        filter_group = QGroupBox("Data Filtering")
        filter_layout = QVBoxLayout()
        filter_group.setLayout(filter_layout)
        
        self.filter_widget = FilterWidget()
        filter_layout.addWidget(self.filter_widget)
        
        filter_button_layout = QHBoxLayout()
        self.add_filter_btn = QPushButton("Add Condition")
        self.add_filter_btn.clicked.connect(self.add_filter_condition)
        self.clear_filter_btn = QPushButton("Clear Conditions")
        self.clear_filter_btn.clicked.connect(self.clear_filter_conditions)
        filter_button_layout.addWidget(self.add_filter_btn)
        filter_button_layout.addWidget(self.clear_filter_btn)
        filter_button_layout.addStretch()
        filter_layout.addLayout(filter_button_layout)
        
        left_layout.addWidget(filter_group)
        splitter.addWidget(left_widget)
        
        # Right side - filter conditions display
        right_widget = QWidget()
        right_layout = QVBoxLayout()
        right_widget.setLayout(right_layout)
        
        # Filter conditions label
        filter_conditions_label = QLabel("Filter Conditions:")
        filter_conditions_label.setStyleSheet("font-weight: bold;")
        right_layout.addWidget(filter_conditions_label)
        
        # Filter conditions container (similar to HISAT2 reads tags)
        self.filter_conditions_container = QFrame()
        self.filter_conditions_layout = QVBoxLayout()
        self.filter_conditions_container.setLayout(self.filter_conditions_layout)
        self.filter_conditions_container.setVisible(False)
        right_layout.addWidget(self.filter_conditions_container)
        
        right_layout.addStretch()
        splitter.addWidget(right_widget)
        
        # Set splitter sizes
        splitter.setSizes([700, 300])
        
    def add_filter_condition(self):
        """Add a filter condition"""
        if not hasattr(self, 'filter_conditions'):
            self.filter_conditions = []
            
        # Get filter values
        column = self.filter_widget.column_combo.currentText()
        operator = self.filter_widget.operator_combo.currentText()
        value = self.filter_widget.value_spinbox.value()
        
        if not column:
            QMessageBox.warning(self, "Warning", "Please select a column to filter")
            return
            
        # Create condition info
        condition_info = {
            'column': column,
            'operator': operator,
            'value': value
        }
        
        self.filter_conditions.append(condition_info)
        
        # Create tag widget (similar to HISAT2 reads tags)
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
        
        # Condition info label
        info_text = f"{column} {operator} {value}"
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
        close_btn.clicked.connect(lambda: self.remove_filter_condition(condition_info, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.filter_conditions_layout.addWidget(tag_widget)
        if not hasattr(self, 'filter_condition_tags'):
            self.filter_condition_tags = []
        self.filter_condition_tags.append((condition_info, tag_widget))
        
        # Show container
        self.filter_conditions_container.setVisible(True)
        
        # Apply filter
        self.apply_filter()
        
    def remove_filter_condition(self, condition_info, tag_widget):
        """Remove a filter condition"""
        if hasattr(self, 'filter_conditions') and condition_info in self.filter_conditions:
            self.filter_conditions.remove(condition_info)
        
        # Remove from tags list
        if hasattr(self, 'filter_condition_tags'):
            self.filter_condition_tags = [(ci, tw) for ci, tw in self.filter_condition_tags if ci != condition_info]
        
        # Remove widget
        self.filter_conditions_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Hide container if no conditions left
        if not hasattr(self, 'filter_conditions') or not self.filter_conditions:
            self.filter_conditions_container.setVisible(False)
            
        # Apply filter
        self.apply_filter()
        
    def clear_filter_conditions(self):
        """Clear all filter conditions"""
        if hasattr(self, 'filter_conditions'):
            self.filter_conditions.clear()
            
        if hasattr(self, 'filter_condition_tags'):
            # Remove all tag widgets
            for condition_info, tag_widget in self.filter_condition_tags:
                self.filter_conditions_layout.removeWidget(tag_widget)
                tag_widget.deleteLater()
            self.filter_condition_tags.clear()
            
        # Hide container
        self.filter_conditions_container.setVisible(False)
        
        # Apply filter (this will reset the table to show all data)
        self.apply_filter()
        
    def apply_filter(self):
        """Apply filter to the data table"""
        if self.current_data is None:
            return
            
        # Start with all data
        filtered_data = self.current_data.copy()
        
        # Apply each filter condition
        if hasattr(self, 'filter_conditions') and self.filter_conditions:
            for condition in self.filter_conditions:
                column = condition['column']
                operator = condition['operator']
                value = condition['value']
                
                # Apply filter based on operator
                if operator == "Greater than or equal to":
                    filtered_data = filtered_data[filtered_data[column] >= value]
                elif operator == "Greater than":
                    filtered_data = filtered_data[filtered_data[column] > value]
                elif operator == "Equal to":
                    filtered_data = filtered_data[filtered_data[column] == value]
                elif operator == "Less than or equal to":
                    filtered_data = filtered_data[filtered_data[column] <= value]
                elif operator == "Less than":
                    filtered_data = filtered_data[filtered_data[column] < value]
        
        # Update data table with filtered data
        self.update_data_table(filtered_data)
        
    def update_data_table(self, data):
        """Update data table with provided data"""
        self.data_table.setRowCount(len(data))
        self.data_table.setColumnCount(len(data.columns))
        self.data_table.setHorizontalHeaderLabels(data.columns.tolist())
        
        for i, row in data.iterrows():
            for j, (col, val) in enumerate(row.items()):
                item = QTableWidgetItem(str(val))
                if isinstance(val, (int, float)):
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.data_table.setItem(i, j, item)
                
        self.data_table.resizeColumnsToContents()
        
    def setup_console_tab(self):
        """Setup console tab similar to HISAT2 plugin"""
        layout = QVBoxLayout()
        self.console_tab.setLayout(layout)
        
        # Console text area
        self.console_text = QTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Courier New", 9))
        # Set dark background for console like HISAT2 plugin
        self.console_text.setStyleSheet("""
            QTextEdit {
                background-color: #272822;
                color: #f8f8f2;
            }
        """)
        layout.addWidget(self.console_text)
        
        # Control buttons
        controls_layout = QHBoxLayout()
        self.clear_console_btn = QPushButton("Clear")
        self.clear_console_btn.clicked.connect(self.clear_console)
        controls_layout.addWidget(self.clear_console_btn)
        controls_layout.addStretch()
        layout.addLayout(controls_layout)
        
    def clear_console(self):
        """Clear console text"""
        self.console_text.clear()
        
    def add_console_message(self, message, msg_type="info"):
        """Add message to console"""
        timestamp = pd.Timestamp.now().strftime("%H:%M:%S")
        if msg_type == "command":
            formatted = f"[{timestamp}] $ {message}"
        elif msg_type == "error":
            formatted = f"[{timestamp}] ERROR: {message}"
        elif msg_type == "warning":
            formatted = f"[{timestamp}] WARNING: {message}"
        else:
            formatted = f"[{timestamp}] {message}"
            
        self.console_text.append(formatted)
        self.console_text.moveCursor(QTextCursor.End)
        
    def new_project(self):
        dialog = NewProjectDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            self.load_projects()
            QMessageBox.information(self, "Success", "Project created successfully!")
            
    def load_projects(self):
        """Load existing projects"""
        self.projects_table.setRowCount(0)
        self.projects = []
        
        # Scan for project.json files in workspace directories
        if self.plugin_path:
            workspace_dirs = []
            # Try to find workspace directories - this is a simplified approach
            # In a real implementation, you might want to store known workspaces somewhere
            trans_dir = Path(self.plugin_path)
            if trans_dir.exists():
                # Look for directories that might contain projects
                for item in trans_dir.parent.parent.iterdir():
                    if item.is_dir() and item.name not in ['bin', 'src', '__pycache__']:
                        workspace_dirs.append(item)
            
            # For each potential workspace, look for projects
            for workspace_dir in workspace_dirs:
                if workspace_dir.exists():
                    for project_dir in workspace_dir.iterdir():
                        if project_dir.is_dir():
                            project_json = project_dir / "project.json"
                            if project_json.exists():
                                try:
                                    with open(project_json, 'r', encoding='utf-8') as f:
                                        project_data = json.load(f)
                                        self.projects.append({
                                            'path': str(project_dir),
                                            'data': project_data
                                        })
                                except Exception as e:
                                    self.add_console_message(f"Loading project failed {project_dir}: {str(e)}", "error")
            
            # Populate the table
            self.projects_table.setRowCount(len(self.projects))
            for i, project in enumerate(self.projects):
                project_data = project['data']
                name_item = QTableWidgetItem(project_data.get('name', 'Unknown'))
                workspace_item = QTableWidgetItem(project_data.get('workspace', 'Unknown'))
                created_item = QTableWidgetItem(project_data.get('created', 'Unknown'))
                
                # Determine status with checkpoint information
                status = project_data.get('status', 'Unknown')
                if status == 'processing':
                    # Check if there's a checkpoint file
                    checkpoint_file = Path(project['path']) / "checkpoint.json"
                    if checkpoint_file.exists():
                        status = 'paused'
                
                status_item = QTableWidgetItem(status)
                
                # Make items non-editable
                for item in [name_item, workspace_item, created_item, status_item]:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                
                self.projects_table.setItem(i, 0, name_item)
                self.projects_table.setItem(i, 1, workspace_item)
                self.projects_table.setItem(i, 2, created_item)
                self.projects_table.setItem(i, 3, status_item)
                self.projects_table.setItem(i, 4, QTableWidgetItem("0%"))
                
    def process_project(self):
        """Process the selected project through the workflow"""
        if self.selected_project_row == -1:
            return
            
        project = self.projects[self.selected_project_row]
        project_path = project['path']
        project_data = project['data']
        
        self.add_console_message(f"Starting project processing: {project_data['name']}")
        self.add_console_message(f"Project path: {project_path}")
        
        # Update UI for processing
        self.process_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Switch to console tab
        self.tab_widget.setCurrentIndex(2)
        
        # Update status to processing
        project_data['status'] = 'processing'
        self.projects_table.item(self.selected_project_row, 3).setText('processing')
        
        # Save updated project info
        try:
            project_json_path = Path(project_path) / "project.json"
            with open(project_json_path, 'w', encoding='utf-8') as f:
                json.dump(project_data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            self.add_console_message(f"Saving project status failed: {str(e)}", "error")
        
        # Start processing thread
        self.process_thread = ProcessThread(project_data, project_path, self.plugin_path)
        self.process_thread.progress.connect(self.on_progress)
        self.process_thread.finished.connect(self.on_process_finished)
        self.process_thread.error.connect(self.on_process_error)
        self.process_thread.console_output.connect(self.add_console_message)
        self.process_thread.start()
        
    def stop_process(self):
        """Stop the current processing"""
        if self.process_thread and self.process_thread.isRunning():
            self.process_thread.is_running = False
            self.process_thread.terminate()
            self.process_thread.wait()
            self.on_process_finished("Processing stopped")
            
    def on_progress(self, message):
        """Handle progress updates"""
        self.add_console_message(f"Progress: {message}", "info")
        self.progress_bar.setFormat(f"Processing: {message}")
        
    def on_process_finished(self, message):
        """Handle process completion"""
        self.process_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(f"Processing completed: {message}", "info")
        
        if self.selected_project_row >= 0:
            # Update status to completed
            project = self.projects[self.selected_project_row]
            project['data']['status'] = 'completed'
            self.projects_table.item(self.selected_project_row, 3).setText('completed')
            
            # Save updated project info
            try:
                project_json_path = Path(project['path']) / "project.json"
                with open(project_json_path, 'w', encoding='utf-8') as f:
                    json.dump(project['data'], f, indent=2, ensure_ascii=False)
                self.add_console_message("Project processing completed and status saved", "info")
            except Exception as e:
                self.add_console_message(f"Saving project status failed: {str(e)}", "error")
        
    def on_process_error(self, error_msg):
        """Handle process errors"""
        self.process_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(f"Processing failed: {error_msg}", "error")
        
        if self.selected_project_row >= 0:
            # Update status to error
            project = self.projects[self.selected_project_row]
            project['data']['status'] = 'error'
            self.projects_table.item(self.selected_project_row, 3).setText('error')
            
    def on_project_selected(self):
        """Handle project selection in the table"""
        selected_items = self.projects_table.selectedItems()
        if selected_items:
            self.selected_project_row = selected_items[0].row()
            # Enable process button if project is unprocessed or paused
            project = self.projects[self.selected_project_row]
            status = project['data'].get('status', 'unknown')
            self.process_btn.setEnabled(status in ['unprocessed', 'pending', 'paused'])
            # Enable import to analysis button if project is completed
            self.import_analysis_btn.setEnabled(status == 'completed')
        else:
            self.selected_project_row = -1
            self.process_btn.setEnabled(False)
            self.import_analysis_btn.setEnabled(False)
            
    def import_to_analysis(self):
        """Import selected project to analysis tab"""
        if self.selected_project_row == -1:
            return
            
        project = self.projects[self.selected_project_row]
        project_data = project['data']
        project_path = project['path']
        
        # Load analysis data
        if self.load_analysis_data(project_path):
            # Switch to analysis tab
            self.tab_widget.setCurrentIndex(1)  # Analysis tab
            QMessageBox.information(self, "Success", f"Project '{project_data['name']}' imported to analysis area")
            
    def load_analysis_data(self, project_path):
        """Load analysis data for visualization"""
        try:
            # Look for results files
            results_dir = os.path.join(project_path, f"{os.path.basename(project_path)}_results")
            
            # Try to load deseq2 results first - prefer filtered results
            deseq2_filtered_file = os.path.join(results_dir, "deseq2_results_filtered.txt")
            deseq2_file = os.path.join(results_dir, "deseq2_results.txt")
            count_file = os.path.join(results_dir, "counts.txt")
            
            data_file = None
            if os.path.exists(deseq2_filtered_file):
                data_file = deseq2_filtered_file
            elif os.path.exists(deseq2_file):
                data_file = deseq2_file
            elif os.path.exists(count_file):
                data_file = count_file
            else:
                # Look for any txt file in results directory
                for file in os.listdir(results_dir):
                    if file.endswith(".txt"):
                        data_file = os.path.join(results_dir, file)
                        break
            
            if not data_file or not os.path.exists(data_file):
                QMessageBox.warning(self, "Warning", "No analyzable data file found")
                return False
                
            # Load data
            df = pd.read_csv(data_file, sep='\t')
            self.current_data = df
            
            # Update data table
            self.update_data_table(df)
            
            # Update filter widget with column names
            self.filter_widget.column_combo.clear()
            self.filter_widget.column_combo.addItems(df.columns.tolist())
            
            # Clear any existing filter conditions
            self.clear_filter_conditions()
            
            return True
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load analysis data: {str(e)}")
            return False
        
    def export_data(self):
        """Export filtered data to CSV"""
        if self.current_data is not None:
            path, _ = QFileDialog.getSaveFileName(self, "导出数据", "", "CSV Files (*.csv)")
            if path:
                self.current_data.to_csv(path, index=False)
                QMessageBox.information(self, "成功", f"数据已导出到: {path}")
                
    def reset_filter(self):
        """Reset all filters"""
        # Implementation would reset filters and show all data
        pass


# Plugin entry point
class TransHubPlugin:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return TransHub(config=self.config, plugin_path=self.plugin_path)