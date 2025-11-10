import os
import json
import pandas as pd
import re
import subprocess
import numpy as np
from pathlib import Path
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidget,
                             QTableWidgetItem, QTabWidget, QDialog, QFormLayout, QLineEdit,
                             QFileDialog, QMessageBox, QLabel, QComboBox, QGroupBox,
                             QHeaderView, QApplication, QCheckBox, QDoubleSpinBox, QFrame,
                             QTextEdit, QSpinBox, QProgressBar, QSplitter)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QTextCursor, QIcon

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


class CommitDialog(QDialog):
    """Dialog for entering commit message"""
    
    def __init__(self, parent=None, title="Commit to History", default_message=""):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(True)
        self.resize(500, 150)
        
        self.commit_message = ""
        self.branch_name = ""
        self.create_new_branch = False
        
        self.setup_ui(default_message)
        
    def setup_ui(self, default_message):
        layout = QVBoxLayout()
        
        # Commit message
        layout.addWidget(QLabel("Commit Message:"))
        self.message_edit = QTextEdit()
        self.message_edit.setMaximumHeight(100)
        self.message_edit.setPlainText(default_message)
        layout.addWidget(self.message_edit)
        
        # Branch name (only for new branch option)
        self.branch_layout = QHBoxLayout()
        self.branch_layout.addWidget(QLabel("New Branch Name:"))
        self.branch_edit = QLineEdit()
        self.branch_layout.addWidget(self.branch_edit)
        self.branch_layout_widget = QWidget()
        self.branch_layout_widget.setLayout(self.branch_layout)
        self.branch_layout_widget.setVisible(False)
        layout.addWidget(self.branch_layout_widget)
        
        # Buttons
        button_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        button_layout.addStretch()
        button_layout.addWidget(ok_btn)
        button_layout.addWidget(cancel_btn)
        layout.addLayout(button_layout)
        
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def set_create_new_branch(self, create_new):
        self.create_new_branch = create_new
        self.branch_layout_widget.setVisible(create_new)
        if create_new:
            self.setWindowTitle("Commit & Create New Branch")
        else:
            self.setWindowTitle("Commit to History")


class RemoteConfigDialog(QDialog):
    """Dialog for configuring remote repository and SSH key"""
    
    def __init__(self, parent=None, current_remote_url=""):
        super().__init__(parent)
        self.setWindowTitle("Configure Remote Repository")
        self.setModal(True)
        self.resize(600, 400)
        
        self.current_remote_url = current_remote_url
        self.ssh_key_generated = False
        
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Remote repository configuration
        remote_group = QGroupBox("Remote Repository")
        remote_layout = QVBoxLayout()
        remote_group.setLayout(remote_layout)
        
        remote_layout.addWidget(QLabel("Remote Repository URL:"))
        self.remote_url_edit = QLineEdit()
        self.remote_url_edit.setPlaceholderText("e.g., https://github.com/user/repo.git or git@github.com:user/repo.git")
        self.remote_url_edit.setText(self.current_remote_url)
        remote_layout.addWidget(self.remote_url_edit)
        
        # SSH Key generation section
        ssh_group = QGroupBox("SSH Key Setup")
        ssh_layout = QVBoxLayout()
        ssh_group.setLayout(ssh_layout)
        
        ssh_info = QLabel("Generate SSH key for secure communication with remote repository:")
        ssh_info.setWordWrap(True)
        ssh_layout.addWidget(ssh_info)
        
        self.generate_key_btn = QPushButton("Generate SSH Key Pair")
        self.generate_key_btn.clicked.connect(self.generate_ssh_key)
        ssh_layout.addWidget(self.generate_key_btn)
        
        self.key_display_area = QTextEdit()
        self.key_display_area.setReadOnly(True)
        self.key_display_area.setVisible(False)
        self.key_display_area.setMaximumHeight(150)
        ssh_layout.addWidget(self.key_display_area)
        
        ssh_note = QLabel("Note: After generating the key, you need to add the public key to your GitHub account.")
        ssh_note.setWordWrap(True)
        ssh_layout.addWidget(ssh_note)
        
        # Buttons
        button_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        button_layout.addStretch()
        button_layout.addWidget(ok_btn)
        button_layout.addWidget(cancel_btn)
        
        layout.addWidget(remote_group)
        layout.addWidget(ssh_group)
        layout.addLayout(button_layout)
        
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def generate_ssh_key(self):
        """Generate SSH key pair using ssh-keygen"""
        try:
            # Create temporary directory for SSH keys
            temp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp_keys")
            os.makedirs(temp_dir, exist_ok=True)
            
            # Generate SSH key pair
            key_path = os.path.join(temp_dir, "yrtools_temp_key")
            
            # Run ssh-keygen command
            cmd = ["ssh-keygen", "-t", "rsa", "-b", "4096", "-f", key_path, "-N", ""]
            result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='ignore')
            
            if result.returncode == 0:
                # Read public key
                pub_key_path = key_path + ".pub"
                if os.path.exists(pub_key_path):
                    with open(pub_key_path, 'r', encoding='utf-8') as f:
                        public_key = f.read()
                    
                    # Display public key
                    self.key_display_area.setVisible(True)
                    self.key_display_area.setPlainText(
                        f"SSH Key Generated Successfully!\n\n"
                        f"Public Key (add this to your GitHub account):\n{public_key}\n\n"
                        f"Private Key Location: {key_path}\n\n"
                        f"Instructions:\n"
                        f"1. Copy the public key above\n"
                        f"2. Go to GitHub Settings -> SSH and GPG keys\n"
                        f"3. Click 'New SSH key'\n"
                        f"4. Paste the public key and save"
                    )
                    self.ssh_key_generated = True
                    QMessageBox.information(self, "Success", "SSH key generated successfully! Please copy the public key to your GitHub account.")
                else:
                    QMessageBox.warning(self, "Warning", "SSH key generated but public key file not found.")
            else:
                QMessageBox.critical(self, "Error", f"Failed to generate SSH key:\n{result.stderr}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to generate SSH key: {str(e)}")


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
        
        # Initialize Git repository
        try:
            subprocess.run(["git", "init"], cwd=project_path, check=True, capture_output=True)
        except subprocess.CalledProcessError:
            pass  # Ignore Git errors for now
        
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
        # self.column_combo.setMinimumWidth(150)
        self.operator_combo = QComboBox()
        self.operator_combo.addItems(["≥", ">", "=", "≤", "<", "∈", "∉"])
        self.value_spinbox = QDoubleSpinBox()
        self.value_spinbox.setRange(-999999, 999999)
        self.value_spinbox.setDecimals(2)
        self.value_spinbox.setValue(0.0)
        
        # layout.addWidget(QLabel("Column:"))
        layout.addWidget(self.column_combo)
        # layout.addWidget(QLabel("Condition:"))
        layout.addWidget(self.operator_combo)
        # layout.addWidget(QLabel("Value:"))
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
                
            # Create gitignore and commit results
            self.create_gitignore_and_commit()
                
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
            self.console_output.emit("Starting FeatureCounts counting workflow...", "info")
            
            # Get FeatureCounts path
            config_path = os.path.join(self.plugin_path, "config.json")
            featurecounts_path = None
            
            if os.path.exists(config_path):
                with open(config_path, 'r', encoding='utf-8') as f:
                    config_data = json.load(f)
                for tool in config_data:
                    if tool.get("name").lower() == "featurecounts":
                        featurecounts_path = os.path.join(self.plugin_path, tool["path"].lstrip("/"))
            
            if not featurecounts_path or not os.path.exists(featurecounts_path):
                self.error.emit("FeatureCounts executable not found")
                return False
                
            # Get annotation file
            gtf_file = self.project_data.get('annotation_file')
            if not gtf_file or not os.path.exists(gtf_file):
                self.error.emit("Annotation file not found or invalid")
                return False
                
            # Find BAM files
            results_dir = os.path.join(self.project_path, f"{self.project_data['name']}_results")
            bam_dir = os.path.join(results_dir, "bam_files")
            
            bam_files = []
            if os.path.exists(bam_dir):
                for file in os.listdir(bam_dir):
                    if file.endswith(".bam"):
                        bam_files.append(os.path.join(bam_dir, file))
                        
            if not bam_files:
                self.error.emit("No BAM files found for counting")
                return False
                
            # Output file
            count_file = os.path.join(results_dir, "counts.txt")
            
            # Build FeatureCounts command
            cmd = [
                featurecounts_path,
                "-a", gtf_file,
                "-o", count_file,
                "-T", "4",  # threads
                "-t", "exon",
                "-g", "gene_id",
                "-p"  # paired-end
            ]
            
            # Add all BAM files
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
                
            # Calculate TPM and FPKM values
            self.console_output.emit("Calculating TPM and FPKM values...", "info")
            if not self._calculate_expression_values(results_dir, count_file, gtf_file):
                self.console_output.emit("Warning: Failed to calculate TPM/FPKM values", "warning")
                
            return True
        except Exception as e:
            self.error.emit(f"Error during FeatureCounts counting: {str(e)}")
            return False
            
    def _calculate_expression_values(self, results_dir, count_file, gtf_file):
        """Calculate TPM and FPKM expression values"""
        try:
            # Load count data
            self.console_output.emit("Loading count data...", "info")
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
            
            # Get gene lengths from GTF file
            self.console_output.emit("Extracting gene lengths from annotation file...", "info")
            gene_lengths = self._get_gene_lengths_from_gtf(gtf_file)
            
            # Get count data (sample columns only)
            count_columns = count_df.columns[1:-4]  # Sample columns (exclude Geneid and last 4 metadata columns)
            count_data = count_df[count_columns]
            gene_ids = count_df.iloc[:, 0]  # Gene IDs
            
            # Ensure count_data contains only numeric values
            count_data = count_data.apply(pd.to_numeric, errors='coerce').fillna(0)
            
            # Match gene lengths with count data
            lengths = []
            valid_genes = []
            valid_count_data = []
            
            for i, gene_id in enumerate(gene_ids):
                gene_id_str = str(gene_id)
                if gene_id_str in gene_lengths:
                    lengths.append(gene_lengths[gene_id_str])
                    valid_genes.append(gene_id_str)
                    valid_count_data.append(count_data.iloc[i])
            
            if not valid_genes:
                self.console_output.emit("No matching genes found between count data and annotation file", "error")
                return False
            
            lengths = np.array(lengths)
            valid_count_data = pd.DataFrame(valid_count_data, columns=count_columns)
            
            # Calculate TPM
            self.console_output.emit("Calculating TPM values...", "info")
            rpk = valid_count_data.div(lengths / 1000, axis=0)
            tpm = rpk.div(rpk.sum(axis=0) / 1e6, axis=1)
            tpm_df = tpm.copy()
            tpm_df.insert(0, "Geneid", valid_genes)

            # remove Chr / Start / End / Strand Columns
            tpm_df = tpm_df.drop(columns=['Chr', 'Start', 'End', 'Strand'])
            
            # Save TPM
            tpm_file = os.path.join(results_dir, "counts_tpm.txt")
            tpm_df.to_csv(tpm_file, sep='\t', index=False)
            
            # Save filtered TPM
            filtered_tpm_df = self._filter_zero_count_rows(tpm_df)
            filtered_tpm_file = os.path.join(results_dir, "counts_tpm_filtered.txt")
            filtered_tpm_df.to_csv(filtered_tpm_file, sep='\t', index=False)
            
            # Calculate FPKM
            self.console_output.emit("Calculating FPKM values...", "info")
            total_counts = count_data.sum(axis=0)  # Total counts per sample
            rpk = valid_count_data.div(lengths / 1000, axis=0)
            fpkm = rpk.div(total_counts / 1e6, axis=1)
            fpkm_df = fpkm.copy()
            fpkm_df.insert(0, "Geneid", valid_genes)

            # remove Chr / Start / End / Strand Columns
            fpkm_df = fpkm_df.drop(columns=['Chr', 'Start', 'End', 'Strand'])
            
            # Save FPKM
            fpkm_file = os.path.join(results_dir, "counts_fpkm.txt")
            fpkm_df.to_csv(fpkm_file, sep='\t', index=False)
            
            # Save filtered FPKM
            filtered_fpkm_df = self._filter_zero_count_rows(fpkm_df)
            filtered_fpkm_file = os.path.join(results_dir, "counts_fpkm_filtered.txt")
            filtered_fpkm_df.to_csv(filtered_fpkm_file, sep='\t', index=False)
            
            self.console_output.emit(f"TPM values saved to: {tpm_file}", "info")
            self.console_output.emit(f"Filtered TPM values saved to: {filtered_tpm_file}", "info")
            self.console_output.emit(f"FPKM values saved to: {fpkm_file}", "info")
            self.console_output.emit(f"Filtered FPKM values saved to: {filtered_fpkm_file}", "info")
            
            return True
        except Exception as e:
            self.console_output.emit(f"Failed to calculate expression values: {str(e)}", "error")
            return False
    
    def _get_gene_lengths_from_gtf(self, gtf_file):
        """Extract gene lengths from GTF annotation file"""
        gene_lengths = {}
        
        try:
            with open(gtf_file, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('#'):
                        continue  # Skip comment lines
                        
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue  # Not a valid GTF line
                        
                    feature_type = parts[2]
                    if feature_type != "exon":
                        continue  # Only process exon features for gene length calculation
                        
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
                        
                        # For gene_id, sum up lengths of all exons
                        if gene_id in gene_lengths:
                            gene_lengths[gene_id] += length
                        else:
                            gene_lengths[gene_id] = length
                    except ValueError:
                        # Skip lines with invalid coordinates
                        continue
                        
        except Exception as e:
            self.console_output.emit(f"Failed to parse annotation file: {str(e)}", "error")
            
        return gene_lengths
    
    def _filter_zero_count_rows(self, df):
        """Filter out rows where all count values are zero"""
        try:
            # Get count columns (exclude Geneid column)
            count_columns = [col for col in df.columns if col != 'Geneid']
            
            # Convert count data to numeric
            count_data = df[count_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
            
            # Identify rows where all counts are zero
            zero_rows = (count_data == 0).all(axis=1)
            
            # Filter out zero-count rows
            filtered_df = df[~zero_rows].copy()
            
            return filtered_df
            
        except Exception as e:
            self.console_output.emit(f"Failed to filter data: {str(e)}", "error")
            return df
            
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
            
    def create_gitignore_and_commit(self):
        """Create .gitignore file and commit results to Git"""
        try:
            results_dir = os.path.join(self.project_path, f"{self.project_data['name']}_results")
            
            # Create .gitignore file
            gitignore_content = """# Ignore everything
*

# But not directories (so empty dirs are tracked, but files in them are still ignored by *)
!*/

# But track key result files
!counts_tpm_filtered.txt
!counts_fpkm_filtered.txt
!deseq2_results_filtered.txt
"""
            
            gitignore_path = os.path.join(self.project_path, ".gitignore")
            with open(gitignore_path, 'w', encoding='utf-8') as f:
                f.write(gitignore_content)
                
            # Try to add and commit with Git
            try:
                # Add all files
                subprocess.run(["git", "add", "."], cwd=self.project_path, check=True, capture_output=True)
                
                # Commit
                subprocess.run(
                    ["git", "commit", "-m", "Initial commit: project configuration and analysis results"], 
                    cwd=self.project_path, 
                    check=True, 
                    capture_output=True
                )
            except subprocess.CalledProcessError:
                pass  # Ignore Git errors for now
        except Exception as e:
            self.console_output.emit(f"Failed to create .gitignore or commit: {str(e)}", "warning")

