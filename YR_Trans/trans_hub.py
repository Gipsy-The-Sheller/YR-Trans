"""
TransHub: Transcriptome Data Management and Analysis Tool
Manages transcriptome projects and datasets with visualization capabilities
"""

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
from .lg_transhub import *

# config PATH for Git
plugin_path = os.path.dirname(os.path.abspath(__file__))
git_path = os.path.join(plugin_path, "bin/git/cmd")
os.environ["PATH"] = f"{git_path};{os.environ['PATH']}"
git_bash_path = os.path.join(plugin_path, "bin/git/bin")
os.environ["PATH"] = f"{git_bash_path};{os.environ['PATH']}"

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
        
        # Create tab widget for analysis types
        self.analysis_tab_widget = QTabWidget()
        layout.addWidget(self.analysis_tab_widget)
        
        # Expression Level Analysis tab
        self.expression_tab = QWidget()
        self.analysis_tab_widget.addTab(self.expression_tab, "Expression Level Analysis")
        self.setup_expression_tab()
        
        # Differential Analysis tab
        self.differential_tab = QWidget()
        self.analysis_tab_widget.addTab(self.differential_tab, "Differential Analysis")
        self.setup_differential_tab()
        
    def setup_expression_tab(self):
        """Setup the expression level analysis tab"""
        layout = QVBoxLayout()
        self.expression_tab.setLayout(layout)
        
        # Create splitter
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)
        
        # Left side - data table and filter controls
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_widget.setLayout(left_layout)
        
        # Data table
        self.expression_table = QTableWidget()
        left_layout.addWidget(self.expression_table)
        
        # Filter section
        filter_group = QGroupBox("Data Filtering")
        filter_layout = QVBoxLayout()
        filter_group.setLayout(filter_layout)
        
        self.expression_filter_widget = FilterWidget()
        filter_layout.addWidget(self.expression_filter_widget)
        
        filter_button_layout = QHBoxLayout()
        self.add_expression_filter_btn = QPushButton("Add Condition")
        self.add_expression_filter_btn.clicked.connect(self.add_expression_filter_condition)
        self.clear_expression_filter_btn = QPushButton("Clear Conditions")
        self.clear_expression_filter_btn.clicked.connect(self.clear_expression_filter_conditions)
        filter_button_layout.addWidget(self.add_expression_filter_btn)
        filter_button_layout.addWidget(self.clear_expression_filter_btn)
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
        self.expression_filter_conditions_container = QFrame()
        self.expression_filter_conditions_layout = QVBoxLayout()
        self.expression_filter_conditions_container.setLayout(self.expression_filter_conditions_layout)
        self.expression_filter_conditions_container.setVisible(False)
        right_layout.addWidget(self.expression_filter_conditions_container)
        
        # Git history controls
        git_group = QGroupBox("Version Control")
        git_layout = QVBoxLayout()
        git_group.setLayout(git_layout)
        
        self.commit_icon = QIcon(os.path.join(self.plugin_path, "YR_Trans", "src", "commit.svg"))
        self.branch_icon = QIcon(os.path.join(self.plugin_path, "YR_Trans", "src", "branch.svg"))
        self.git_icon    = QIcon(os.path.join(self.plugin_path, "YR_Trans", "src", "git.svg"))

        self.commit_btn = QPushButton("Commit to History")
        self.commit_btn.clicked.connect(self.commit_to_history)
        self.commit_btn.setIcon(self.commit_icon)
        self.commit_new_branch_btn = QPushButton("Commit && New Branch")
        self.commit_new_branch_btn.clicked.connect(self.commit_and_new_branch)
        self.commit_new_branch_btn.setIcon(self.branch_icon)
        self.manage_history_btn = QPushButton("Manage History")
        self.manage_history_btn.clicked.connect(self.manage_history)
        self.manage_history_btn.setIcon(self.git_icon)
        
        # New collaboration buttons
        self.rollback_btn = QPushButton("Scroll Back to Last History")
        self.rollback_btn.clicked.connect(self.rollback_to_last_commit)
        
        # Separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        
        self.configure_remote_btn = QPushButton("Configure Remote Repository")
        self.configure_remote_btn.clicked.connect(self.configure_remote_repository)
        
        self.sync_btn = QPushButton("Sync to Remote")
        self.sync_btn.clicked.connect(self.sync_with_remote)
        
        git_layout.addWidget(self.commit_btn)
        git_layout.addWidget(self.commit_new_branch_btn)
        git_layout.addWidget(self.manage_history_btn)
        git_layout.addWidget(self.rollback_btn)
        git_layout.addWidget(separator)
        git_layout.addWidget(self.configure_remote_btn)
        git_layout.addWidget(self.sync_btn)
        
        right_layout.addWidget(git_group)
        right_layout.addStretch()
        splitter.addWidget(right_widget)
        
        # Set splitter sizes
        splitter.setSizes([700, 300])
        
    def setup_differential_tab(self):
        """Setup the differential analysis tab"""
        layout = QVBoxLayout()
        self.differential_tab.setLayout(layout)
        
        # Create splitter
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)
        
        # Left side - data table and filter controls
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_widget.setLayout(left_layout)
        
        # Data table
        self.differential_table = QTableWidget()
        left_layout.addWidget(self.differential_table)
        
        # Filter section
        filter_group = QGroupBox("Data Filtering")
        filter_layout = QVBoxLayout()
        filter_group.setLayout(filter_layout)
        
        self.differential_filter_widget = FilterWidget()
        filter_layout.addWidget(self.differential_filter_widget)
        
        filter_button_layout = QHBoxLayout()
        self.add_differential_filter_btn = QPushButton("Add Condition")
        self.add_differential_filter_btn.clicked.connect(self.add_differential_filter_condition)
        self.clear_differential_filter_btn = QPushButton("Clear Conditions")
        self.clear_differential_filter_btn.clicked.connect(self.clear_differential_filter_conditions)
        filter_button_layout.addWidget(self.add_differential_filter_btn)
        filter_button_layout.addWidget(self.clear_differential_filter_btn)
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
        self.differential_filter_conditions_container = QFrame()
        self.differential_filter_conditions_layout = QVBoxLayout()
        self.differential_filter_conditions_container.setLayout(self.differential_filter_conditions_layout)
        self.differential_filter_conditions_container.setVisible(False)
        right_layout.addWidget(self.differential_filter_conditions_container)
        
        # Git history controls
        git_group = QGroupBox("Version Control")
        git_layout = QVBoxLayout()
        git_group.setLayout(git_layout)
        
        self.commit_btn2 = QPushButton("Commit to History")
        self.commit_btn2.clicked.connect(self.commit_to_history)
        self.commit_btn2.setIcon(self.commit_icon)
        self.commit_new_branch_btn2 = QPushButton("Commit && New Branch")
        self.commit_new_branch_btn2.clicked.connect(self.commit_and_new_branch)
        self.commit_new_branch_btn2.setIcon(self.branch_icon)
        self.manage_history_btn2 = QPushButton("Manage History")
        self.manage_history_btn2.clicked.connect(self.manage_history)
        self.manage_history_btn2.setIcon(self.git_icon)
        
        # New collaboration buttons
        self.rollback_btn2 = QPushButton("Scroll Back to Last History")
        self.rollback_btn2.clicked.connect(self.rollback_to_last_commit)
        
        # Separator line
        separator2 = QFrame()
        separator2.setFrameShape(QFrame.HLine)
        separator2.setFrameShadow(QFrame.Sunken)
        
        self.configure_remote_btn2 = QPushButton("Configure Remote Repository")
        self.configure_remote_btn2.clicked.connect(self.configure_remote_repository)
        
        self.sync_btn2 = QPushButton("Sync to Remote")
        self.sync_btn2.clicked.connect(self.sync_with_remote)
        
        git_layout.addWidget(self.commit_btn2)
        git_layout.addWidget(self.commit_new_branch_btn2)
        git_layout.addWidget(self.manage_history_btn2)
        git_layout.addWidget(self.rollback_btn2)
        git_layout.addWidget(separator2)
        git_layout.addWidget(self.configure_remote_btn2)
        git_layout.addWidget(self.sync_btn2)
        
        right_layout.addWidget(git_group)
        right_layout.addStretch()
        splitter.addWidget(right_widget)
        
        # Set splitter sizes
        splitter.setSizes([700, 300])
        
    def add_expression_filter_condition(self):
        """Add a filter condition for expression data"""
        if not hasattr(self, 'expression_filter_conditions'):
            self.expression_filter_conditions = []
            
        # Get filter values
        column = self.expression_filter_widget.column_combo.currentText()
        operator = self.expression_filter_widget.operator_combo.currentText()
        value = self.expression_filter_widget.value_spinbox.value()
        
        if not column:
            QMessageBox.warning(self, "Warning", "Please select a column to filter")
            return
            
        # Create condition info
        condition_info = {
            'column': column,
            'operator': operator,
            'value': value
        }
        
        self.expression_filter_conditions.append(condition_info)
        
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
        close_btn.clicked.connect(lambda: self.remove_expression_filter_condition(condition_info, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.expression_filter_conditions_layout.addWidget(tag_widget)
        if not hasattr(self, 'expression_filter_condition_tags'):
            self.expression_filter_condition_tags = []
        self.expression_filter_condition_tags.append((condition_info, tag_widget))
        
        # Show container
        self.expression_filter_conditions_container.setVisible(True)
        
        # Apply filter
        self.apply_expression_filter()
        
    def remove_expression_filter_condition(self, condition_info, tag_widget):
        """Remove an expression filter condition"""
        if hasattr(self, 'expression_filter_conditions') and condition_info in self.expression_filter_conditions:
            self.expression_filter_conditions.remove(condition_info)
        
        # Remove from tags list
        if hasattr(self, 'expression_filter_condition_tags'):
            self.expression_filter_condition_tags = [(ci, tw) for ci, tw in self.expression_filter_condition_tags if ci != condition_info]
        
        # Remove widget
        self.expression_filter_conditions_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Hide container if no conditions left
        if not hasattr(self, 'expression_filter_conditions') or not self.expression_filter_conditions:
            self.expression_filter_conditions_container.setVisible(False)
            
        # Apply filter
        self.apply_expression_filter()
        
    def clear_expression_filter_conditions(self):
        """Clear all expression filter conditions"""
        if hasattr(self, 'expression_filter_conditions'):
            self.expression_filter_conditions.clear()
            
        if hasattr(self, 'expression_filter_condition_tags'):
            # Remove all tag widgets
            for condition_info, tag_widget in self.expression_filter_condition_tags:
                self.expression_filter_conditions_layout.removeWidget(tag_widget)
                tag_widget.deleteLater()
            self.expression_filter_condition_tags.clear()
            
        # Hide container
        self.expression_filter_conditions_container.setVisible(False)
        
        # Apply filter (this will reset the table to show all data)
        self.apply_expression_filter()
        
    def apply_expression_filter(self):
        """Apply filter to the expression data table"""
        if not hasattr(self, 'current_expression_data') or self.current_expression_data is None:
            return
            
        # Start with all data
        filtered_data = self.current_expression_data.copy()
        
        # Apply each filter condition
        if hasattr(self, 'expression_filter_conditions') and self.expression_filter_conditions:
            for condition in self.expression_filter_conditions:
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
        self.update_expression_table(filtered_data)
        
    def add_differential_filter_condition(self):
        """Add a filter condition for differential data"""
        if not hasattr(self, 'differential_filter_conditions'):
            self.differential_filter_conditions = []
            
        # Get filter values
        column = self.differential_filter_widget.column_combo.currentText()
        operator = self.differential_filter_widget.operator_combo.currentText()
        value = self.differential_filter_widget.value_spinbox.value()
        
        if not column:
            QMessageBox.warning(self, "Warning", "Please select a column to filter")
            return
            
        # Create condition info
        condition_info = {
            'column': column,
            'operator': operator,
            'value': value
        }
        
        self.differential_filter_conditions.append(condition_info)
        
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
        close_btn.clicked.connect(lambda: self.remove_differential_filter_condition(condition_info, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.differential_filter_conditions_layout.addWidget(tag_widget)
        if not hasattr(self, 'differential_filter_condition_tags'):
            self.differential_filter_condition_tags = []
        self.differential_filter_condition_tags.append((condition_info, tag_widget))
        
        # Show container
        self.differential_filter_conditions_container.setVisible(True)
        
        # Apply filter
        self.apply_differential_filter()
        
    def remove_differential_filter_condition(self, condition_info, tag_widget):
        """Remove a differential filter condition"""
        if hasattr(self, 'differential_filter_conditions') and condition_info in self.differential_filter_conditions:
            self.differential_filter_conditions.remove(condition_info)
        
        # Remove from tags list
        if hasattr(self, 'differential_filter_condition_tags'):
            self.differential_filter_condition_tags = [(ci, tw) for ci, tw in self.differential_filter_condition_tags if ci != condition_info]
        
        # Remove widget
        self.differential_filter_conditions_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Hide container if no conditions left
        if not hasattr(self, 'differential_filter_conditions') or not self.differential_filter_conditions:
            self.differential_filter_conditions_container.setVisible(False)
            
        # Apply filter
        self.apply_differential_filter()
        
    def clear_differential_filter_conditions(self):
        """Clear all differential filter conditions"""
        if hasattr(self, 'differential_filter_conditions'):
            self.differential_filter_conditions.clear()
            
        if hasattr(self, 'differential_filter_condition_tags'):
            # Remove all tag widgets
            for condition_info, tag_widget in self.differential_filter_condition_tags:
                self.differential_filter_conditions_layout.removeWidget(tag_widget)
                tag_widget.deleteLater()
            self.differential_filter_condition_tags.clear()
            
        # Hide container
        self.differential_filter_conditions_container.setVisible(False)
        
        # Apply filter (this will reset the table to show all data)
        self.apply_differential_filter()
        
    def apply_differential_filter(self):
        """Apply filter to the differential data table"""
        if not hasattr(self, 'current_differential_data') or self.current_differential_data is None:
            return
            
        # Start with all data
        filtered_data = self.current_differential_data.copy()
        
        # Apply each filter condition
        if hasattr(self, 'differential_filter_conditions') and self.differential_filter_conditions:
            for condition in self.differential_filter_conditions:
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
        self.update_differential_table(filtered_data)
        
    def update_expression_table(self, data):
        """Update expression data table with provided data"""
        self.expression_table.setRowCount(len(data))
        self.expression_table.setColumnCount(len(data.columns))
        self.expression_table.setHorizontalHeaderLabels(data.columns.tolist())
        
        for i, row in data.iterrows():
            for j, (col, val) in enumerate(row.items()):
                # Format numeric values to 2 decimal places
                if isinstance(val, (int, float)) and not pd.isna(val):
                    formatted_val = f"{val:.2f}"
                    item = QTableWidgetItem(formatted_val)
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                else:
                    item = QTableWidgetItem(str(val))
                    if isinstance(val, (int, float)):
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.expression_table.setItem(i, j, item)
                
        self.expression_table.resizeColumnsToContents()
        
    def update_differential_table(self, data):
        """Update differential data table with provided data"""
        self.differential_table.setRowCount(len(data))
        self.differential_table.setColumnCount(len(data.columns))
        self.differential_table.setHorizontalHeaderLabels(data.columns.tolist())
        
        for i, row in data.iterrows():
            for j, (col, val) in enumerate(row.items()):
                # Format numeric values to 2 decimal places
                if isinstance(val, (int, float)) and not pd.isna(val):
                    formatted_val = f"{val:.2f}"
                    item = QTableWidgetItem(formatted_val)
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                else:
                    item = QTableWidgetItem(str(val))
                    if isinstance(val, (int, float)):
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.differential_table.setItem(i, j, item)
                
        self.differential_table.resizeColumnsToContents()
        
    def setup_console_tab(self):
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
        
    def new_project(self):
        """Create a new project"""
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
            self.add_console_message("Processing stopped by user", "info")
        
    def import_to_analysis(self):
        """Import selected project to analysis tab"""
        if self.selected_project_row == -1:
            return
            
        project = self.projects[self.selected_project_row]
        project_data = project['data']
        project_path = project['path']
        
        # Load analysis data
        if self.load_analysis_data(project_path):
            self.tab_widget.setCurrentIndex(1)
            QMessageBox.information(self, "Success", "Project data imported successfully!")
            
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
        
    def on_progress(self, message):
        """Handle progress updates"""
        self.add_console_message(message, "info")
        
    def on_process_finished(self, message):
        """Handle process completion"""
        self.process_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        self.add_console_message(f"Processing completed: {message}", "info")
        QMessageBox.information(self, "Success", "Processing completed successfully!")
        
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
            
    def load_analysis_data(self, project_path, exptype='tpm'):
        """Load analysis data for visualization"""
        try:
        # if True:
            # Look for results files
            results_dir = os.path.join(project_path, f"{os.path.basename(project_path)}_results")
            
            # Load expression level data (TPM/FPKM files)
            expression_files = []
            if os.path.exists(results_dir):
                for file in os.listdir(results_dir):
                    if ((file.endswith("_tpm.txt") or file.endswith("_fpkm.txt") or 
                        file.endswith("_tpm_filtered.txt") or file.endswith("_fpkm_filtered.txt"))):
                        expression_files.append(os.path.join(results_dir, file))
            

            # Load differential analysis data (DESeq2 results)
            deseq2_filtered_file = os.path.join(results_dir, "deseq2_results_filtered.txt")
            deseq2_file = os.path.join(results_dir, "deseq2_results.txt")
            
            # Load count data as fallback
            count_file = os.path.join(results_dir, "counts.txt")

            # Load expression data (prefer filtered versions)
            expression_data_file = None
            for file in expression_files:
                if exptype in file and "filtered" in file:
                    expression_data_file = file
                    break
                
            if not expression_data_file and os.path.exists(count_file):
                expression_data_file = count_file
                
            # Load differential data (prefer filtered DESeq2 results)
            differential_data_file = None
            if os.path.exists(deseq2_filtered_file):
                differential_data_file = deseq2_filtered_file
            elif os.path.exists(deseq2_file):
                differential_data_file = deseq2_file
            elif os.path.exists(count_file):
                differential_data_file = count_file
                
            # Load expression data if available
            if expression_data_file and os.path.exists(expression_data_file):
                try:
                # if True:
                    with open(count_file, 'r', encoding='utf-8', errors='ignore') as f:
                        lines = f.readlines()
            
                    # Find the header line (starts with Geneid)
                    header_idx = 0
                    for i, line in enumerate(lines):
                        if line.startswith("Geneid"):
                            header_idx = i
                            break
                    expression_df = pd.read_csv(expression_data_file, sep='\t', skiprows=header_idx)
                    self.current_expression_data = expression_df
                    self.update_expression_table(expression_df)
                    
                    # Update filter widget with column names
                    self.expression_filter_widget.column_combo.clear()
                    self.expression_filter_widget.column_combo.addItems(expression_df.columns.tolist())
                    
                    # Clear any existing filter conditions
                    self.clear_expression_filter_conditions()
                except Exception as e:
                    QMessageBox.warning(self, "Warning", f"Failed to load expression data: {str(e)}")
            else:
                # Clear expression table
                self.expression_table.setRowCount(0)
                self.expression_table.setColumnCount(0)
                self.current_expression_data = None
            
            # Load differential data if available
            if differential_data_file and os.path.exists(differential_data_file):
                try:
                    differential_df = pd.read_csv(differential_data_file, sep='\t')
                    self.current_differential_data = differential_df
                    self.update_differential_table(differential_df)
                    
                    # Update filter widget with column names
                    self.differential_filter_widget.column_combo.clear()
                    self.differential_filter_widget.column_combo.addItems(differential_df.columns.tolist())
                    
                    # Clear any existing filter conditions
                    self.clear_differential_filter_conditions()
                except Exception as e:
                    QMessageBox.warning(self, "Warning", f"Failed to load differential data: {str(e)}")
            else:
                # Clear differential table
                self.differential_table.setRowCount(0)
                self.differential_table.setColumnCount(0)
                self.current_differential_data = None
                
            # If neither data type is available, show warning
            if not (expression_data_file and os.path.exists(expression_data_file)) and \
               not (differential_data_file and os.path.exists(differential_data_file)):
                QMessageBox.warning(self, "Warning", "No analyzable data file found")
                return False
                
            return True
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load analysis data: {str(e)}")
            return False
        
    def export_data(self):
        """Export filtered data to CSV"""
        if self.current_data is not None:
            path, _ = QFileDialog.getSaveFileName(self, "Export Data", "", "CSV Files (*.csv)")
            if path:
                self.current_data.to_csv(path, index=False)
                QMessageBox.information(self, "Success", f"Data exported to: {path}")
                
    def reset_filter(self):
        """Reset all filters"""
        # Implementation would reset filters and show all data
        pass
        
    def commit_to_history(self):
        """Commit current state to Git history"""
        dialog = CommitDialog(self, "Commit to History", "Analysis results update")
        dialog.set_create_new_branch(False)
        if dialog.exec_() == QDialog.Accepted:
            commit_message = dialog.message_edit.toPlainText().strip()
            if not commit_message:
                QMessageBox.warning(self, "Warning", "Please enter a commit message")
                return
                
            # Get current project path
            if self.selected_project_row >= 0:
                project = self.projects[self.selected_project_row]
                project_path = project['path']
                
                try:
                    # Stage all changes
                    subprocess.run(["git", "add", "."], cwd=project_path, check=True, capture_output=True)
                    
                    # Commit
                    subprocess.run(
                        ["git", "commit", "-m", commit_message], 
                        cwd=project_path, 
                        check=True, 
                        capture_output=True
                    )
                    
                    QMessageBox.information(self, "Success", "Changes committed successfully!")
                except subprocess.CalledProcessError as e:
                    QMessageBox.critical(self, "Error", f"Failed to commit changes: {e}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
            else:
                QMessageBox.warning(self, "Warning", "No project selected")
                
    def commit_and_new_branch(self):
        """Commit current state and create a new branch"""
        dialog = CommitDialog(self, "Commit & Create New Branch", "Analysis results update")
        dialog.set_create_new_branch(True)
        if dialog.exec_() == QDialog.Accepted:
            commit_message = dialog.message_edit.toPlainText().strip()
            branch_name = dialog.branch_edit.text().strip()
            
            if not commit_message:
                QMessageBox.warning(self, "Warning", "Please enter a commit message")
                return
                
            if not branch_name:
                QMessageBox.warning(self, "Warning", "Please enter a branch name")
                return
                
            # Get current project path
            if self.selected_project_row >= 0:
                project = self.projects[self.selected_project_row]
                project_path = project['path']
                
                try:
                    # Stage all changes
                    subprocess.run(["git", "add", "."], cwd=project_path, check=True, capture_output=True)
                    
                    # Create and switch to new branch
                    subprocess.run(
                        ["git", "checkout", "-b", branch_name], 
                        cwd=project_path, 
                        check=True, 
                        capture_output=True
                    )
                    
                    # Commit
                    subprocess.run(
                        ["git", "commit", "-m", commit_message], 
                        cwd=project_path, 
                        check=True, 
                        capture_output=True
                    )
                    
                    QMessageBox.information(self, "Success", f"Changes committed and new branch '{branch_name}' created!")
                except subprocess.CalledProcessError as e:
                    QMessageBox.critical(self, "Error", f"Failed to commit and create branch: {e}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
            else:
                QMessageBox.warning(self, "Warning", "No project selected")
                
    def manage_history(self):
        """Open Git GUI to manage history"""
        # Get current project path
        if self.selected_project_row >= 0:
            project = self.projects[self.selected_project_row]
            project_path = project['path']

            # get the absolute path
            project_path_abs = os.path.abspath(project_path)
            
            try:
                # Try to run gitk
                cwd = os.getcwd()
                os.chdir(project_path_abs)
                # subprocess.Popen(["gitk"], cwd=project_path)
                subprocess.Popen("gitk .")
                os.chdir(cwd)
            except FileNotFoundError:
                QMessageBox.critical(self, "Error", "Git GUI (gitk) not found. Please ensure Git is properly installed.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to open Git GUI: {str(e)}")
        else:
            QMessageBox.warning(self, "Warning", "No project selected")
            
    def rollback_to_last_commit(self):
        """Rollback to the last commit"""
        if self.selected_project_row >= 0:
            project = self.projects[self.selected_project_row]
            project_path = project['path']
            
            # Confirm with user before rollback
            reply = QMessageBox.question(
                self, 
                "Confirm Rollback", 
                "This will discard all changes since the last commit. Are you sure you want to proceed?",
                QMessageBox.Yes | QMessageBox.No, 
                QMessageBox.No
            )
            
            if reply == QMessageBox.Yes:
                try:
                    # Perform git reset --hard to rollback to last commit
                    result = subprocess.run(
                        ["git", "reset", "--hard", "HEAD~1"], 
                        cwd=project_path, 
                        capture_output=True, 
                        text=True,
                        encoding='utf-8',
                        errors='ignore'
                    )
                    
                    if result.returncode == 0:
                        QMessageBox.information(self, "Success", "Rolled back to the previous commit successfully!")
                        
                        # Reload data and refresh tables
                        self.load_analysis_data(project_path)
                        
                        # Clear all filters
                        self.clear_expression_filter_conditions()
                        self.clear_differential_filter_conditions()
                    else:
                        QMessageBox.critical(self, "Error", f"Failed to rollback: {result.stderr}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"An error occurred during rollback: {str(e)}")
        else:
            QMessageBox.warning(self, "Warning", "No project selected")
            
    def configure_remote_repository(self):
        """Configure remote repository"""
        if self.selected_project_row >= 0:
            project = self.projects[self.selected_project_row]
            project_path = project['path']
            project_data = project['data']
            
            # Get current remote URL if exists
            current_remote_url = ""
            try:
                result = subprocess.run(
                    ["git", "remote", "get-url", "origin"], 
                    cwd=project_path, 
                    capture_output=True, 
                    text=True,
                    encoding='utf-8',
                    errors='ignore'
                )
                if result.returncode == 0:
                    current_remote_url = result.stdout.strip()
            except Exception:
                pass  # Ignore errors if no remote is set
            
            # Show dialog to configure remote
            dialog = RemoteConfigDialog(self, current_remote_url)
            if dialog.exec_() == QDialog.Accepted:
                remote_url = dialog.remote_url_edit.text().strip()
                
                if remote_url:
                    try:
                        # Remove existing remote if it exists
                        subprocess.run(
                            ["git", "remote", "remove", "origin"], 
                            cwd=project_path, 
                            capture_output=True
                        )
                    except Exception:
                        pass  # Ignore errors if remote doesn't exist
                    
                    # Add new remote
                    try:
                        result = subprocess.run(
                            ["git", "remote", "add", "origin", remote_url], 
                            cwd=project_path, 
                            capture_output=True, 
                            text=True,
                            encoding='utf-8',
                            errors='ignore'
                        )
                        
                        if result.returncode == 0:
                            # Save remote URL to project data
                            project_data['remote_url'] = remote_url
                            try:
                                project_json_path = Path(project_path) / "project.json"
                                with open(project_json_path, 'w', encoding='utf-8') as f:
                                    json.dump(project_data, f, indent=2, ensure_ascii=False)
                                QMessageBox.information(self, "Success", "Remote repository configured successfully!")
                            except Exception as e:
                                QMessageBox.warning(self, "Warning", f"Remote added to Git but failed to save to project file: {str(e)}")
                        else:
                            QMessageBox.critical(self, "Error", f"Failed to configure remote repository: {result.stderr}")
                    except Exception as e:
                        QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
                else:
                    # If no URL provided, just remove existing remote
                    try:
                        subprocess.run(
                            ["git", "remote", "remove", "origin"], 
                            cwd=project_path, 
                            capture_output=True
                        )
                        # Remove remote URL from project data
                        if 'remote_url' in project_data:
                            del project_data['remote_url']
                            try:
                                project_json_path = Path(project_path) / "project.json"
                                with open(project_json_path, 'w', encoding='utf-8') as f:
                                    json.dump(project_data, f, indent=2, ensure_ascii=False)
                            except Exception as e:
                                QMessageBox.warning(self, "Warning", f"Failed to update project file: {str(e)}")
                        QMessageBox.information(self, "Success", "Remote repository removed successfully!")
                    except Exception:
                        QMessageBox.information(self, "Info", "No remote repository configured.")
        else:
            QMessageBox.warning(self, "Warning", "No project selected")
            
    def sync_with_remote(self):
        """Sync with remote repository (pull & push)"""
        if self.selected_project_row >= 0:
            project = self.projects[self.selected_project_row]
            project_path = project['path']
            
            try:
                # Check if remote is configured
                result = subprocess.run(
                    ["git", "remote"], 
                    cwd=project_path, 
                    capture_output=True, 
                    text=True,
                    encoding='utf-8',
                    errors='ignore'
                )
                
                if result.returncode != 0 or not result.stdout.strip():
                    QMessageBox.warning(self, "Warning", "No remote repository configured. Please configure a remote repository first.")
                    return
                    
                # First, pull changes from remote
                self.add_console_message("Pulling changes from remote repository...", "info")
                pull_result = subprocess.run(
                    ["git", "pull"], 
                    cwd=project_path, 
                    capture_output=True, 
                    text=True,
                    encoding='utf-8',
                    errors='ignore'
                )
                
                if pull_result.returncode != 0:
                    # Handle pull errors (especially merge conflicts)
                    if "conflict" in pull_result.stderr.lower() or "conflict" in pull_result.stdout.lower():
                        QMessageBox.critical(self, "Merge Conflict", 
                                           "Merge conflict detected during pull. Please resolve conflicts manually using Git tools.")
                        self.add_console_message("Merge conflict detected during pull. Please resolve manually.", "error")
                        return
                    else:
                        QMessageBox.warning(self, "Pull Warning", f"Failed to pull changes: {pull_result.stderr}")
                        self.add_console_message(f"Failed to pull changes: {pull_result.stderr}", "warning")
                
                self.add_console_message("Pull completed successfully", "info")
                
                # Stage all changes
                subprocess.run(["git", "add", "."], cwd=project_path, capture_output=True)
                
                # Check if there are changes to commit
                status_result = subprocess.run(
                    ["git", "status", "--porcelain"], 
                    cwd=project_path, 
                    capture_output=True, 
                    text=True,
                    encoding='utf-8',
                    errors='ignore'
                )
                
                if status_result.stdout.strip():
                    # Commit changes with a default message
                    commit_result = subprocess.run(
                        ["git", "commit", "-m", "Auto-commit before sync"], 
                        cwd=project_path, 
                        capture_output=True, 
                        text=True,
                        encoding='utf-8',
                        errors='ignore'
                    )
                    
                    if commit_result.returncode != 0 and "nothing to commit" not in commit_result.stderr:
                        QMessageBox.warning(self, "Commit Warning", f"Failed to commit changes: {commit_result.stderr}")
                        self.add_console_message(f"Failed to commit changes: {commit_result.stderr}", "warning")
                
                # Push changes to remote
                self.add_console_message("Pushing changes to remote repository...", "info")
                push_result = subprocess.run(
                    ["git", "push", "origin", "HEAD"], 
                    cwd=project_path, 
                    capture_output=True, 
                    text=True,
                    encoding='utf-8',
                    errors='ignore'
                )
                
                if push_result.returncode == 0:
                    QMessageBox.information(self, "Success", "Sync with remote repository completed successfully!")
                    self.add_console_message("Push completed successfully", "info")
                else:
                    # Handle common push errors
                    error_output = push_result.stderr.lower()
                    if "permission denied" in error_output or "authentication failed" in error_output:
                        QMessageBox.critical(self, "Authentication Error", 
                                           "Authentication failed. Please check your SSH keys or credentials.")
                    elif "rejected" in error_output:
                        QMessageBox.warning(self, "Push Rejected", 
                                          "Push was rejected. You may need to pull and merge changes first.")
                    else:
                        QMessageBox.critical(self, "Push Error", f"Failed to push changes: {push_result.stderr}")
                    self.add_console_message(f"Failed to push changes: {push_result.stderr}", "error")
                    
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred during sync: {str(e)}")
                self.add_console_message(f"Sync error: {str(e)}", "error")
        else:
            QMessageBox.warning(self, "Warning", "No project selected")


# Plugin entry point
class TransHubPlugin:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
        
    def run(self):
        return TransHub(config=self.config, plugin_path=self.plugin_path)