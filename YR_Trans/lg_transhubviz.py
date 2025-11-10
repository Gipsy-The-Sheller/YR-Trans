from PyQt5.QtWidgets import (QWidget, QTabWidget, QVBoxLayout, QHBoxLayout, 
                             QPushButton, QFileDialog, QLabel, QComboBox, 
                             QCheckBox, QDoubleSpinBox, QLineEdit, QGroupBox, 
                             QFrame, QTextEdit, QMessageBox, QDialog)
from PyQt5.QtCore import Qt
import os
import subprocess
import pandas as pd
from .viz import VolcanoPlot

class VizWidget(QWidget):
    """
    Abstract visualization widget with Parameters and Results tabs.
    
    Args:
        data: Optional data matrix to use for visualization (default: None)
        project_path: Path to the project directory for git operations (default: None)
        plot_type: Type of plot to generate (default: 'volcano')
    """
    
    def __init__(self, data=None, project_path=None, plot_type='volcano'):
        super().__init__()
        self.data = data
        self.project_path = project_path
        self.plot_type = plot_type
        
        # Initialize UI components
        self.init_ui()
        
        # Connect signals
        self.connect_signals()
        
        # Set up default parameters based on plot type
        self.setup_default_parameters()
        
        # If data is provided, load it
        if data is not None:
            self.load_data(data)
            
    def init_ui(self):
        """Initialize the user interface."""
        # Main layout
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Tab widget for Parameters and Results
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Parameters tab
        self.parameters_tab = QWidget()
        self.tab_widget.addTab(self.parameters_tab, "Parameters")
        self.setup_parameters_tab()
        
        # Results tab
        self.results_tab = QWidget()
        self.tab_widget.addTab(self.results_tab, "Results")
        self.setup_results_tab()
        
    def setup_parameters_tab(self):
        """Set up the parameters tab."""
        layout = QVBoxLayout()
        self.parameters_tab.setLayout(layout)
        
        # Input file selection
        self.input_file_label = QLabel("Input Matrix:")
        layout.addWidget(self.input_file_label)
        
        self.input_file_lineedit = QLineEdit()
        layout.addWidget(self.input_file_lineedit)
        
        self.select_file_btn = QPushButton("Select File")
        self.select_file_btn.clicked.connect(self.select_input_file)
        layout.addWidget(self.select_file_btn)
        
        # Plot parameters
        self.plot_params_group = QGroupBox("Plot Parameters")
        self.plot_params_layout = QVBoxLayout()
        self.plot_params_group.setLayout(self.plot_params_layout)
        layout.addWidget(self.plot_params_group)
        
        # Design theme
        self.design_theme_label = QLabel("Design Theme:")
        self.design_theme_combo = QComboBox()
        self.design_theme_combo.addItems(['linedraw', 'classic', 'bw', 'dark', 'light', 'minimal', 'seaborn', 'gray', 'void', 'xkcd'])
        self.design_theme_combo.setCurrentIndex(0)  # Default to linedraw
        
        design_theme_layout = QHBoxLayout()
        design_theme_layout.addWidget(self.design_theme_label)
        design_theme_layout.addWidget(self.design_theme_combo)
        self.plot_params_layout.addLayout(design_theme_layout)
        
        # Color scheme
        self.color_scheme_label = QLabel("Color Scheme:")
        self.color_scheme_combo = QComboBox()
        self.color_scheme_combo.addItems(['discrete', 'gradient'])
        self.color_scheme_combo.setCurrentIndex(0)  # Default to discrete
        
        color_scheme_layout = QHBoxLayout()
        color_scheme_layout.addWidget(self.color_scheme_label)
        color_scheme_layout.addWidget(self.color_scheme_combo)
        self.plot_params_layout.addLayout(color_scheme_layout)
        
        # Additional parameters based on plot type
        self.setup_additional_parameters()
        
        # Plot button
        self.plot_btn = QPushButton("Plot")
        self.plot_btn.clicked.connect(self.generate_plot)
        layout.addWidget(self.plot_btn)
        
        # Hide file selection if data is provided
        if self.data is not None:
            self.input_file_lineedit.setVisible(False)
            self.select_file_btn.setVisible(False)
            self.input_file_label.setVisible(False)
            
    def setup_additional_parameters(self):
        """Set up additional parameters based on plot type."""
        pass
            
    def setup_results_tab(self):
        """Set up the results tab."""
        layout = QVBoxLayout()
        self.results_tab.setLayout(layout)
        
        # Results display area
        self.results_display = QTextEdit()
        self.results_display.setReadOnly(True)
        layout.addWidget(self.results_display)
        
        # Add a separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        layout.addWidget(separator)
        
        # Export button
        self.export_btn = QPushButton("Export Image")
        self.export_btn.clicked.connect(self.export_image)
        layout.addWidget(self.export_btn)
        
    def connect_signals(self):
        """Connect signals and slots."""
        pass
        
    def setup_default_parameters(self):
        """Set up default parameters based on plot type."""
        if self.plot_type == 'volcano':
            # Set default values for volcano plot
            self.design_theme_combo.setCurrentText('linedraw')
            self.color_scheme_combo.setCurrentText('discrete')
            self.x_threshold_spinbox.setValue(1.0)
            self.y_threshold_spinbox.setValue(0.05)
            self.threshold_indicator_checkbox.setChecked(True)
            self.under_threshold_color_checkbox.setChecked(True)
            self.color_distribution_convert_checkbox.setChecked(True)
            
    def select_input_file(self):
        """Open file dialog to select input file."""
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        file_dialog.setNameFilter("Text files (*.txt);;All files (*.*)")
        
        if file_dialog.exec_():
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                self.input_file_lineedit.setText(selected_files[0])
                
    def load_data(self, data):
        """Load data into the widget."""
        self.data = data
        
    def generate_plot(self):
        """Generate the plot based on current parameters."""
        if self.data is None:
            QMessageBox.warning(self, "No Data", "Please provide data or select an input file.")
            return
            
        # Create a temporary copy of the data for plotting
        plot_data = self.data.copy()
        
        # Perform any necessary data preprocessing
        # For volcano plot, we need log2FC and pvalue columns
        if self.plot_type == 'volcano':
            # Ensure required columns exist
            required_columns = ['log2FoldChange', 'pvalue']
            if not all(col in plot_data.columns for col in required_columns):
                QMessageBox.warning(self, "Missing Columns", 
                                  f"Data must contain columns: {required_columns}")
                return
                
            # Extract column indices
            log2fc_idx = plot_data.columns.get_loc('log2FoldChange')
            pvalue_idx = plot_data.columns.get_loc('pvalue')
            
            # Generate the plot using VolcanoPlot function
            try:
                # Prepare color scheme
                if self.color_scheme_combo.currentText() == 'discrete':
                    color_scheme = {
                        'up': self.color_scheme_discrete_up.text(),
                        'down': self.color_scheme_discrete_down.text(),
                        'no-DEGs': self.color_scheme_discrete_noDEGs.text()
                    }
                    color_scheme = ['discrete', color_scheme]
                else:
                    # Gradient color scheme
                    color_scheme = ['gradient', ['#f57f74', '#82cc5e']]
                    
                # Generate the plot
                plot = VolcanoPlot(
                    data=plot_data,
                    GeneID=0,
                    log2FC=log2fc_idx,
                    pvalue=pvalue_idx,
                    theme=self.design_theme_combo.currentText(),
                    color_scheme=color_scheme,
                    color_distribution_convert=self.color_distribution_convert_checkbox.isChecked(),
                    x_threshold=self.x_threshold_spinbox.value(),
                    y_threshold=self.y_threshold_spinbox.value(),
                    threshold_indicator=self.threshold_indicator_checkbox.isChecked(),
                    under_threshold_color=self.under_threshold_color_checkbox.isChecked(),
                    trimmode="none",
                    alt=True
                )
                
                # Display the plot in the results tab
                self.display_plot(plot)
                
                # Commit to git if project path is provided
                if self.project_path:
                    self.commit_to_git(plot_data)
                    
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to generate plot: {str(e)}")
                return
                
    def display_plot(self, plot):
        """Display the generated plot in the results tab."""
        # Convert plot to image and display
        # This would require additional code to convert plotnine plots to images
        # For now, just show a message
        self.results_display.setText(f"Plot generated successfully!\n\n{plot}")
        
    def commit_to_git(self, data):
        """Commit the current data to git with a new branch."""
        try:
            # Create a new branch named after the plot type
            branch_name = f"{self.plot_type}_plot"
            
            # Change to the project directory
            original_dir = os.getcwd()
            os.chdir(self.project_path)
            
            # Create a new branch
            subprocess.run(["git", "checkout", "-b", branch_name], check=True, capture_output=True)
            
            # Stage the data file
            # Assuming the data is saved in a file named after the plot type
            data_filename = f"{branch_name}.txt"
            data.to_csv(data_filename, sep='\t', index=False)
            subprocess.run(["git", "add", data_filename], check=True, capture_output=True)
            
            # Commit with a message describing the data
            commit_message = f"Generated {self.plot_type} plot with filtered data"
            subprocess.run(["git", "commit", "-m", commit_message], check=True, capture_output=True)
            
            # Switch back to original branch
            subprocess.run(["git", "checkout", "main"], check=True, capture_output=True)
            
            # Return to original directory
            os.chdir(original_dir)
            
            QMessageBox.information(self, "Success", f"Data committed to git branch '{branch_name}'")
            
        except subprocess.CalledProcessError as e:
            QMessageBox.critical(self, "Git Error", f"Failed to commit to git: {e}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to commit to git: {str(e)}")
            
    def export_image(self):
        """Export the generated plot as an image."""
        # Implementation would depend on how the plot is stored
        # For now, just show a message
        QMessageBox.information(self, "Export", "Image export functionality not implemented yet.")


class VolcanoViz(VizWidget):
    def __init__(self, data=None, project_path=None):
        super().__init__(data, project_path, plot_type='volcano')
        self.plot_type = 'volcano'

    def setup_additional_parameters(self):
        """Set up additional parameters based on plot type."""
        if self.plot_type == 'volcano':
            # Thresholds
            self.x_threshold_label = QLabel("X Threshold:")
            self.x_threshold_spinbox = QDoubleSpinBox()
            self.x_threshold_spinbox.setValue(1.0)
            self.x_threshold_spinbox.setDecimals(2)
            
            x_threshold_layout = QHBoxLayout()
            x_threshold_layout.addWidget(self.x_threshold_label)
            x_threshold_layout.addWidget(self.x_threshold_spinbox)
            self.plot_params_layout.addLayout(x_threshold_layout)
            
            self.y_threshold_label = QLabel("Y Threshold:")
            self.y_threshold_spinbox = QDoubleSpinBox()
            self.y_threshold_spinbox.setValue(0.05)
            self.y_threshold_spinbox.setDecimals(6)
            
            y_threshold_layout = QHBoxLayout()
            y_threshold_layout.addWidget(self.y_threshold_label)
            y_threshold_layout.addWidget(self.y_threshold_spinbox)
            self.plot_params_layout.addLayout(y_threshold_layout)
            
            # Checkboxes
            self.threshold_indicator_checkbox = QCheckBox("Threshold Indicator Line")
            self.threshold_indicator_checkbox.setChecked(True)
            
            self.under_threshold_color_checkbox = QCheckBox("Under Threshold Coloring")
            self.under_threshold_color_checkbox.setChecked(True)
            
            self.color_distribution_convert_checkbox = QCheckBox("Color Distribution Conversion")
            self.color_distribution_convert_checkbox.setChecked(True)
            
            self.plot_params_layout.addWidget(self.threshold_indicator_checkbox)
            self.plot_params_layout.addWidget(self.under_threshold_color_checkbox)
            self.plot_params_layout.addWidget(self.color_distribution_convert_checkbox)
            
            # Color scheme options
            self.color_scheme_discrete_label = QLabel("Discrete Colors:")
            self.color_scheme_discrete_up = QLineEdit("#FF0000")  # Red for up-regulated
            self.color_scheme_discrete_down = QLineEdit("#0000FF")  # Blue for down-regulated
            self.color_scheme_discrete_noDEGs = QLineEdit("#ADAAAB")  # Gray for no DEGs
            
            color_scheme_discrete_layout = QHBoxLayout()
            color_scheme_discrete_layout.addWidget(self.color_scheme_discrete_label)
            color_scheme_discrete_layout.addWidget(self.color_scheme_discrete_up)
            color_scheme_discrete_layout.addWidget(self.color_scheme_discrete_down)
            color_scheme_discrete_layout.addWidget(self.color_scheme_discrete_noDEGs)
            self.plot_params_layout.addLayout(color_scheme_discrete_layout)


class VolcanoDialog(QDialog):
    def __init__(self, parent=None, data=None, project_path=None):
        super().__init__(parent)
        self.volcano_widget = VolcanoViz(data, project_path)
        self.setWindowTitle("Volcano Plot")
        self.resize(800, 600)  # 设置合适的窗口大小
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.volcano_widget)
        self.setLayout(self.layout)