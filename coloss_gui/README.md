# COLOSS GUI

A modern graphical user interface for COLOSS (Coupled-channel calculations for nuclear reactions).

## Features

- **Intuitive Input Interface**: Easy-to-use forms for all COLOSS parameters
- **Real-time Output Monitoring**: Live display of calculation output
- **Integrated Plotting**: Automatic visualization of results
- **Modern Design**: Clean, professional interface with light and dark themes
- **File Management**: Save and load input files
- **Smart Executable Detection**: Automatically finds COLOSS executable

## Installation

### Prerequisites

- Python 3.8 or higher
- COLOSS compiled and ready (run `make` in the main COLOSS directory)

### Setup

1. Install required Python packages:

```bash
cd coloss_gui
pip install -r requirements.txt
```

Alternatively, install packages individually:
```bash
pip install PySide6 matplotlib numpy scipy
```

## Running the GUI

### From the command line:

```bash
cd coloss_gui
python main.py
```

Or make it executable and run directly:
```bash
chmod +x main.py
./main.py
```

### From the parent directory:

```bash
python coloss_gui/main.py
```

## Usage

### 1. Setting Up a Calculation

- **General Parameters**: Configure mesh points, angular steps, and Coulomb wave function settings
- **System Parameters**: Define projectile and target (Z, mass, name), angular momentum range, and energy
- **Optical Potential**: Set volume and surface potential parameters
- **Nonlocal Potential**: Enable/disable nonlocal potential corrections

### 2. Running Calculations

1. Adjust parameters in the input panel (left side)
2. Click "Generate Input File" or use File → Save to save your input
3. Click "Run" in the toolbar or use Run → Run COLOSS (Ctrl+R)
4. Monitor output in the "Output Log" tab
5. View results in the "Plot" tab when complete

### 3. Loading Example Files

Use File → Open to load existing input files from the `test/` directory:
- `6Li208Pb.in` - 6Li + 208Pb example
- `alpha28Si.in` - α + 28Si example
- And more...

### 4. Themes

Switch between light and dark themes via View menu:
- View → Light Theme
- View → Dark Theme

## File Structure

```
coloss_gui/
├── main.py              # Application entry point
├── main_window.py       # Main window with menus and layout
├── input_panel.py       # Input parameter forms
├── plot_widget.py       # Results visualization
├── log_widget.py        # Output log display
├── runner.py            # COLOSS execution manager
├── styles.py            # UI theming
├── path_utils.py        # Executable and path detection
├── requirements.txt     # Python dependencies
└── README.md           # This file
```

## Keyboard Shortcuts

- **Ctrl+N**: New file
- **Ctrl+O**: Open file
- **Ctrl+S**: Save file
- **Ctrl+R**: Run COLOSS
- **Ctrl+C**: Stop calculation
- **Ctrl+Q**: Quit application

## Troubleshooting

### COLOSS executable not found

If the GUI can't find the COLOSS executable:

1. Make sure you've compiled COLOSS:
   ```bash
   cd /path/to/COLOSS
   make
   ```

2. Set the `COLOSS_EXE` environment variable:
   ```bash
   export COLOSS_EXE=/path/to/COLOSS/COLOSS
   ```

### Missing Python packages

Install all required packages:
```bash
pip install -r requirements.txt
```

### GUI doesn't start

Make sure you have Python 3.8+ and PySide6:
```bash
python --version
pip show PySide6
```

## Development

The GUI is built with:
- **PySide6**: Qt6 bindings for Python
- **Matplotlib**: Plotting and visualization
- **NumPy**: Numerical data handling

## License

This GUI follows the same license as COLOSS.

## Support

For issues or questions about the GUI, please refer to the main COLOSS documentation and repository.
