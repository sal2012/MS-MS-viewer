 MS/MS Viewer - Final Project

This repository hosts the `Final_project_presentation.py` Python script, developed as part of a class project aimed at analyzing mass spectrometry (MS/MS) data to determine the match quality of peptides to spectra.

## Project Overview

The MS/MS Viewer is designed to take mzXML files along with a scan number and peptide sequence as input, compute the peptide's b-ion and y-ion m/z values, and match these with peaks in the data. The output is an annotated figure/plot that visually represents these matches.

### Key Features

- **Input Handling**: Accepts mzXML files, scan numbers, and peptide sequences through the command line.
- **Ion Computation**: Calculates b-ion and y-ion m/z values from peptide sequences.
- **Peak Matching**: Identifies peaks that match the computed m/z values and annotates them.
- **Visualization**: Generates a plot showing the matched peaks with labels.

### Modules Used

- System Module
- XML Parsing Module
- Array Handling Module
- Base64 Encoding Module
- Data Visualization Module

## Getting Started

### Prerequisites

Ensure you have Python and necessary libraries installed. The project uses several Python libraries which you can install using pip:

```bash
pip install numpy pandas matplotlib biopython

## Challenges and Limitations

- **Decimal Precision**: The script does not use exact decimal numbers for monoisotopic masses of amino acids, which may affect accuracy.
- **Ion Production**: Not all possible ions produced from peptide fragmentation are accounted for in the calculations.
- **Empty Lists**: Occasionally, the matched list for y-ions may be empty due to bugs in the script.

## Contact

For more information or queries regarding this project, please contact:
- **Salwa Elsaadawy**
- Email: [sme76@georgetown.edu]
