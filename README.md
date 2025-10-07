# 🔥 Flush End-Plate Beam-to-Column Joints – Parametric Study under Seismic and Fire Loading

This repository provides a **fully documented Abaqus CAE Python script** for simulating and performing **parametric finite element studies** on **Flush End-Plate (FEP) steel beam-to-column connections** subjected to **monotonic, cyclic (seismic), and fire loading**.  
The scripts form the computational foundation of an ongoing research paper titled:

> **“Numerical Investigation of Monotonic Seismic and Fire Performance of Flush End-Plate Beam-to-Column Joints.”**  
> *Manuscript under preparation (2025).*

---

## 🧭 Overview

- 🔁 Automated geometry generation, meshing, and partitioning  
- ⚙️ Parametric variation of plate thickness, bolt diameter, and joint spacing  
- 🔥 Coupled **thermal–mechanical simulation** capability for fire exposure  
- 🌊 Cyclic and monotonic loading amplitude control for seismic analysis  
- 📊 Automated extraction of load–displacement, moment–rotation, and temperature histories  
- 🧩 Modularized functions for direct reuse in other Abaqus studies  

---

## 📂 Repository Structure

```
FEP_Parametric_Study/
│
├── P1_FEP_ParametricStudy.py     # Main Abaqus Python script
├── LICENSE                       # MIT license
├── .gitignore                    # Ignore system / result files
├── README.md                     # Project documentation
├── /src/                         # Helper modules (if any)
├── /results/                     # Output ODB, CSV, and figures
└── /models/                      # Input/CAE models
```

---

## ⚙️ Usage

### Requirements
- **Abaqus CAE** (tested on versions 2021–2024)  
- Python **3.x** through the Abaqus scripting interface  
- Libraries: `os`, `numpy`, `csv`, `json`

### Running the Model
```bash
abaqus cae noGUI=P1_FEP_ParametricStudy.py
```

### Output Files
- `*.inp` and `*.cae` generated for each parametric case  
- Batch job submission for multiple runs  
- Summary file: `results/summary_results.csv` containing  
  - Ultimate load  
  - Rotation at failure  
  - Maximum temperature  
  - Connection damage state  

---

## 👨‍💻 Authors

**Md. Ibrahim Kholil¹***  
¹ Postgraduate Student, Department of Civil Engineering,  
Khulna University of Engineering & Technology (KUET), Bangladesh  
📧 engikholil@gmail.com · ORCID: [0000-0002-3349-8496](https://orcid.org/0000-0002-3349-8496)

**Muhammad Harunur Rashid²**  
² Professor, Department of Civil Engineering,  
Khulna University of Engineering & Technology (KUET), Bangladesh  
📧 mhrashid@ce.kuet.ac.bd · ORCID: [0000-0002-3656-5025](https://orcid.org/0000-0002-3656-5025)

**Aziz Ahmed³**  
³ Senior Lecturer, School of Civil, Mining, Environmental and Architectural Engineering,  
University of Wollongong, Australia  
📧 aziza@uow.edu.au · ORCID: [0000-0001-9707-2606](https://orcid.org/0000-0001-9707-2606)

---

## 🧾 License

Released under the **MIT License**.  
You are free to use, modify, and distribute the code with appropriate credit.

---

## 🧩 Related Publication (Under Processing)

This codebase supports the manuscript:

> **“Numerical Investigation of Monotonic Seismic and Fire Performance of Flush End-Plate Beam-to-Column Joints.”**  
> *Md. Ibrahim Kholil, Muhammad Harunur Rashid, and Aziz Ahmed*  

The developed Abaqus framework has been used to:
- Calibrate **moment–rotation** and **failure mechanisms** under cyclic/seismic loading.  
- Simulate **thermal degradation and residual capacity** under ISO-834 fire exposure.  
- Conduct **parametric studies** varying plate thickness, bolt strength, and joint geometry.

Once accepted and indexed, the publication DOI and full citation will be updated here.

---

## 🔍 DOI and Citation

This repository will be archived on **Zenodo** for citation and permanent DOI assignment.

**DOI (to be updated after Zenodo upload):**  
[https://doi.org/10.5281/zenodo.xxxxxxxx](https://doi.org/10.5281/zenodo.xxxxxxxx)

If you use this repository, please cite as:

```bibtex
@software{Kholil2025_FEP_ParametricStudy,
  author       = {Kholil, Md. Ibrahim and Rashid, Muhammad Harunur and Ahmed, Aziz},
  title        = {Flush End-Plate Beam-to-Column Joints – Parametric Study under Seismic and Fire Loading},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.xxxxxxxx},
  url          = {https://doi.org/10.5281/zenodo.xxxxxxxx}
}
```

---

## 🏗️ Repository Maintenance

**Maintainer:** Md. Ibrahim Kholil  
📧 engikholil@gmail.com  

Future updates will include:
- Additional parameter cases (bolt grade, edge distances, etc.)  
- Post-processing scripts for extracting moment–rotation envelopes  
- Fire duration vs. capacity degradation data  

---

*© 2025 Md. Ibrahim Kholil, Muhammad Harunur Rashid, and Aziz Ahmed.*
