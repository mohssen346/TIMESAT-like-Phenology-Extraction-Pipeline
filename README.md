<div align="center">

# 🌿 TIMESAT-like Phenology Extraction Pipeline 🌿

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg?logo=python)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Build](https://img.shields.io/badge/status-active-success?logo=github)](#)
[![Made with](https://img.shields.io/badge/made%20with-💚%20Nature-green)](#)

_A modern Python pipeline for vegetation time series & phenology extraction_

</div>

---

## ✨ Key Highlights

> 🌱 Analyze vegetation time series from multi-band satellite imagery (e.g., Sentinel-2).  
> ⚡ Extract phenological metrics like SOS, EOS, Peak, Amplitude, etc.  
> 🌀 Advanced smoothing, function fitting, & noise handling.  
> 🔄 Detect multiple growing seasons.  
> 📊 Output CSV results + beautiful plots.  

---

## 🌟 Features

- 🌍 **Vegetation Indices** → NDVI, SAVI, MSAVI2, GNDVI, RVI  
- 🧹 **Time Series Smoothing** → Double Savitzky-Golay + outlier removal  
- 🎯 **Function Fitting** → Asymmetric Gaussian & Double Logistic models  
- 🌾 **Multiple Seasons Detection** → Peak finding with configurable thresholds  
- 📌 **Phenology Metrics** → SOS, EOS, Peak, LOS, Rates, Integrals, Harvest DOY  
- 🗺️ **Per-Region Processing** → Labeled raster analysis with filtering  
- 🖥️ **Parallel Processing** → Multiprocessing for large datasets  
- 🛠️ **Preprocessing Tools** → Deduplication + rasterization  
- 📑 **Outputs** → CSVs + PNG plots  
- ⚙️ **Configurable Parameters** → Central `config.py`  
- ✅ **Quality Control** → NaN handling, gap interpolation, validation  

---

## 🆚 Comparison with TIMESAT

| 🌿 Feature                        | ⏳ TIMESAT (Fortran)            | 🐍 This Python Code             |
|-----------------------------------|---------------------------------|---------------------------------|
| Double smoothing                  | ❌                              | ✅                              |
| Seasonal amplitude                | ✅                              | ✅                              |
| Asymmetric Gaussian fitting       | ✅                              | ✅                              |
| Double Logistic fitting           | ✅                              | ✅                              |
| Savitzky-Golay smoothing          | ✅                              | ✅                              |
| Multiple seasons detection        | ✅                              | ✅                              |
| Outlier removal (MAD/Z-score)     | ❌                              | ✅                              |
| Per-region processing             | ❌                              | ✅                              |
| Plotting & visualization          | ❌                              | ✅                              |
| Parallel processing               | ❌                              | ✅                              |
| Quality control                   | ❌                              | ✅                              |
| Irregular time-steps support      | ❌ (v4 in dev)                  | ❌                              |
| Input TIFF multi-band             | ❌                              | ✅                              |
| Output CSV/Plots                  | ❌                              | ✅                              |
| Fast processing large data        | ✅                              | ✅        |

> 💡 _This project keeps TIMESAT’s strengths while adding modern Python flexibility._

---

## ⚙️ Installation

### 📋 Prerequisites
- Python ≥ 3.8  
- Install dependencies:  
```bash
pip install -r requirements.txt
```

📦 Sample `requirements.txt`:
```
numpy
pandas
scipy
geopandas
rasterio
rioxarray
xarray
matplotlib
```

### 🚀 Setup
```bash
git clone https://github.com/yourusername/timesat-like-pipeline.git
cd timesat-like-pipeline
```

Then configure paths in `config.py`.  

---

## 🛠️ Usage

### 🔧 Preprocessing
```bash
python Preprocess.py
```
✨ Cleans duplicates + rasterizes polygons.  

### 🔬 Main Processing
```bash
python run.py
```
- Computes indices 🌿  
- Smooths + fits models 📈  
- Extracts phenology metrics 🌱  
- Saves CSVs + PNGs 💾  

**Outputs**:  
- `phenology_per_region_tile<TILE>.csv`  
- `phenology_ALL_TILES_COMPLETE.csv`  
- Time series CSVs + Plots per region  

---

## ⚙️ Configuration (`config.py`)

- 📂 **Paths** → Input, labels, outputs  
- 🌈 **Bands** → Sentinel-2 (e.g., B4=RED, B8=NIR)  
- 🌀 **Smoothing** → Savitzky-Golay params  
- 🌾 **Seasons** → Peak prominence, distance, amplitude  
- 📈 **Fitting** → Gaussian / Logistic  
- 🍂 **Phenology** → Amplitude thresholds, harvest ratio  
- 🖥️ **Processing** → Min region size, parallel CPUs  
- 💾 **Output** → Save CSVs & plots  


---

## 📜 License

```
MIT License

Copyright (c) 2025 Mohsen Forouzandeh

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```


---

## 🙏 Acknowledgments

- Inspired by TIMESAT (Jönsson & Eklundh, 2004)  🙏
- Powered by **NumPy, SciPy, GeoPandas, Rasterio** 💚  

📩 Questions? → [fmohssen161@gmail.com]  

---

👨‍💻  Developed by: Mohsen Forouzandeh  |  🏛️ University of Tehran

---
<div align="center">

✨ _“Bringing satellite data to life through open-source phenology tools.”_ ✨

</div>
