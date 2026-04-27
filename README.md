# snowMapper ❄️🏔️🛰️🗺️ *(...Coming 30 April 2026!)*

[![Journal: The Cryosphere](https://img.shields.io/badge/Journal-The_Cryosphere-teal.svg)](https://doi.org/10.5194/tc-20-2209-2026)
[![EGU Preprint](https://img.shields.io/badge/EGUsphere-Preprint-white.svg)]([https://doi.org/10.5194/tc-20-2209-2026](https://egusphere.copernicus.org/preprints/2026/egusphere-2026-327/))
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.17663731-blue.svg)](https://doi.org/10.5281/zenodo.17663731)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### ⁉️ Will you attend EGU26?
### ✅ Come to [CR6.5](https://meetingorganizer.copernicus.org/EGU26/EGU26-347.html) to learn more about snowMapper!

## *A Python GEE package for mapping & reconstructing high-res daily snow cover.*

<img width="9376" height="6163" alt="snowMapper_schematic_github" src="https://github.com/user-attachments/assets/8d63cc20-d60f-4c73-b596-e64b8c8e864f" />

### ➡️ Fully-configurable architecture for optimised results
snowMapper is a **physics-informed machine learning (PIML) model** developed to generate daily reconstructions of snow cover across complex mountain terrain. The system integrates *in situ* measurements with gridded meteorological inputs and incorporates binary snow presence/absence observations derived from high-resolution satellite imagery. Its modular design enables users to tailor configurations to specific study sites, producing daily snow-cover maps at spatial resolutions typically ranging from 20 to 100 meters. 

### ➡️ Fully-configurable architecture for optimised results
The workflow includes a preprocessing pipeline compatible with imagery from Landsat 4–9 and Sentinel-2; multiple terrain- and land-cover–based masking options (i.e., forest, glaciers, surface water, elevation, urban areas); five fully-configurable schemes for converting satellite reflectance to binary snow cover (*more to come soon!*); and a quasi–physically based downscaling of coarse resolution meteorological inputs. 

### ➡️ Leveraging phyics-infomed machine learning for snow cover reconstruction 
Snow-cover reconstruction is accomplished through two sequential, configurable gap-filling procedures: an initial decision-tree step followed by a PIML classifier. The classifier can be trained either with local field observations or with *in situ* data originating from other regions, allowing the model to operate in a fully physics-informed mode in the absence of a local monitoring network. A built-in evaluation module compares model outputs with satellite-derived snow cover, providing accuracy metrics directly within the final product. 

### ➡️ Spatiotemporal aggregations, built for large-scale climatological analysis
Optional aggregation routines allow fractional snow-cover metrics to be generated across temporal and spatial scales. 

### ➡️ Cloud-based processing
The system operates entirely on Google Earth Engine via its Python API, reducing dependence on local data storage and eliminating local computational demands.

## Get started with snowMapper!
### 🌐 Choose one of two environment files based on your operating system.

For MacOs:
```
conda env create -f environment_macos.yaml
```
For Windows:
```
conda env create -f environment_windows.yaml
```
Activate environment:
```
conda activate snowmapper
```
```
import snowmapper as sm
```
