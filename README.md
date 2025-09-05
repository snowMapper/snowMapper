# snowMapper ❄️🏔️🛰️🗺️ (...Coming soon!)

*A Python GEE package for mapping & reconstructing high-res daily snow cover.*

SnowMapper, is a modular, open-access, physics-informed, machine-learning-driven model for reconstructing snow cover. It is (a) trained on in-situ data, (b) forced by reanalysis-derived meteorological conditions, and (c) assimilated with binary snow cover from high-resolution satellite imagery from all Landsat missions (i.e. 1984 onwards) and Sentinel-2. The model output is daily snow cover, or monthly means (i.e. monthly snow cover fraction), at a scale of 100 m, and it is equipped with an integrated validation scheme providing the accuracy of the model run as metadata in the final output. The entire process is performed using the cloud-based resources of Google Earth Engine through the Python API. It is optimal for reconstructing monthly snow cover fraction across a mountain or catchment.

Choose one of two environment files based on your operating system.

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
