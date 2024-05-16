## Scripts to quantify neuron reconstruction

Quantification between QDyeFinder and Imaris-auto trace using a simple **Intersection over Union (IoU) metric**

First, Imaris Filament traces are converted and exported to .swc using a custom-written plugin in Python.
Then the mat files and the swc are loaded into the `comparing_reconstructions.ipynb` script for generating Fig 12(b) - comparision with Imaris Autotracing and calculating overlap percentages.
