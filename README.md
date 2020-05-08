# CHA₂DS₂-VASc Score Detection

This repository contains python scripts used in the study at The University of Utah Department of Cardiovascular Medicine to detect CHA₂DS₂-VASc scores in the EHR using natural language processing.

## Directory Structure

```
chadsvasc-score-detection
|
+-- src
|	|
|	+-- data
|	|
|	+-- results
|	|
|	+-- train.py
|
|
+-- README.md
```

## Method

### Usage

To extract CHA₂DS₂-VASc scores from the training metadata, run the following command on CLI:

```
python train.py
```