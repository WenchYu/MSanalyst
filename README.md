# MSanalyst
This repository contains the original source code of MSanalyst. 
![MSanalystlogo](MSanalyst_logo.jpg)

# Installation
You can use the [online version](https://msanalyst.net/) or installing the stand-alone version from this repository
```bash
# Downloading msanalyst library using command or manually
wget --no-check-certificate 'https://drive.google.com/file/d/1w6HF3w1KIJlTz_QaVqqtN1BzkGDhDgzw/view?usp=sharing'

# Cloning MSanalyst repository
git clone git@github.com:WenchYu/MSanalyst.git && cd MSanalyst
unzip msdb.zip -d ./ && rm msdb.zip

# Creating environment
conda create -n msanalyst python=3.8
conda activate msanalyst
pip install -r requirements.txt
```


# Quick start
## Preprocess
Before applying MSanalyst, raw mass spectrometry (MS) data should be converted using [MSconvert](https://mzmine.github.io/mzmine_documentation/data_conversion.html). 
It is recommended to use [MZmine-based untargeted LC-MS workflow](https://mzmine.github.io/mzmine_documentation/workflows/lcmsworkflow/lcms-workflow.html) 
to generate the `quant.csv` and `mgf` files as inputs. 

## Module usage
Here we briefly introduce the command of MSanalyst quick start:
Using `-h` for help messages in MSanalyst: 

- `mn.py` 
for Default analysis workflow of MSanalyst.

```bash
python mn.py  -q ./example/example_quant.csv -m ./example/example.mgf -o ./example/
```

- `re-networking.py` for quick re-analysis of the results generated by `main.py` command.

```bash
python reanalysing.py -m ./example/example.mgf -q ./example/example_quant.csv -scm neutral_loss -scs 0.5 -scp 4
```

- `ms1search.py`
for single quick ms<sup>1</sup> searching

```bash
python ms1search.py -qms1 227.234
```

- `ms2search.py` 
for single ms<sup>2</sup> searching

```bash
python ms2search.py -m ./example/319.mgf
```

- `customizing.py` 
for generating customized mass library

```bash
python customizing.py -m ./example/exampleDB.mgf -li ./example/exampleDB.xlsx -o ./msdb/
python mn.py  -q ./example/example_quant.csv -m ./example/example.mgf -o ./example/ -e1f ./msdb/exampleDB_ms1.csv -e2f ./msdb/exampleDB_ms2.json
```

- `merging.py`
for merging different molecular networks

```bash
python mn.py  -q ./example/example_quant.csv -m ./example/example.mgf -o ./example/
python reanalysing.py -m ./example/example.mgf -q ./example/example_quant.csv -scm neutral_loss -scs 0.5 -scp 4
python merging.py -mn1 ./example/example_quant_result/example_modified_cosine_0.7_5.graphml -mn2 ./example/example_quant_result/example_neutral_loss_0.5_4.graphml -o ./example/
```

## Documentation
Please see the following links for detailed instructions, parameter usage and more information.
- [Online documentation](https://msanalyst.net/a/about) 
- [Video tutorials](https://msanalyst.net/a/about)
- [MSanalyst library](https://drive.google.com/file/d/1w6HF3w1KIJlTz_QaVqqtN1BzkGDhDgzw/view)

# 
MSanalyst was carefully tested, bugs may appear. Don't hesitate to contact us or open an issue. 
For other inquiries, please email us and we will be happy to answer you.
