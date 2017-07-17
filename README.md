# 
# The Computer-Assisted-Structure-Elucidation Kit (CASEkit)
 
Copyright 2017 Christoph Steinbeck

License: MIT, see doc/mit.license

## Introduction

This project hosts various Java classes for teaching and research dealing with spectral data in chemistry and metabolomics.
This project depends on the Chemistry Development Project (CDK), hosted under http://cdk.github.io/
Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is available though our [Github site](http://cdk.github.io/cdk/).

## Download Spectra Source code

This assumes that you have git working on your system and you have initialised your local repository. Refer to https://help.github.com/articles/set-up-git/ for more

Then, downloading spectra is just a matter of

```bash
$ git clone https://github.com/steinbeck/casedk.git
```

## Compiling

Compiling the library is performed with Apache Maven and requires Java 1.7 or later:

```bash
spectra/$ mvn package
```
will create an all-in-one-jar under ./target

## Usage

### Shift Prediction with HOSE codes

The following classes are to demonstrate the prediction of Carbon-13 NMR spectra with HOSE codes. They only demonstrate the basic working principle and ignore, for example, stereochemistry, which can lead to large errors in, for example, the prediction of shifts for E/Z configurations of double bonds. If you want to know more about the details and a sophisticated system implementing them, please refer to Schutz V, Purtuc V, Felsinger S, Robien W (1997) CSEARCH-STEREO: A new generation of NMR database systems allowing three-dimensional spectrum prediction. Fresenius Journal of Analytical Chemistry 359:33–41. doi: 10.1007/s002160050531.

#### NMRShiftDBSDFParser

Take the NMRShiftDB SDF with assigned spectra (download from help section of NMRShiftDB.org) and produces a Tab-separated file with HOSE codes and assigned shift values. This file can then be read by HOSECodePredictor and SimilarityRanker. 

```bash
usage: java -jar spectra.jar casekit.NMRShiftDBSDFParser -i <arg> -o <arg> [-v] '[-d <arg>]' [-m <arg>]
Generates a table of HOSE codes and assigned shifts from an NMRShiftDB SDF
file from http://nmrshiftdb.nmr.uni-koeln.de/portal/js_pane/P-Help.

 -i,--infile <arg>       filename of NMRShiftDB SDF with spectra
                         (required)
 -o,--outfile <arg>      filename of generated HOSE code table (required)
 -v,--verbose            generate messages about progress of operation
 -d,--picdir <arg>       store pictures in given directory
 -m,--maxspheres <arg>   maximum sphere size up to which to generate HOSE
                         codes
```
