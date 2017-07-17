# spectra
# 
 
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
$ git clone https://github.com/steinbeck/spectra.git
```

## Compiling

Compiling the library is performed with Apache Maven and requires Java 1.7 or later:

```bash
spectra/$ mvn package
```
will create an all-in-one-jar under ./target
