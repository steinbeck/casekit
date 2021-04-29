#        

# The Computer-Assisted-Structure-Elucidation Kit (CASEkit)

Copyright 2017 Christoph Steinbeck

License: MIT, see doc/mit.license

## Introduction

This project depends on the Chemistry Development Project (CDK), hosted under https://cdk.github.io/
Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is
available though our [Github site](http://cdk.github.io/cdk/).

## Download Source code

This assumes that you have git working on your system and you have initialised your local repository.

Then, downloading casekit is just a matter of

```bash
$ git clone https://github.com/michaelwenk/casekit
```

## Compiling

Compiling the library is performed with Apache Maven and requires Java 1.7 or later:

```bash
cd casekit
mvn clean package
```

will create an all-in-one-jar under ./target which you can use in your Java project.




