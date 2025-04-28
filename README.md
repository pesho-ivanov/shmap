# Map-SHmap - Fast and accurate sketch-based read mapper

![Test Status](https://github.com/pesho-ivanov/sweepmap/actions/workflows/test.yml/badge.svg)

## Description

Map-SHmap is an algorithm for sketch-based read mapping of genomic sequences.

## Usage

Requirements: C++2a

```
git clone --recurse-submodules git@github.com:pesho-ivanov/shmap.git
cd shmap
make
./release/shmap -h
```

## Dependencies

* [ankerl/unordered_dense](https://github.com/martinus/unordered_dense) -- fast hashmap
* [cmd_line_parser](https://github.com/jermp/cmd_line_parser) -- command line parser

## Evaluation
* [PBSIM](https://github.com/pfaucon/PBSIM-PacBio-Simulator) -- long read simulator
* [samtools](https://github.com/samtools/samtools) -- labeling reads
* [paftools](https://github.com/RBGKew/pypaftol/blob/master/paftools_tutorial.md) -- evaluating .paf
