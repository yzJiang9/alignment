# Alignment

Protein sequences alignment using the Smith-Waterman algorithm with an affine gap penalty.

## Getting Started

Here is an example of how to install the package and use it.

### Installtion

```sh
pip install git+https://github.com/yzJiang9/alignment.git
```

### Usage

```py
from alignment import smith_waterman

inputFile = 'path/To/Input'
scoreFile = 'path/To/ScoreFile'
outputFile = 'path/To/Output'

smith_waterman.runSW(inputFile, scoreFile, outputFile, openGap = -2, extGap = -1)
```

`blosum62.txt` is provided and used as a substitution matrix in demos. Demos can be found under `examples`.

An output file consists of three parts as follows.
* Sequences: two sequences in the input file. **The input file should only contain two sequences, each line for one.**
* Scoring matrix: matrix calculated by the Smith-Waterman algorithm (by default, gap opening penealty is -2 and gap extension penalty is -1).

<p align="center">
  <img src="https://user-images.githubusercontent.com/96606552/161417120-3ea54cdc-4d65-4bf4-b922-75b0a9b80f97.png", alt="Smith-Waterman Algorithm" width='80% height='80%'/>
</p>
 
 * Best local alignments: two sequences with aligned results: `|` for a match and `-` for a gap.

## License

Distributed under the GPLv3 License. See `LICENSE` for more information.

## Contact

Yunzhe Jiang - yunzhe.jiang@yale.edu.

Course - [CB&B752: Biomedical Data Science: Mining and Modeling](http://cbb752b22.gersteinlab.org/) given by Prof. Mark Gerstein at Yale Univerisity.
