# clusfps
clusfps is just a script that I wrote on the top of RDKit (2015_03_1) to do compound structure clustering.
This script is under the terms of the  WTFPL(Do What the Fuck You Want to Public License).


RDKit
---------
RDKit is a collection of cheminformatics and machine-learning software written in C++ and Python. Detail information can be found at [Github](https://github.com/rdkit/rdkit)

### RDKit on windows

1. cp the whole RDKit package to a folder like d:/RDKit_2015_03_1

2. Add environment variables RBASE=d:/RDKit_2015_03_1

3. Add environment variables  PYTHONPATH:%RBASE%

4. Add environment variables  PATH:%RDBASE%/lib

### RDKit on Linux

Detail installation guild can be found at [RDKit install guild](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md)

clusfps clustering process
---------
Two steps are included in clusfps: fingerprint generation and cluster calculation. Three type of fingerprint can be calculated:  Topological Fingerprints, MACCS Keys and Morgan Fingerprints. Morgan Fingerprints generation needs users to define radius through command-line. Cluster algorithm include Butina and Murtagh. When using Butina, a cutoff is needed , which means elements within this range of each other are considered to be neighbors. When using Murtagh, the number of clusters should be pre-defined. Example can be found at [Example](https://github.com/kaiwang0112006/clusfps/blob/master/example/example.txt)

<p align="center">
  <img src="https://github.com/kaiwang0112006/clusfps/blob/master/guild.png?raw=true" alt="clusfps Architecture"/>
</p>

