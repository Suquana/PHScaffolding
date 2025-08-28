# PHScaffolding
Scaffolding method based on Pore-C reads and hypergraph approach
# Requirements
python 3.8.12,numpy 1.24.4,scipy 1.10.1,falign,g++
# Installation
## 1.Create a virtual environment
`conda create -n PHScaffolding python=3.8`

`conda activate PHScaffolding`
## 2.Install required packages
`conda install numpy,scipy`
## 3.Clone the repository
`git clone https://github.com/Suquana/PHScaffolding`

`cd PHScaffolding`
## 4.run
`cd code`

Modify the contig names in the contig file to facilitate assembly (requires modifying the file paths inside).

`python gaiming.py`

Use Falign to align the Pore-C data with the contig file to generate a `.map` file.

`g++ smap.cpp -o samp`

Generate a contig length file, with the input being a contig file.

`python line.py`

Process the alignment file.

`./smap path/xxx.map path/smap.txt`

Modify the paths in tiqu.py to generate an extraction file, with the input file being a `.map` file.

`python tiqu.py path/xxx.map path/tiqu.txt`

Generate a hypergraph file by modifying the file paths within it. Set the input path to samp.txt, and output a hypergraph file. 

`pyhton chap2.py`

Hypergraph clustering partition, modify the file paths inside, requires a hypergraph file as input, and the juzhen output can be ignored.

`python chaobian5.py`

Perform directional sorting, which requires modifying the file paths in the code. The inputs needed are a contig length file, a hypergraph clustering file, and an `tiqu.txt` file.

`python binpai4.py`

Connect to form scaffolds.

`python lianjie1.py path/binpai4.out scaffolds.fasta`
