# PHScaffolding
Scaffolding method based on Pore-C reads and hypergraph approach
# Requirements
python 3.8.12,numpy 1.24.4,scipy 1.10.1,g++
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

Use Falign to align the Pore-C data with the contig file to generate a paf file.You need to ensure that you already have the aligned paf files and the original contig files.(The usage of Falign can be found at: https://github.com/xiaochuanle/Falign)

Run the code

`python main.py --paf path/map.paf --contig path/contig.fasta --output_dir path/output -q 50.0 -r 1.0 --min_w 20.0 --drop_w 0.3 --gap 500`

After sequentially entering the PAF file path, contig file path, and output directory path, the assembly results will be generated.

Use `python main.py -h` to view command line arguments.

# 5.Dataset

The dataset and test data can be found and downloaded at the following link: https://github.com/Suquana/PHScaffolding-Dataset​
