# Pocket-Definition
This scripts determine the KLIFS binding pocket numbering of the given Kinase. 

## Required Software:
 - Clustal Omega
 - EMBOSS
 - BLAST+

 **Note:** The following installation guide for *Mac* requires the package manager [Homebrew](https://brew.sh/). You can install [Homebrew](https://brew.sh/) using:
 ```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
 
 ### Clustal Omega installation:
 #### Linux
 ```bash
 sudo apt install clustalo
```
 #### Mac
  ```bash
# install argtable first
brew install argtable

# download clustal omega source code and unpack it
wget -qO- http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz > clustal-omega-1.2.1.tar.gz
tar zxvf clustal-omega-1.2.1.tar.gz

# install Clustal Omega
cd clustal-omega-1.2.1
./configure CFLAGS='-I/opt/homebrew/include' LDFLAGS='-L/opt/homebrew/lib'
make
make install
```

### EMBOSS
#### Linux
```bash
sudo apt-get install emboss
```
#### Mac 
  ```bash
brew install brewsci/bio/emboss
```

### Blast+
#### Linux
```bash
sudo apt-get install ncbi-blast+
```
#### Mac
  ```bash
brew install blast
```

## Usage
```bash
python3 pocket_definition.py [-h] -f FAMILY -n NAME -p PATH_TO_PDB_FILE
```
`FAMILY` should be replaced by the kinase's family (e.g., `PKA`), and `NAME` by the kinase's name (e.g., `Prkaca`).
