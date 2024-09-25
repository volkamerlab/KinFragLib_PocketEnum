# Required Software:
 - Clustal Omega
 - EMBOSS
 - BLAST+

 ## Clustal Omega installation:
 ### Linux
 ```bash
 sudo apt install clustalo
```
 ### Mac
  - install argtable first
  ```bash
brew install argtable
```
  - download clustal omega source code and unpack it:
  ```bash
wget -qO- http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz > clustal-omega-1.2.1.tar.gz
tar zxvf clustal-omega-1.2.1.tar.gz
```
  - installation
  ```bash
  cd clustal-omega-1.2.1
  ./configure CFLAGS='-I/opt/homebrew/include' LDFLAGS='-L/opt/homebrew/lib'
  make
  make install
```

## EMBOSS
### Linux
```bash
sudo apt-get install emboss
```
### Mac 
  ```bash
brew install brewsci/bio/emboss
```

## Blast+
### Linux
```bash
sudo apt-get install ncbi-blast+
```
### Mac
  ```bash
brew install blast
```

# Usage
```bash
python3 pocket_definition.py [-h] -f FAMILY -n NAME -p PDB_STRUCTURE
```
