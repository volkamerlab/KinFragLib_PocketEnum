# Required Software:
 - Clustal Omega
 - EMBOSS
 - BLAST+

 ## Clustal Omega installation:
 ### Linux
 ```sudo apt install clustalo```
 ### Mac
  - install argtable first
  ```brew install argtable```
  - download clustal omega source code and unpack it:
  ```wget -qO- http://www.clustal.org/omega/clustal-omega-1.2.1.tar.gz > clustal-omega-1.2.1.tar.gz
      tar zxvf clustal-omega-1.2.1.tar.gz```
  - installation
  ```
  cd clustal-omega-1.2.1
  ./configure CFLAGS='-I/opt/homebrew/include' LDFLAGS='-L/opt/homebrew/lib'
  make
  make install```

## EMBOSS
### Linux
```sudo apt-get install emboss```
### Mac 
  ```brew install brewsci/bio/emboss```

## Blast+
### Linux
```sudo apt-get install ncbi-blast+```
### Mac
  ```brew install blast```

# Usage
```bash
python3 pocket_definition.py [-h] -f FAMILY -n NAME -p PDB_STRUCTURE```