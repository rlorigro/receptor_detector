# receptor_detector
Solution to sequence labeling problem


## Installation
```
git clone https://github.com/rlorigro/receptor_detector.git
cd receptor_detector/
mkdir build
cd build/
cmake ..
make -j 8
```

## Usage


```
Options:
  -h,--help                   Print this help message and exit
  -o,--output_directory TEXT REQUIRED
                              Path to directory where output should be written
  -r,--referece TEXT REQUIRED Path to tsv containing receptor and linker sequences
  -i,--input_reads TEXT REQUIRED
                              Path to tsv containing reads to be split and labeled
  -t,--n_threads UINT REQUIRED
                              Maximum number of threads to use
```

For example:
```
./classify_receptors \
-r /home/ryan/data/broad_app/linker_and_receptor_segments.tsv \
-i /home/ryan/data/broad_app/unlabelled_reads.tsv \
-o test/ \
-t 8
```

Results are dumped to files in the output directory.

## Example output

### Final results
![image](https://user-images.githubusercontent.com/28764332/189575269-fa09aff9-8e2e-4d50-888d-3611b27dc76b.png)

### Intermediate information
![image](https://user-images.githubusercontent.com/28764332/189575145-6be96bfd-1ebd-4ead-970b-df5c90abe9ce.png)


