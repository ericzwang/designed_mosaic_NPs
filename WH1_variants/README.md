# Design of RBD variants using the Wuhan-Hu-1 RBD as a reference

You can uncompress the CSV files in order to get all the data files produced from the code.
```bash
tar -xf csv_files.tar.gz
```

Here are steps for running the code:

(1) Create the class CSV files
```bash
python create_class_csv_files.py
```

(2) Create the paired CSV files. Since this is computationally expensive, this is parallelized and run on a cluster.
```bash
./create_pair_csv_files.scr
```

(3) Concatenate the CSV files into `all.csv`
```bash
./make_all.scr
```

(4) Create input tensors for neural network evaluation
```bash
python create_input_tensors.py
```

(5) Evaluate the neural network to predict expressions
```bash
python predict_expressions.py
```

(6) Evaluate solubility with aggrescan. Having aggrescan installed in your environment is required here.
```bash
cd aggrescan
python make_directories.py
sh submit.sh
cd ..
```

(7) Concatenate solubility scores
```bash
./make_soluble_indices.scr
```

(8) Print final sequences
```bash
python print_sequences.py
```
