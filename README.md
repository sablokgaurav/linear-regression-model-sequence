# linear-regression-model-sequence
I coded this linear regression based training model based on the sequence features across the sequences. It has two arguments, just train the model or train and predict the model. You can give a short motif file, or any sequences and corresponding expression and it will estimate the linear regression using the genomic based characteristics and scipy.

```python 
expressionLinear("/Users/gauravsablok/Desktop/CodeCheck/fasta_sample_datasets/sample1.fasta", 
                      "/Users/gauravsablok/Desktop/CodeCheck/csv_test_datasets/test_coverage.csv", 
                                                                               arg_type="train_model")
[0     0.397802
 1     0.397802
 2     0.397736
 3     0.396936
 4     0.397736
 5     0.397202
 6     0.397402
 7     0.397069
 8     0.397802
 9     0.397802
 10    0.397736]
 ````                                                                               
Gaurav \
Academic Staff Member \
Bioinformatics \
Institute for Biochemistry and Biology \
University of Potsdam \
Potsdam,Germany
