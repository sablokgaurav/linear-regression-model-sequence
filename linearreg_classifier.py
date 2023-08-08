from scipy import stats
import pandas as pd
def expressionLinear(sequence_file, \
                        expression_file, \
                               prediction_file = None, \
                                           arg_type = None):
    """
    a sequence motif training or a sequeunce logo motif trainer based 
    on the surrounding genomic content of the sequences. 
    __sequence_file: the fasta sequences on which you want to train
    __expression_file: a gene expression file for the corresponding genes
    __prediction_file : the fasta sequences on which you want to predict
    __arg_type : either just train or train and predict. In case of the 
    predict, the prediction file is needed.
    """                                       
    if arg_type == "train_model":
        sequence_file_train_read = list(filter(None,[x.strip() \
                                         for x in open(sequence_file).readlines()]))
        sequence_train_dict = {}
        for i in sequence_file_train_read:
            if i.startswith(">"):
                path = i.strip()
                if i not in sequence_train_dict:
                    sequence_train_dict[i] = ""
                    continue
            sequence_train_dict[path] += i.strip()
        ids = list(sequence_train_dict.keys())
        sequences = list(sequence_train_dict.values())
        sequence_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequence"})
        sequence_dataframe["expression"] = pd.read_csv(expression_file)
        sequence_dataframe["genomic_content"] = sequence_dataframe["sequence"]. \
                                  apply(lambda n: (n.count("G")+n.count("C"))/len(n))
        X = sequence_dataframe["expression"]
        y = sequence_dataframe["genomic_content"]
        slope, intercept, r, p, std_err = stats.linregress(X,y)
        linear_model = list(map(lambda n: slope * X + intercept,X))
        return linear_model
    if arg_type == "predict" and prediction_file not None:
        sequence_file_train_read = list(filter(None,[x.strip() \
                                           for x in open(sequence_file).readlines()]))
        sequence_train_dict = {}
        for i in sequence_file_train_read:
            if i.startswith(">"):
                path = i.strip()
                if i not in sequence_train_dict:
                    sequence_train_dict[i] = ""
                    continue
            sequence_train_dict[path] += i.strip()
        ids = list(sequence_train_dict.keys())
        sequences = list(sequence_train_dict.values())
        sequence_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequence"})
        sequence_dataframe["expression"] = pd.read_csv(expression_file)
        sequence_dataframe["genomic_content"] = sequence_dataframe["sequence"]. \
                                  apply(lambda n: (n.count("G")+n.count("C"))/len(n))
        X = sequence_dataframe["expression"]
        y = sequence_dataframe["genomic_content"]
        slope, intercept, r, p, std_err = stats.linregress(X,y)
        linear_model = list(map(lambda n: slope * X + intercept,X))
        prediction_file_predict_read = list(filter(None,[x.strip() \
                                             for x in open(prediction_file).readlines()]))
        prediction_train_dict = {}
        for i in prediction_file_predict_read:
            if i.startswith(">"):
                path = i.strip()
                if i not in prediction_train_dict:
                    prediction_train_dict[i] = ""
                    continue
            prediction_train_dict[path] += i.strip()
        prediction_ids = list(prediction_train_dict.keys())
        prediction_sequences = list(prediction_train_dict.values())
        prediction_dataframe = pd.DataFrame([(i,j) for i,j in 
                                 zip(prediction_ids, prediction_sequences)]). \
                                            rename(columns = {0: "ids", 1: "sequence"})
        prediction_dataframe["genomic_content"] = prediction_dataframe["sequence"]. \
                                    apply(lambda n: (n.count("G") + n.count("C"))/len(n))
        prediction_values = list(map(lambda n: slope * X + intercept, prediction_dataframe["genomic_content"]))
        return prediction_values                            
