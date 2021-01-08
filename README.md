# SSBER
the source code of SSBER :removing batch effetcs

required R version > 3.6

The shared cell type labels can be got by SciBet (https://github.com/PaulingLiu/scibet).

[SciBet trains a classifying model based on data with known cell types and predicts cell type for query data. It gives probability distributions of known cell types if the cell type of query data being contained in training data, otherwise it pops out “unassign”]

SSBER consider the cell as a type label if the maximum probability of the label in this cell > 0.8, otherwise the cell is considered as  “unshared”.
