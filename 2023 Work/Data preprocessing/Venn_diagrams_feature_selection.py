from matplotlib_venn import venn3
from matplotlib import pyplot as plt


infile_rna_30 = open('2023 Work\Data\Preprocessed_data\Feature selection\RNA_30.txt', "r")
RNA_30 = infile_rna_30.read().split("\n")
infile_rna_15 = open('2023 Work\Data\Preprocessed_data\Feature selection\RNA_15.txt', "r")
RNA_15 = infile_rna_15.read().split("\n")
infile_rna_30_niv = open('2023 Work\Data\Preprocessed_data\Feature selection\RNA_30_niv.txt', "r")
RNA_30_niv = infile_rna_30_niv.read().split("\n")
venn3(subsets = [set(RNA_30), set(RNA_15), set(RNA_30_niv)], set_colors=('blue', 'green', 'yellow'), set_labels=('RNA_30', 'RNA_15', 'RNA_30_niv'))
plt.show()