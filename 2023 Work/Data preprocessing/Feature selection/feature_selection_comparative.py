
biogrid_found_genes = "2023 Work\Graph networks\Biogrid\Data\genes_found.txt"
infile = open(biogrid_found_genes, "r")
all_genes_found = infile.read().split("\n")
infile.close()

comparative_genes = "2023 Work/Data/Preprocessed_data/Feature selection/SFS/LR_CV_0/Nuevos datos RNA/300_features.txt"
infile = open(comparative_genes, "r")
comparative_genes = infile.read().split("\n")
infile.close()

included = []
for g in comparative_genes:
    if g in all_genes_found:
        print(g)
        included.append(g)

print("Number of biogrid genes: " + str(len(all_genes_found)))
print("Number of selection genes to compare: " + str(len(comparative_genes)))

print("Included genes: " + str(len(included)))
print(included)