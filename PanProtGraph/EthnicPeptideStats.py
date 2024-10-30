import sys
from collections import defaultdict

def process_file(file_path):
    # Initialize data structures
    class_to_peptides = defaultdict(set)
    peptide_to_classes = defaultdict(set)

    # Read the file
    with open(file_path, 'r') as file:
        for line in file:
            class_name, peptide = line.strip().split(";")
            class_to_peptides[class_name].add(peptide)
            peptide_to_classes[peptide].add(class_name)

    # Count peptides shared by different numbers of classes
    shared_counts = defaultdict(int)
    for classes in peptide_to_classes.values():
        shared_counts[len(classes)] += 1

    # Print peptides shared by 1, 2, 3, ... classes
    tot_count = 0
    for num_classes, count in sorted(shared_counts.items()):
        tot_count += count
        print(f"{num_classes}:{count}")
    print(tot_count)

    # output the ethnic group, from which each peptide is 
    with open("ethnic_groups_by_peptides.txt", "w") as f:
        for peptide, classes in peptide_to_classes.items():
            f.write(f"{peptide}:{','.join(classes)}\n")

# Call the function with your file path
process_file(sys.argv[1])

