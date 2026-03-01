#!/usr/bin/env python3

# Script that shows each fasta file and each annotation
# file present in the directories campy_fas and campy_ann 
# respectively.

import os

count = 1

directories = [
    "campy_fas",
    "campy_ann"
]

for directory in directories:
    for file in os.listdir(directory):
        filepath = os.path.join(directory, file)

        print(f"File {count}: {file}")
        count += 1

        with open(filepath, "r") as f:
            for line in f:
                if line.startswith(">"):
                    print(line.strip())
                    print("\n")
                    break
