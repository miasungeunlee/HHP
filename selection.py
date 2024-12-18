#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Python script to process files based on conditions
# Date: 20240321

# Open the file containing the list of filenames
with open('Name-list.txt', 'r') as file_list:
    # Iterate through each filename in the list
    for filename in file_list:
        filename = filename.strip()  # Remove trailing newline characters
        
        # Open the file for reading
        with open(filename, 'r') as file:
            # Read the first line and extract the second field
            line1 = file.readline().split()
            if len(line1) < 2:
                print(f"Error: Insufficient data in {filename}")
                continue
            a1 = int(line1[1])
            
            # Read the second line, extract the second field, and multiply by 3
            line2 = file.readline().split()
            if len(line2) < 2:
                print(f"Error: Insufficient data in {filename}")
                continue
            b1 = int(line2[1]) * 3
            
            # Read all lines and remove duplicates based on the first field
            lines = []
            seen = set()
            for line in file:
                fields = line.split()
                if len(fields) < 1:
                    print(f"Error: Invalid line in {filename}")
                    continue
                if fields[0] not in seen:
                    lines.append(line)
                    seen.add(fields[0])
            
            # Check the condition
            if a1 > b1:
                with open('New_' + filename, 'w') as new_file:
                    # Write the first line to the output file
                    new_file.write(' '.join(line1) + '\n')
                    # Write the rest of the lines
                    new_file.writelines(lines)

