#!/bin/bash

# Directory to scan (current directory by default)
directory="${1:-.}"

# Loop through all files in the directory
find "$directory" -type f | while read -r file; do
    # Check for the presence of Cyrillic characters (range: \u0400-\u04FF)
    grep -n '[А-Яа-яЁё]' "$file" | while read -r line; do
        # Extract line number and the Cyrillic characters found
        line_number=$(echo "$line" | cut -d: -f1)
        cyrillic_chars=$(echo "$line" | cut -d: -f2-)
        
        # Print file path, line number, and the Cyrillic characters found
        echo "File: $file, Line: $line_number"
        echo "Cyrillic characters: $cyrillic_chars"
        echo "----------------------------------------"
    done
done
