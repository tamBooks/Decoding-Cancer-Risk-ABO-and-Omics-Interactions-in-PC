#!/bin/bash

# Path to the VCF file
VCF_FILE="dups_pairs.recode.vcf"

# Path to the text file containing the sample pairs
PAIR_FILE="dups.txt"

# Loop over each line (pair) in the text file
while IFS=' ' read -r sample1 sample2; do
    echo "Processing pair: $sample1 and $sample2"

    # Use awk to extract the column numbers for both samples
    col1=$(awk -F'\t' -v sample="$sample1" '{for(i=1;i<=NF;i++) if ($i==sample) print i}' "$VCF_FILE")
    col2=$(awk -F'\t' -v sample="$sample2" '{for(i=1;i<=NF;i++) if ($i==sample) print i}' "$VCF_FILE")

    # Check if both columns were found
    if [ -z "$col1" ] || [ -z "$col2" ]; then
        echo "Error: One or both samples not found in VCF for pair: $sample1, $sample2"
        continue
    fi

    # Use awk to extract the relevant columns and create a new file for each pair
    awk -v c1="$col1" -v c2="$col2" 'BEGIN {OFS="\t"} NR == 1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $(c1), $(c2)} NR > 1 {gsub(/:.*/, "", $(c1)); gsub(/:.*/, "", $(c2)); print $1, $2, $3, $4, $5, $6, $7, $8, $9, $(c1), $(c2)}' "$VCF_FILE" > "${sample1}_${sample2}.vcf"

    echo "Created file: ${sample1}_${sample2}.vcf"
done < "$PAIR_FILE"

# Loop over each pair VCF file created
for vcf in *_*.vcf; do
    # Extract the GT (genotype) field and save it to a new output file
    awk '
    BEGIN { OFS="\t" }
    /^#/ { print; next }  # Print header lines as they are
    {
        for(i=10; i<=NF; i++) {  # Start from the 10th column (first sample)
            split($i, geno, ":");  # Split each sample field by ":"
            $i = geno[1];  # Replace the sample field with only the genotype
        }
        print
    }' "$vcf" > "filtered_${vcf}"  # Save the filtered result to a new file
done

for vcf in filtered_*.vcf; do
    # Use sed to skip the first 57 lines and then cut the desired fields
    #sed '1,77d; /^[[:space:]]*$/d' "$vcf" | cut -f10-11 > "${vcf}.txt"
    sed '1,136d; /^[[:space:]]*$/d' "$vcf" | cut -f10-11 > "${vcf}.txt"
done

for vcf in *.txt; do
    grep -v '\./\.' "$vcf" > "no_missing_${vcf}"
done

for vcf in filtered_*.txt; do
    head -1 "$vcf" > "with_missing_${vcf}"
    grep '\./\.' "$vcf" >> "with_missing_${vcf}"
done

for file in with_missing_*.txt; do
    echo -n "$file: "  # Print the file name
    awk '$1 == "./." && $2 != "./." || $1 != "./." && $2 == "./."' "$file" | wc -l
done

for file in with_missing_*.txt; do
    echo "$file:"
    echo -n "Column 1 (./.): "
    awk '$1 == "./."' "$file" | wc -l  # Count for the first column
    echo -n "Column 2 (./.): "
    awk '$2 == "./."' "$file" | wc -l  # Count for the second column
    echo
done


for file in with_missing_*.txt; do
    echo "$file:"
    
    echo -n "Column 1 (./.) but Column 2 is not (./.): "
    awk '$1 == "./." && $2 != "./."' "$file" | wc -l  # Count where column 1 is `./.` but column 2 is not

    echo -n "Column 2 (./.) but Column 1 is not (./.): "
    awk '$2 == "./." && $1 != "./."' "$file" | wc -l  # Count where column 2 is `./.` but column 1 is not

    echo
done

for file in `ls no_missing_filtered_*`; do awk 'BEGIN {FS=" "} {if ($1 != $2) print $0}' $file ; done
for file in `ls no_missing_filtered_*`; do awk 'BEGIN {FS=" "} {if ($1 != $2) print $0}' $file | wc -l ; done

vcftools --gzvcf merged_panGen2024.norm.vcf.gz --indv 207769820036_R01C01 --indv 208025880084_R01C01 --indv 208025880086_R02C02 --indv 208025880090_R03C01 --indv 208025880092_R04C02 --indv 208025880109_R05C01 --indv 208025880110_R06C02 --indv 208025880114_R07C01 --indv 208025880116_R08C01 --indv 208026130054_R09C02 --indv 208026130069_R10C01 --recode --recode-INFO-all --out /local/tsaid/trabajo/samples_10_only
for i in {1..9}; do     for j in $(seq $((i+1)) 10); do ; done        

awk -v col1="$i" -v col2="$j" 'BEGIN {FS=" "} {if ($col1 != $col2) print $0}' sample10_disagrement.txt | wc -l ; done; done


for file in no_missing_filtered_*; do
    # Print the current file being processed
    echo "Processing file: $file"
    
    # Run the awk command and count the results
    result=$(awk 'BEGIN {FS=" "} {if ($1 != $2) print $0}' "$file" | wc -l)
    
    # Display the result with the file name
    echo "File: $file - Count: $result"
done


for file in with_missing_*.txt; do
    # Get the column names from the first line (header) of the file
    column1_name=$(head -1 "$file" | awk '{print $1}')
    column2_name=$(head -1 "$file" | awk '{print $2}')

    # Display the file name
    echo "File: $file"
    
    # Display counts for column 1
    echo -n "Column '$column1_name' (./.): "
    awk '$1 == "./."' "$file" | wc -l

    # Display counts for column 2
    echo -n "Column '$column2_name' (./.): "
    awk '$2 == "./."' "$file" | wc -l
    
    echo
done

for file in filtered_*.txt; do
    # Print the current file being processed
    echo "Processing file: $file"
    echo -n "Non-matching rows (e.g., './. 1/0' or '1/0 ./.'): "
    awk '$1 != $2 && (($1 == "./." && $2 != "./.") || ($2 == "./." && $1 != "./."))' "$file" | wc -l
    echo
done

for file in filtered_*.txt; do
    # Print the current file being processed
    echo "Processing file: $file"
    
    # Count non-matching rows where values in the columns are different
    non_matching_count=$(awk '$1 != $2 && (($1 == "./." && $2 != "./.") || ($2 == "./." && $1 != "./."))' "$file" | wc -l)
    
    # Count how many times column 1 is './.' while column 2 is not
    col1_non_matching=$(awk '$1 == "./." && $2 != "./."' "$file" | wc -l)
    
    # Count how many times column 2 is './.' while column 1 is not
    col2_non_matching=$(awk '$2 == "./." && $1 != "./."' "$file" | wc -l)

    echo "Non-matching rows: $non_matching_count"
    echo "Column 1 ('./.' while Column 2 is different): $col1_non_matching"
    echo "Column 2 ('./.' while Column 1 is different): $col2_non_matching"
    
    # Determine which column has more non-matching
    if [ "$col1_non_matching" -gt "$col2_non_matching" ]; then
        echo "Column 1 has more non-matching cases."
    elif [ "$col2_non_matching" -gt "$col1_non_matching" ]; then
        echo "Column 2 has more non-matching cases."
    else
        echo "Both columns have the same number of non-matching cases."
    fi
    
    echo
done
