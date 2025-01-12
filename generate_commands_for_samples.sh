#!/bin/bash

# Paths to input files and script
metadata_file=$1
form_data_template=$2
generate_command_script="/home/ark/MAB/bin/breseq-local/generate_breseq_command.sh"
generate_command_script_nanopore="/home/ark/MAB/bin/breseq-local/generate_breseq_command_nanopore.sh"
platform=$4
id=$5
# Output directory for generated commands
output_dir=$3
mkdir -p "$output_dir"

# Convert the metadata file to Unix line endings (in case of Windows formatting issues)
awk '{ sub("\r$", ""); print }' "$metadata_file" | tail -n +2 | while IFS=, read -r sample_name read1 read2 rest_of_line; do
    # Skip empty lines
    [[ -z "$sample_name" ]] && continue

    # Create a new form-data.txt for this sample
    sample_form_data="${output_dir}/form-data-${sample_name}.txt"
    cp "$form_data_template" "$sample_form_data"

    # Append the sample-specific read libraries to form-data.txt
    echo "Input File: /home/ark/MAB/breseq/${id}/${read1} /home/ark/MAB/breseq/${id}/${read2}" >> "$sample_form_data"

    # Append the sample-specific name parameter to form-data.txt
    echo "name: breseq_${sample_name}" >> "$sample_form_data"

    # Generate the breseq command for this sample
    echo "Generating command for sample: $sample_name"
    if [ "$platform" == "Illumina" ]; then
        bash "$generate_command_script" "$sample_form_data" "breseq_${sample_name}" "${id}" > "${output_dir}/command-${sample_name}.txt"
        echo "Command for sample $sample_name saved to ${output_dir}/command-${sample_name}.txt"
    elif [ "$platform" == "Nanopore" ]; then
        bash "$generate_command_script_nanopore" "$sample_form_data" "${sample_name}" "${id}" > "${output_dir}/command-${sample_name}.txt"
        echo "Command for sample $sample_name saved to ${output_dir}/command-${sample_name}.txt"
    else
        echo 'platform not recognized'
    fi

done

echo "All commands generated successfully in $output_dir."

