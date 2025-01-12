#!/bin/bash

# Input file (passed as the first argument)
input_file="$1"
output_name="$2"
output_id="$3"

# Initialize the breseq command
breseq_command="breseq"

# Arrays to store references and input files
reference_files=()
junction_reference_files=()
input_files=()

# Read form-data.txt line by line
while IFS= read -r line || [[ -n "$line" ]]; do
    # Extract key and value
    key=$(echo "$line" | cut -d':' -f1 | xargs)     # Trim whitespace
    value=$(echo "$line" | cut -d':' -f2- | xargs)  # Trim whitespace

    # Handle key-value pairs
    case "$key" in
        "Reference")
            reference_files+=("$value")
            ;;
        "Field")
            breseq_command+=" --genbank-field-for-seq-id $value"
            ;;
        "Junction Reference")
            junction_reference_files+=("$value")
            ;;
        "name")
            breseq_command+=" -n $value"
            ;;
        "limit-fold-coverage")
            [[ "$value" != "OFF" ]] && breseq_command+=" -l $value"
            ;;
        "read-min-length")
            breseq_command+=" --read-min-length $value"
            ;;
        "read-max-same-base-fraction")
            breseq_command+=" --read-max-same-base-fraction $value"
            ;;
        "read-max-N-fraction")
            breseq_command+=" --read-max-N-fraction $value"
            ;;
        "minimum-mapping-quality")
            breseq_command+=" --minimum-mapping-quality $value"
            ;;
        "base-quality-cutoff")
            breseq_command+=" --base-quality-cutoff $value"
            ;;
        "require-match-length")
            breseq_command+=" --require-match-length $value"
            ;;
        "require-match-fraction")
            breseq_command+=" --require-match-fraction $value"
            ;;
        "maximum-read-mismatches")
            breseq_command+=" --maximum-read-mismatches $value"
            ;;
        "deletion-coverage-propagation-cutoff")
            breseq_command+=" --deletion-coverage-propagation-cutoff $value"
            ;;
        "deletion-coverage-seed-cutoff")
            breseq_command+=" --deletion-coverage-seed-cutoff $value"
            ;;
        "junction-indel-split-length")
            breseq_command+=" --junction-indel-split-length $value"
            ;;
        "junction-alignment-pair-limit")
            breseq_command+=" --junction-alignment-pair-limit $value"
            ;;
        "junction-minimum-candidates")
            breseq_command+=" --junction-minimum-candidates $value"
            ;;
        "junction-maximum-candidates")
            breseq_command+=" --junction-maximum-candidates $value"
            ;;
        "junction-candidate-length-factor")
            breseq_command+=" --junction-candidate-length-factor $value"
            ;;
        "junction-minimum-candidate-pos-hash-score")
            breseq_command+=" --junction-minimum-candidate-pos-hash-score $value"
            ;;
        "junction-score-cutoff")
            breseq_command+=" --junction-score-cutoff $value"
            ;;
        "junction-minimum-pos-hash-score")
            breseq_command+=" --junction-minimum-pos-hash-score $value"
            ;;
        "junction-minimum-side-match")
            breseq_command+=" --junction-minimum-side-match $value"
            ;;
        "junction-minimum-pr-no-read-start-per-position")
            breseq_command+=" --junction-minimum-pr-no-read-start-per-position $value"
            ;;
        "consensus-score-cutoff")
            breseq_command+=" --consensus-score-cutoff $value"
            ;;
        "consensus-frequency-cutoff")
            breseq_command+=" --consensus-frequency-cutoff $value"
            ;;
        "polymorphism-score-cutoff")
            breseq_command+=" --polymorphism-score-cutoff $value"
            ;;
        "polymorphism-frequency-cutoff")
            breseq_command+=" --polymorphism-frequency-cutoff $value"
            ;;
        "Input File")
            input_files+=("$value")
            ;;
    esac

    # Handle flags (lines without values)
    case "$key" in
        "flag")
            case "$value" in
                "polymorphism-prediction")
                    breseq_command+=" -p"
                    ;;
                "aligned-sam")
                    breseq_command+=" --aligned-sam"
                    ;;
                "targeted-sequencing")
                    breseq_command+=" -t"
                    ;;
                "quality-score-trim")
                    breseq_command+=" --quality-score-trim"
                    ;;
                "junction-allow-suboptimal-matches")
                    breseq_command+=" --junction-allow-suboptimal-matches"
                    ;;
                "polymorphism-no-indels")
                    breseq_command+=" --polymorphism-no-indels"
                    ;;
                "skip-RA-MC-prediction")
                    breseq_command+=" --skip-RA-MC-prediction"
                    ;;
                "skip-JC-prediction")
                    breseq_command+=" --skip-JC-prediction"
                    ;;
                "skip-MC-prediction")
                    breseq_command+=" --skip-MC-prediction"
                    ;;
                "cnv")
                    breseq_command+=" --cnv"
                    ;;
                "cnv-ignore-redundant")
                    breseq_command+=" --cnv-ignore-redundant"
                    ;;
            esac
            ;;
    esac
done < "$input_file"

# Add reference files to the command
for ref in "${reference_files[@]}"; do
    breseq_command+=" -r /home/ark/MAB/breseq/$output_id/$ref"
done

breseq_command+=" -o /home/ark/MAB/breseq/completed/$output_id/$output_name"

breseq_command+=" -j 24 --brief-html-output"

# Add junction-only reference files to the command
for junction_ref in "${junction_reference_files[@]}"; do
    breseq_command+=" -s $junction_ref"
done

# Add input files to the command
breseq_command+=" ${input_files[*]}"

# Output the constructed command
echo "$breseq_command"

# Optionally, run the command (uncomment the following line to enable running)
# eval "$breseq_command"
