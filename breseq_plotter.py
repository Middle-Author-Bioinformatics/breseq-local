#!/usr/bin/env python3
import os
import subprocess
import json
from bs4 import BeautifulSoup
import pandas as pd
import argparse

# Global list to track downloaded folders
downloaded_folders = []

def run_samtools_command(output_dir):
    bam_file = os.path.join(output_dir, "data", "reference.bam")
    coverage_file = os.path.join(output_dir, "data", "coverage.txt")

    if not os.path.exists(bam_file):
        print(f"Error: reference.bam not found in {bam_file}")
        return None

    # Run samtools depth command
    command = f"samtools depth -a {bam_file} > {coverage_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running samtools: {result.stderr}")
        return None
    else:
        print(f"Samtools depth ran successfully.")
    return coverage_file


def extract_mutations(output_dir):
    html_file_path = os.path.join(output_dir, "output", "index.html")

    if not os.path.exists(html_file_path):
        print(f"Error: index.html not found in {output_dir}")
        return

    with open(html_file_path, "r", encoding="utf-8") as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    mutation_table = None
    for table in soup.find_all("table"):
        header = table.find("th", string="Predicted mutations")
        if header:
            mutation_table = table
            break

    if not mutation_table:
        print("Error: Could not find the mutation table.")
        return

    data = []
    for row in mutation_table.find_all("tr", class_="normal_table_row"):
        columns = row.find_all("td")
        if len(columns) >= 6:
            entry = {
                "position": columns[1].get_text(strip=True),
                "mutation": columns[2].get_text(strip=True),
                "annotation": columns[3].get_text(strip=True),
                "gene": columns[4].get_text(strip=True),
                "description": columns[5].get_text(strip=True)
            }
            data.append(entry)

    json_file_path = os.path.join(output_dir, "mutation_predictions.json")
    with open(json_file_path, "w", encoding="utf-8") as json_file:
        json.dump(data, json_file, indent=4, ensure_ascii=False)


def calculate_coverage_averages(coverage_file, output_dir):
    averages_file = os.path.join(output_dir, "averages.csv")

    if not os.path.exists(coverage_file):
        print(f"Error: coverage.txt not found in {coverage_file}")
        return

    df = pd.read_csv(coverage_file, sep="\t", header=None, names=["ID", "Index", "Value"])
    averages = []
    chunk_size = 1000

    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i + chunk_size]
        avg = round(chunk["Value"].mean(), 2)
        averages.append((i + 1, i + len(chunk), avg))

    averages_df = pd.DataFrame(averages, columns=["Start Row", "End Row", "Average Value"])
    averages_df.to_csv(averages_file, index=False)
    print(f"Coverage averages saved to {averages_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulate paired-end Illumina reads from a genome FASTA file with optional mutations.")
    parser.add_argument("breseq", type=str, help="Breseq output directory")
    parser.add_argument("output", type=str, help="Output directory for processed files")
    args = parser.parse_args()

    breseq_folder = args.breseq
    output_dir = args.output

    print(f"Local folder path: {breseq_folder}")

    # Mutation file extraction
    mutation_file = os.path.join(output_dir, "mutation_predictions.json")
    extract_mutations(breseq_folder)

    # Coverage file processing
    coverage_file = run_samtools_command(breseq_folder)

    # Averages file creation
    averages_file = os.path.join(output_dir, "averages.csv")
    calculate_coverage_averages(coverage_file, output_dir)







