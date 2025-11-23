import boto3
import os
import subprocess
import json
from bs4 import BeautifulSoup
import pandas as pd
from send_email import send_email_without_attachment
from gen_presign_url import generate_presigned_url, shorten_url
import argparse

# Initialize S3 client
s3_client = boto3.client('s3')


def extract_form_data(folder_path):
    form_file = os.path.join(folder_path, "form-data.txt")
    email = reference = accession = poly = None
    if os.path.exists(form_file):
        with open(form_file, "r") as f:
            for line in f:
                if line.startswith("Email"):
                    email = line.strip().split(" ", 1)[1]
                elif line.startswith("ReferenceFile"):
                    reference = line.strip().split(" ")[1]
                # elif line.startswith("ContigsFile"):
                #     contigs = line.strip().split(" ")[1]
                elif line.startswith("Accession"):
                    accession = line.strip().split(" ")[1]
                elif line.startswith("Polymorphic"):
                    poly = " ".join(line.strip().split(" ", 2)[1:])
                elif line.startswith("ForwardReads"):
                    fwd = line.strip().split(" ")[1]
                elif line.startswith("ReverseReads"):
                    rev = line.strip().split(" ")[1]

    if reference != "N/A":
        referenceFile = os.path.join(folder_path, reference)
        # contigsFile = os.path.join(folder_path, contigs)

    elif accession != "N/A":
        os.system(f"/home/ark/MAB/bin/breseq-local/bit2local.sh -a {accession} -o {folder_path}")
        reference = accession + ".gb"
        referenceFile = os.path.join(folder_path, reference)

        contigs = accession + ".fa"
        contigsFile = os.path.join(folder_path, contigs)
    else:
        referenceFile = "None"
        contigsFile = "None"

    return email, referenceFile, poly, fwd, rev


def load_seen_folders(log_path):
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            return set(line.strip() for line in f)
    return set()


def append_seen_folder(log_path, folder):
    with open(log_path, 'a') as f:
        f.write(folder + '\n')


def list_folders_in_bucket(bucket_name):
    paginator = s3_client.get_paginator('list_objects_v2')
    response_iterator = paginator.paginate(Bucket=bucket_name)

    folders = set()
    for page in response_iterator:
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("form-data.txt"):
                parts = key.split('/')
                print(f"Processing key: {key}")
                if len(parts) >= 3:
                    user = parts[0]
                    subfolder = parts[1]
                    if user in ['ark', 'vaughn.cooper']:
                        print(user)
                        folders.add(f"{user}/{subfolder}/")
    return sorted(folders)


def list_user(bucket_name):
    paginator = s3_client.get_paginator('list_objects_v2')
    response_iterator = paginator.paginate(Bucket=bucket_name)
    user = 'unknown'

    for page in response_iterator:
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("form-data.txt"):
                parts = key.split('/')
                print(f"Processing key: {key}")
                if len(parts) >= 3:
                    user = parts[0]
    return user


def download_s3_folder(bucket_name, s3_folder, local_dir):
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    paginator = s3_client.get_paginator('list_objects_v2')
    response_iterator = paginator.paginate(Bucket=bucket_name, Prefix=s3_folder)

    for page in response_iterator:
        if 'Contents' in page:
            for obj in page['Contents']:
                s3_file_path = obj['Key']
                local_file_path = os.path.join(local_dir, os.path.relpath(s3_file_path, s3_folder))
                local_file_dir = os.path.dirname(local_file_path)

                if not os.path.exists(local_file_dir):
                    os.makedirs(local_file_dir)

                s3_client.download_file(bucket_name, s3_file_path, local_file_path)

    downloaded_folders.append(s3_folder)


def find_fastq_files(folder_path):
    fastq_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path)
                   if f.endswith('.fastq') or f.endswith('.fastq.gz')]
    return fastq_files


def run_breseq_command(folder_path, fwd, rev, output_dir, poly, gbk_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fwd = os.path.join(folder_path, fwd)
    rev = os.path.join(folder_path, rev)

    fastq_files_str = fwd + " " + rev
    if poly == "clonal":
        command = f"breseq -l 60 -t -j 12 -o {output_dir} -r {gbk_file} {fastq_files_str}"
    else:
        command = f"breseq --polymorphism-prediction -l 60 -t -j 12 -o {output_dir} -r {gbk_file} {fastq_files_str}"

    full_command = ['/home/ark/miniconda3/bin/conda', 'run', '-n', 'breseq_env', 'bash', '-c', command]

    result = subprocess.run(full_command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running breseq: {result.stderr}")
    else:
        print(f"Breseq ran successfully for folder: {folder_path}")


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
    """
    Parse breseq output/index.html and write mutation_predictions.json

    We now:
      * treat column 1 (second <td>) as seq_id
      * write "seq_id" into the JSON
      * keep backward-compatibility with 6-column tables (no seq_id)
    """
    html_file_path = os.path.join(output_dir, "output", "index.html")

    if not os.path.exists(html_file_path):
        print(f"Error: index.html not found at {html_file_path}")
        return

    with open(html_file_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    mutation_table = None
    for table in soup.find_all("table"):
        if table.find("th", string="Predicted mutations"):
            mutation_table = table
            break

    if mutation_table is None:
        print("Error: Could not find mutation table.")
        return

    rows = mutation_table.find_all("tr", class_="normal_table_row")
    print(f"Found {len(rows)} mutation rows")

    data = []

    for row in rows:
        cols = row.find_all("td")
        n = len(cols)
        if n < 6:
            print("Skipping row with too few columns:", n)
            continue

        # Column layout in your breseq HTML:
        # n == 6 : evidence | position | mutation | annotation | gene | description
        # n >= 7: evidence | seq_id  | position | mutation   | annotation | gene | description
        evidence = cols[0].get_text(strip=True)

        if n == 6:
            seq_id = None
            pos = cols[1].get_text(strip=True).replace(",", "")
            mut = cols[2].get_text(strip=True)
            annotation = cols[3].get_text(strip=True)
            gene = cols[4].get_text(strip=True)
            description = cols[5].get_text(strip=True)
        else:
            seq_id = cols[1].get_text(strip=True)
            pos = cols[2].get_text(strip=True).replace(",", "")
            mut = cols[3].get_text(strip=True)
            annotation = cols[4].get_text(strip=True)
            gene = cols[5].get_text(strip=True)
            description = cols[6].get_text(strip=True)

        data.append({
            "evidence": evidence,
            "seq_id": seq_id,
            "position": pos,
            "mutation": mut,
            "annotation": annotation,
            "gene": gene,
            "description": description,
        })

    out_json = os.path.join(output_dir, "mutation_predictions.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)

    print(f"Wrote {len(data)} mutations to {out_json}")


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


def upload_directory_to_s3(bucket_name, s3_folder, local_directory):
    s3_client = boto3.client('s3')

    for root, dirs, files in os.walk(local_directory):
        for filename in files:
            local_file_path = os.path.join(root, filename)

            # Build the S3 key path relative to the directory being uploaded
            relative_path = os.path.relpath(local_file_path, local_directory)
            s3_key = os.path.join(s3_folder, relative_path).replace("\\", "/")

            print(f"Uploading {local_file_path} to s3://{bucket_name}/{s3_key}")
            s3_client.upload_file(local_file_path, bucket_name, s3_key)


def upload_file_to_s3(bucket_name, s3_folder, local_file):
    s3_key = os.path.join(s3_folder, os.path.basename(local_file))
    s3_client.upload_file(local_file, bucket_name, s3_key)
    print(f"Uploaded {local_file} to s3://{bucket_name}/{s3_key}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process breseq output & upload to S3.")
    parser.add_argument("-i", "--input", required=True, help="Breseq output folder")
    parser.add_argument("-b", "--bucket", required=True, help="S3 bucket name")
    parser.add_argument("-u", "--user_folder", required=True, help="User prefix in S3")
    args = parser.parse_args()


    print(f"Processing output directory: {args.input}")

    # Mutation file extraction
    mutation_file = os.path.join(args.input, "mutation_predictions.json")
    extract_mutations(args.input)
    if os.path.exists(mutation_file):
        print(f"Mutation file exists: {mutation_file}")
    else:
        print(f"Mutation file does not exist: {mutation_file}")

    # Coverage file processing
    # coverage_file = os.path.join(output_dir, "data", "coverage.txt")
    coverage_file = run_samtools_command(args.input)
    if coverage_file and os.path.exists(coverage_file):
        print(f"Coverage file exists: {coverage_file}")
    else:
        print(f"Coverage file does not exist: {coverage_file}")

    # Averages file creation
    averages_file = os.path.join(args.input, "averages.csv")
    if coverage_file:
        calculate_coverage_averages(coverage_file, args.input)
        if os.path.exists(averages_file):
            print(f"Averages file exists: {averages_file}")
        else:
            print(f"Averages file does not exist: {averages_file}")

    # Averages file creation
    bam_file = os.path.join(args.input, "data", "reference.bam")
    if os.path.exists(bam_file):
        print(f"BAM file exists: {bam_file}")
    else:
        print(f"BAM file does not exist: {bam_file}")

    # Averages file creation
    bai_file = os.path.join(args.input, "data", "reference.bam.bai")
    if os.path.exists(bai_file):
        print(f"BAM file exists: {bai_file}")
    else:
        print(f"BAM file does not exist: {bai_file}")

    # Attempt to upload files to S3
    print("Starting upload to S3...")
    if os.path.exists(mutation_file):
        print(f"Uploading mutation file: {mutation_file}")
        upload_file_to_s3(args.bucket, args.user_folder, mutation_file)
    else:
        print(f"Mutation file not found, skipping upload: {mutation_file}")

    if os.path.exists(coverage_file):
        print(f"Uploading coverage file: {coverage_file}")
        upload_file_to_s3(args.bucket, args.user_folder, coverage_file)
    else:
        print(f"Coverage file not found, skipping upload: {coverage_file}")

    if os.path.exists(averages_file):
        print(f"Uploading averages file: {averages_file}")
        upload_file_to_s3(args.bucket, args.user_folder, averages_file)
    else:
        print(f"Averages file not found, skipping upload: {averages_file}")

    if os.path.exists(bam_file):
        print(f"Uploading BAM file: {bam_file}")
        upload_file_to_s3(args.bucket, args.user_folder, bam_file)
    else:
        print(f"BAM file not found, skipping upload: {bam_file}")

    if os.path.exists(bai_file):
        print(f"Uploading BAI file: {bai_file}")
        upload_file_to_s3(args.bucket, args.user_folder, bai_file)
    else:
        print(f"BAI file not found, skipping upload: {bai_file}")

    contigsFile = os.path.join(args.input, "data", "reference.fasta")
    # s3_path_contigs = f"s3://{bucket_name}/{s3_folder}"
    if os.path.exists(contigsFile):
        print(f"Preparing and uploading reference files based on: {contigsFile}")
        subprocess.run(["samtools", "faidx", contigsFile], check=True)
        contigsIndex = contigsFile + ".fai"
        upload_file_to_s3(args.bucket, args.user_folder, contigsFile)
        upload_file_to_s3(args.bucket, args.user_folder, contigsIndex)
    else:
        print(f"Contigs file not found, skipping upload: {contigsFile}")