import boto3
import os
import subprocess
import json
from bs4 import BeautifulSoup
import pandas as pd
from send_email import send_email_without_attachment
from gen_presign_url import generate_presigned_url, shorten_url

# Initialize S3 client
s3_client = boto3.client('s3')

# Global list to track downloaded folders
downloaded_folders = []

log_file_path = '/home/ark/MAB/breseq/processed_folders.log'

def prepare_and_upload_reference(contigsFile: str, original_s3_path: str):
    """
    Takes the provided contigsFile, renames it to reference.fa,
    creates a reference.fa.fai index, and uploads both files
    back to the S3 bucket where the contigs came from.

    Parameters
    ----------
    contigsFile : str
        Local path to the contigs FASTA file that already exists on disk.
    original_s3_path : str
        Full S3 URI where the contigs were originally located, e.g.:
        's3://mybucket/userX/job123/contigs.fa'
    """

    s3 = boto3.client("s3")

    # -------------------------------------------------------
    # Derive bucket + prefix from original S3 path
    # -------------------------------------------------------
    if not original_s3_path.startswith("s3://"):
        raise ValueError("original_s3_path must start with s3://")

    no_prefix = original_s3_path.replace("s3://", "")
    bucket = no_prefix.split("/")[0]
    key_parts = no_prefix.split("/")[1:-1]   # everything except last element
    prefix = "/".join(key_parts)
    if prefix != "":
        prefix += "/"   # add trailing slash

    # -------------------------------------------------------
    # Create local reference files
    # -------------------------------------------------------
    reference_file = "reference.fa"
    index_file = "reference.fa.fai"

    # Copy/rename the contigs file
    subprocess.run(["cp", contigsFile, reference_file], check=True)

    # Create the FASTA index
    subprocess.run(["samtools", "faidx", reference_file], check=True)

    # -------------------------------------------------------
    # Upload both files to the same bucket/prefix
    # -------------------------------------------------------
    s3.upload_file(reference_file, bucket, prefix + "reference.fa")
    s3.upload_file(index_file, bucket, prefix + "reference.fa.fai")

    print(f"Uploaded reference.fa and reference.fa.fai to s3://{bucket}/{prefix}")

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
                elif line.startswith("ContigsFile"):
                    contigs = line.strip().split(" ")[1]
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
        contigsFile = os.path.join(folder_path, contigs)

    elif accession != "N/A":
        os.system(f"/home/ark/MAB/bin/breseq-local/bit2local.sh -a {accession} -o {folder_path}")
        reference = accession + ".gb"
        referenceFile = os.path.join(folder_path, reference)

        contigs = accession + ".fa"
        contigsFile = os.path.join(folder_path, contigs)
    else:
        referenceFile = "None"
        contigsFile = "None"

    return email, referenceFile, contigsFile, poly, fwd, rev

def load_seen_folders(log_path):
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            return set(line.strip() for line in f)
    return set()

def append_seen_folder(log_path, user, folder):
    with open(log_path, 'a') as f:
        f.write(user + "/" + folder + '\n')

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
                    if user in ['rkdgarber']:
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

    # # Define mutations to ignore
    # ignored_mutations = {
    #     ("45881", "+G"),
    #     ("985333", "+C"),
    #     ("3447984", "(C)5→3"),
    #     ("3447984", "(C)₅→3"),  # Unicode subscript 5 version if encountered
    #     ("3,447,984", "(C)5→3"),
    #     ("3,447,984", "(C)₅→3"),  # With commas
    # }

    data = []
    for row in mutation_table.find_all("tr", class_="normal_table_row"):
        columns = row.find_all("td")
        if len(columns) >= 6:
            pos = columns[1].get_text(strip=True).replace(",", "")
            mut = columns[2].get_text(strip=True).replace("₅", "5")  # Normalize subscript
            # if (pos, mut) in ignored_mutations:
            #     continue

            entry = {
                "position": pos,
                "mutation": mut,
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
    bucket_name = 'breseqbucket'
    local_base_dir = '/home/ark/MAB/breseq/data'
    folders = list_folders_in_bucket(bucket_name)
    user = list_user(bucket_name)
    output_dir = '/home/ark/MAB/breseq/results'

    seen_folders = load_seen_folders(log_file_path)

    new_folders = [f for f in folders if f not in seen_folders]

    for s3_folder in new_folders:
        append_seen_folder(log_file_path, user, s3_folder)

    # Debugging: Check if folders are retrieved
    print(f"Folders found in bucket: {folders}")
    print(f"New folders to process: {new_folders}")

    for s3_folder in new_folders:
        print(f"Processing S3 folder: {s3_folder}")

        local_folder = os.path.join(local_base_dir, s3_folder)
        print(f"Local folder path: {local_folder}")

        download_s3_folder(bucket_name, s3_folder, local_folder)

        # Extract form data and send notification email
        email, referenceFile, contigsFile, poly, fwd, rev = extract_form_data(local_folder)
        if email:
            subject = f"Data received for your variant-calling analysis"
            body = (
                f"Hi,\n\n"
                "We have received your sequencing data.\n"
                "You will recieve another email once the results are ready.\n\n"
                "If you have any questions, feel free to reach out.\n\n"
                "Best regards,\nMAB Team\n"
            )
            sender_email = "binfo@midauthorbio.com"
            send_email_without_attachment(
                sender_email=sender_email,
                recipient_email=email,
                subject=subject,
                body=body
            )
        else:
            print(f"No valid form-data.txt found or missing email/name in: {local_folder}")


        # fastq_files = find_fastq_files(local_folder)
        print(f"FASTQ files found: {fwd, rev}")

        # Define output directory
        output_dir = os.path.join(output_dir, s3_folder)
        # os.system(f"mkdir -p {output_dir}")
        print(f"Output directory path: {output_dir}")

        # Skip breseq if output_dir already exists
        if not os.path.exists(output_dir):
            print("Output directory does not exist. Running Breseq.")
            if fwd and rev:
                run_breseq_command(local_folder, fwd, rev, output_dir, poly, referenceFile)
        else:
            print("Output directory already exists. Skipping Breseq.")

        outtar = ''
        # Proceed to mutation extraction and coverage calculations
        if os.path.exists(output_dir):
            print(f"Processing output directory: {output_dir}")

            # Mutation file extraction
            mutation_file = os.path.join(output_dir, "mutation_predictions.json")
            extract_mutations(output_dir)
            if os.path.exists(mutation_file):
                print(f"Mutation file exists: {mutation_file}")
            else:
                print(f"Mutation file does not exist: {mutation_file}")

            # Coverage file processing
            # coverage_file = os.path.join(output_dir, "data", "coverage.txt")
            coverage_file = run_samtools_command(output_dir)
            if coverage_file and os.path.exists(coverage_file):
                print(f"Coverage file exists: {coverage_file}")
            else:
                print(f"Coverage file does not exist: {coverage_file}")

            # Averages file creation
            averages_file = os.path.join(output_dir, "averages.csv")
            if coverage_file:
                calculate_coverage_averages(coverage_file, output_dir)
                if os.path.exists(averages_file):
                    print(f"Averages file exists: {averages_file}")
                else:
                    print(f"Averages file does not exist: {averages_file}")

            # Averages file creation
            bam_file = os.path.join(output_dir, "data", "reference.bam")
            if os.path.exists(bam_file):
                print(f"BAM file exists: {bam_file}")
            else:
                print(f"BAM file does not exist: {bam_file}")

            # Averages file creation
            bai_file = os.path.join(output_dir, "data", "reference.bam.bai")
            if os.path.exists(bai_file):
                print(f"BAM file exists: {bai_file}")
            else:
                print(f"BAM file does not exist: {bai_file}")

            # Attempt to upload files to S3
            print("Starting upload to S3...")
            if os.path.exists(output_dir):
                print(f"Uploading output dir: {output_dir.rstrip('/')}")
                # s3_folder = s3_folder.rstrip('/')
                # compress the output folder below
                os.system(f"tar -czf {output_dir.rstrip('/')}.tar.gz -C {output_dir.rstrip('/')} .")
                outtar = f"{output_dir.rstrip('/')}.tar.gz"
                upload_file_to_s3(bucket_name, s3_folder, outtar)
            else:
                print(f"Output directory not found, skipping upload: {output_dir}")

            if os.path.exists(mutation_file):
                print(f"Uploading mutation file: {mutation_file}")
                upload_file_to_s3(bucket_name, s3_folder, mutation_file)
            else:
                print(f"Mutation file not found, skipping upload: {mutation_file}")

            if os.path.exists(coverage_file):
                print(f"Uploading coverage file: {coverage_file}")
                upload_file_to_s3(bucket_name, s3_folder, coverage_file)
            else:
                print(f"Coverage file not found, skipping upload: {coverage_file}")

            if os.path.exists(averages_file):
                print(f"Uploading averages file: {averages_file}")
                upload_file_to_s3(bucket_name, s3_folder, averages_file)
            else:
                print(f"Averages file not found, skipping upload: {averages_file}")

            if os.path.exists(bam_file):
                print(f"Uploading BAM file: {bam_file}")
                upload_file_to_s3(bucket_name, s3_folder, bam_file)
            else:
                print(f"BAM file not found, skipping upload: {bam_file}")

            if os.path.exists(bai_file):
                print(f"Uploading BAI file: {bai_file}")
                upload_file_to_s3(bucket_name, s3_folder, bai_file)
            else:
                print(f"BAI file not found, skipping upload: {bai_file}")

            s3_path_contigs = f"s3://{bucket_name}/{s3_folder}contigs.fa"
            if os.path.exists(contigsFile):
                print(f"Preparing and uploading reference files based on: {contigsFile}")
                prepare_and_upload_reference(contigsFile, s3_path_contigs)
            else:
                print(f"Contigs file not found, skipping upload: {contigsFile}")

        else:
            print(f"Output directory not found, skipping further processing for: {s3_folder}")

        # Cleanup: delete the local folder after processing (commented out for testing)
        # import shutil
        # shutil.rmtree(local_folder)

        # Generate presigned URLs for download links
        download_links = []
        for file_path in [outtar]:
            if os.path.exists(file_path):
                key = os.path.join(s3_folder, os.path.basename(file_path))
                url = generate_presigned_url(bucket_name, key, expiration=86400)
                if url:
                    short_url = shorten_url(url)
                    download_links.append(f"{os.path.basename(file_path)}: {short_url}")

        if email and download_links:
            subject = f"Results ready for your analysis!"
            body = (
                    f"Hi,\n\n"
                    "Your sequencing data has been processed. You can access your results at the links below:\n\n"
                    # "Web Viewer: https://evolvingstem.midauthorbio.com\n\n"
                    "📥 Downloadable Files:\n" +
                    "\n".join(download_links) + "\n\n"
                                                "If you have any questions, feel free to reach out.\n\n"
                                                "Best regards,\nMAB Team"
            )
            send_email_without_attachment(
                sender_email="binfo@midauthorbio.com",
                recipient_email=email,
                subject=subject,
                body=body
            )

        print(f"Completed processing for folder: {s3_folder}")

    print("All folders processed.")

