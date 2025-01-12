#!/usr/bin/env python3
import boto3
import os
import logging
from botocore.exceptions import NoCredentialsError, PartialCredentialsError
import subprocess

# S3 bucket name and folder prefix
BUCKET_NAME = 'breseqbucket'
FOLDER_PREFIX = 'breseq-'
LOG_FILE = '/home/ark/MAB/breseq/processed_folders.log'

# Initialize S3 client (using default credentials)
s3 = boto3.client('s3')

def get_processed_folders():
    """Read the list of already processed folders from the log file."""
    if not os.path.exists(LOG_FILE):
        return set()
    with open(LOG_FILE, 'r') as log:
        return set(line.strip() for line in log)

def log_processed_folder(folder_name):
    """Log the processed folder to the log file."""
    with open(LOG_FILE, 'a') as log:
        log.write(f"{folder_name}\n")

def check_for_new_folders():
    """Check the S3 bucket for new folders containing form-data.txt."""
    processed_folders = get_processed_folders()

    try:
        response = s3.list_objects_v2(Bucket=BUCKET_NAME, Prefix=FOLDER_PREFIX, Delimiter='/')
        if 'CommonPrefixes' not in response:
            print("No new folders found.")
            return None

        for prefix_obj in response['CommonPrefixes']:
            folder_name = prefix_obj['Prefix'].rstrip('/')
            if folder_name in processed_folders:
                continue  # Skip already processed folders

            # Check if form-data.txt exists in the folder
            form_data_key = f"{folder_name}/form-data.txt"
            try:
                s3.head_object(Bucket=BUCKET_NAME, Key=form_data_key)
                print(f"Found new folder ready for processing: {folder_name}")
                log_processed_folder(folder_name)
                return folder_name
            except s3.exceptions.ClientError:
                # form-data.txt does not exist yet
                continue

        print("No new folders ready for processing.")
        return None

    except (NoCredentialsError, PartialCredentialsError) as e:
        print(f"Error: {e}")
        return None

def download_folder(bucket_name, folder_name, local_dir):
    """Download all files from a specific folder in an S3 bucket to a local directory."""
    paginator = s3.get_paginator('list_objects_v2')
    try:
        for page in paginator.paginate(Bucket=bucket_name, Prefix=folder_name):
            if 'Contents' not in page:
                print(f"No files found in folder {folder_name}.")
                return
            for obj in page['Contents']:
                key = obj['Key']
                local_path = os.path.join(local_dir, os.path.relpath(key, folder_name))
                os.makedirs(os.path.dirname(local_path), exist_ok=True)
                s3.download_file(bucket_name, key, local_path)
                print(f"Downloaded {key} to {local_path}")
    except Exception as e:
        print(f"Error downloading folder {folder_name}: {e}")

def process_folder(folder_name):
    """Process the folder by downloading its contents and performing analysis."""
    local_dir = f"/home/ark/MAB/breseq/{folder_name}"
    print(f"Processing folder: {folder_name}")
    download_folder(BUCKET_NAME, folder_name, local_dir)

    # Run breseq.sh script with folder_name as argument
    breseq_script_path = '/home/ark/MAB/bin/breseq-local/breseq.sh'
    try:
        subprocess.run([breseq_script_path, folder_name], check=True)
        print(f"Successfully processed folder {folder_name} using breseq.sh.")
    except subprocess.CalledProcessError as e:
        print(f"Error running breseq.sh on folder {folder_name}: {e}")

def main():
    """Main function to check for new folders and trigger processing."""
    folder_name = check_for_new_folders()
    if folder_name:
        process_folder(folder_name)

if __name__ == "__main__":
    main()
