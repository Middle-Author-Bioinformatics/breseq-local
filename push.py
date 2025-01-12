#!/usr/bin/env python3
import os
import boto3
import argparse
from botocore.exceptions import NoCredentialsError

# Initialize S3 client
s3 = boto3.client('s3', region_name='us-east-2')

# Function to upload a single file to S3
def upload_file_to_s3(local_file_path, bucket_name, s3_key):
    try:
        s3.upload_file(local_file_path, bucket_name, s3_key)
        print(f"Uploaded {local_file_path} to bucket {bucket_name} as {s3_key}")
    except NoCredentialsError:
        print("Credentials not available")

# Function to upload an entire directory to S3
def upload_directory_to_s3(local_dir_path, bucket_name, s3_prefix):
    for root, dirs, files in os.walk(local_dir_path):
        for file in files:
            # Construct the full local path
            local_file_path = os.path.join(root, file)
            # Construct the relative path to maintain directory structure in S3
            relative_path = os.path.relpath(local_file_path, local_dir_path)
            s3_key = os.path.join(s3_prefix, relative_path)
            upload_file_to_s3(local_file_path, bucket_name, s3_key)

# Main function
def main():
    parser = argparse.ArgumentParser(description='Upload files or directories to S3')
    parser.add_argument('--bucket', type=str, required=True, help='S3 bucket name')
    parser.add_argument('--output_key', type=str, required=True, help='Output key in S3 (file or directory prefix)')
    parser.add_argument('--source', type=str, required=True, help='Local file or directory to upload')
    args = parser.parse_args()

    source_path = args.source
    bucket_name = args.bucket
    output_key = args.output_key

    if os.path.isdir(source_path):
        print(f"Uploading directory {source_path} to bucket {bucket_name} with prefix {output_key}")
        upload_directory_to_s3(source_path, bucket_name, output_key)
    elif os.path.isfile(source_path):
        print(f"Uploading file {source_path} to bucket {bucket_name} as {output_key}")
        upload_file_to_s3(source_path, bucket_name, output_key)
    else:
        print(f"Source path {source_path} is neither a file nor a directory. Please check the input.")

if __name__ == "__main__":
    main()
