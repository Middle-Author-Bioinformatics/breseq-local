#!/usr/bin/env python3
import os
import boto3
import argparse
from botocore.exceptions import NoCredentialsError
from botocore.config import Config
import pyshorteners

# Initialize S3 client
# s3 = boto3.client('s3', aws_access_key_id='YOUR_ACCESS_KEY',
#     aws_secret_access_key='YOUR_SECRET_KEY',
#     aws_session_token='YOUR_SESSION_TOKEN', region_name='us-east-2')

s3 = boto3.client('s3', region_name='us-east-2')

def shorten_url(long_url):
    try:
        shortener = pyshorteners.Shortener()
        short_url = shortener.tinyurl.short(long_url)
        return short_url
    except Exception as e:
        print(f"Error shortening URL: {e}")
        return long_url


def generate_presigned_url(bucket_name, object_key, expiration=86400):
    # s3_client = boto3.client('s3')
    s3_client = boto3.client('s3', config=Config(signature_version='s3v4'), region_name='us-east-2') # Enforce AWS Signature Version 4
    try:
        url = s3_client.generate_presigned_url(
            'get_object',
            Params={'Bucket': bucket_name, 'Key': object_key},
            ExpiresIn=expiration
        )
        return url
    except Exception as e:
        print(f"Error generating presigned URL: {e}")
        return None

# Main function
def main():
    parser = argparse.ArgumentParser(description='Upload files or directories to S3')
    parser.add_argument('--bucket', type=str, required=True, help='S3 bucket name')
    parser.add_argument('--key', type=str, required=True, help='Output key in S3 (file or directory prefix)')
    parser.add_argument('--expiration', type=str, required=True, help='Presigned URL expiration in seconds')
    args = parser.parse_args()

    expiration = args.expiration
    bucket_name = args.bucket
    output_key = args.key

    presigned_url = generate_presigned_url(bucket_name, output_key, expiration)
    # Example: After generating the pre-signed URL
    short_url = shorten_url(presigned_url)
    print(short_url)

if __name__ == "__main__":
    main()