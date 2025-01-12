#!/usr/bin/env python3
import boto3

# Initialize SES client
ses = boto3.client('ses', region_name='us-east-2')

def is_email_verified(email_address):
    """
    Checks if the email address is already verified in Amazon SES.
    """
    try:
        response = ses.list_verified_email_addresses()
        verified_emails = response.get('VerifiedEmailAddresses', [])
        return email_address in verified_emails
    except Exception as e:
        print(f"Failed to fetch verified email addresses.")
        print(e)
        return False

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Check if an email address is verified in Amazon SES.")
    parser.add_argument(
        '--email',
        type=str,
        required=True,
        help="Email address to check."
    )
    args = parser.parse_args()

    email_address = args.email
    if is_email_verified(email_address):
        print("verified")
    else:
        print("unverified")

if __name__ == "__main__":
    main()