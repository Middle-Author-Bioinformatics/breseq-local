#!/usr/bin/env python3
import boto3
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import argparse
import os

# Initialize SES client
ses = boto3.client('ses', region_name='us-east-2')

# Function to send email with attachment using SES
def send_email_with_attachment(sender_email, recipient_email, subject, body, attachment_path):
    """
    Sends an email with an attachment using Amazon SES.
    """
    try:
        # Create the email
        msg = MIMEMultipart()
        msg['From'] = sender_email
        msg['To'] = recipient_email
        msg['Subject'] = subject

        # Add body text
        msg.attach(MIMEText(body, 'plain'))

        # Add the attachment
        with open(attachment_path, 'rb') as attachment:
            part = MIMEBase('application', 'octet-stream')
            part.set_payload(attachment.read())
            encoders.encode_base64(part)
            part.add_header(
                'Content-Disposition',
                f'attachment; filename={os.path.basename(attachment_path)}',
            )
            msg.attach(part)

        # Convert the email to raw format
        raw_message = {
            'Data': msg.as_string(),
        }

        # Send the email via SES
        response = ses.send_raw_email(
            Source=sender_email,
            Destinations=[recipient_email],
            RawMessage=raw_message
        )
        print("Email sent successfully!")
        print(response)

    except Exception as e:
        print(f"Failed to send email: {e}")

# Function to send email without attachment using SES
def send_email_without_attachment(sender_email, recipient_email, subject, body):
    """
    Sends an email with an attachment using Amazon SES.
    """
    try:
        # Create the email
        msg = MIMEMultipart()
        msg['From'] = sender_email
        msg['To'] = recipient_email
        msg['Subject'] = subject

        # Add body text
        msg.attach(MIMEText(body, 'plain'))

        # Convert the email to raw format
        raw_message = {
            'Data': msg.as_string(),
        }

        # Send the email via SES
        response = ses.send_raw_email(
            Source=sender_email,
            Destinations=[recipient_email],
            RawMessage=raw_message
        )
        print("Email sent successfully!")
        print(response)

    except Exception as e:
        print(f"Failed to send email: {e}")

# Main function for command-line interface
def main():
    """
    Parses command-line arguments and sends email using SES.
    """
    parser = argparse.ArgumentParser(description="Send email with attachment using Amazon SES.")
    parser.add_argument('--sender', type=str, required=True, help="Sender's email address (must be verified in SES).")
    parser.add_argument('--recipient', type=str, required=True, help="Recipient's email address.")
    parser.add_argument('--subject', type=str, required=True, help="Subject of the email.")
    parser.add_argument('--body', type=str, required=True, help="Body text of the email.")
    parser.add_argument('--attachment', type=str, required=False, help="Path to the attachment file.")
    args = parser.parse_args()

    if args.attachment is None:
        send_email_without_attachment(
            sender_email=args.sender,
            recipient_email=args.recipient,
            subject=args.subject,
            body=args.body,
        )
    else:
        send_email_with_attachment(
            sender_email=args.sender,
            recipient_email=args.recipient,
            subject=args.subject,
            body=args.body,
            attachment_path=args.attachment
        )

if __name__ == "__main__":
    main()
