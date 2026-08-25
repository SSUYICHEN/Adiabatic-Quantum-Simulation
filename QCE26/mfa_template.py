"""
Exchange a long-term AWS key for a temporary 12-hour MFA session.

Setup: replace the MFA device ARN below with your own (find it in the
AWS IAM console under your user's "Security credentials" tab).
Then run:  python3 mfa.py   and paste the printed export lines.
"""
import boto3

MFA_DEVICE_ARN = "arn:aws:iam::<YOUR_ACCOUNT_ID>:mfa/<YOUR_DEVICE_NAME>"

sts = boto3.client("sts", region_name="us-east-1")
resp = sts.get_session_token(
    SerialNumber=MFA_DEVICE_ARN,
    TokenCode=input("Enter current 6-digit MFA code: "),
    DurationSeconds=43200,
)
c = resp["Credentials"]
print("export AWS_ACCESS_KEY_ID=" + c["AccessKeyId"])
print("export AWS_SECRET_ACCESS_KEY=" + c["SecretAccessKey"])
print("export AWS_SESSION_TOKEN=" + c["SessionToken"])