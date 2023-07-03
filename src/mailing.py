import smtplib
from email.message import EmailMessage

def send_mail(content: str, subject: str, to: str = "kathakabu@gmail.com"):
    # set your email and password
    # please use App Password
    email_address = "hamsterfliege@gmail.com"
    email_password = "qlhwsqnotattidgg"

    # create email
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = email_address
    msg['To'] = to
    msg.set_content(content)

    # send email
    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
        smtp.login(email_address, email_password)
        smtp.send_message(msg)