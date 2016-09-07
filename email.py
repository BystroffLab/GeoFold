import smtplib
msg = "TESTINGINGGININGIN"
me = 'walcob@rpi.edu'
you = 'benjamin.walcott@gmail.com'
s = smtplib.smtp('localhost')
s.sendmail(me,you,msg)
s.quit()
