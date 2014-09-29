from __future__ import print_function, division

import os
import wx

from monty.os.path import which
from monty.string import is_string


__all__ = [
    "showErrorMessage",
    "askUser",
    "showLicense",
]


# Helper functions.
def _straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def askUser(parent, message):
    """Open dialog with message, return user's answer."""
    ask = wx.MessageDialog(parent, message)
    answer = ask.ShowModal() == wx.ID_OK
    ask.Destroy()
    return answer


class ErrorDialog(wx.MessageDialog):
    def __init__(self, parent, message):
        super(ErrorDialog, self).__init__(parent, message=message, caption='Error Message',
                                          style=wx.YES_NO | wx.CANCEL | wx.NO_DEFAULT | wx.ICON_ERROR | wx.STAY_ON_TOP)

def showErrorMessage(parent, message=None):
    """
    Open a `MessageDialog` with an error message.
    If message is None, the python traceback is used.
    """
    if message is None: message = _straceback()

    message += "\n\n Do you want to send a bug report?"
    dialog = ErrorDialog(parent, message)

    # Send mail if the user clicked YES.
    if dialog.ShowModal() == wx.ID_YES:
        mail = SendMail(parent)
        mail.setSender(_user_at_host())
        mail.setSubject("Bug report")
        mail.setBody(message)
        mail.ShowModal()
        mail.Destroy()

    dialog.Destroy()


def showLicense(parent=None, codename=None):
    codename = "Abipy" if codename is None else codename

    license_text = """%(codename)s is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

%(codename)s is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have
received a copy of the GNU General Public License along with the code;
if not, write to the Free Software Foundation, Inc., 59 Temple Place,
Suite 330, Boston, MA  02111-1307  USA""" % {"codename": codename}

    dialog = License(parent, license_text)
    dialog.ShowModal()
    dialog.Destroy()


class License(wx.Dialog):
    def __init__(self, parent, license_text,  **kwargs):
        wx.Dialog.__init__ (self, parent, id=-1, title="License")

        vsizer = wx.BoxSizer( wx.VERTICAL )

        text = wx.TextCtrl( self, -1, license_text, style=wx.TE_MULTILINE | wx.TE_READONLY )
        vsizer.Add(text, 0, wx.ALL|wx.EXPAND, 5 )

        vsizer.Add(wx.Button( self, wx.ID_OK), 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        self.SetSizerAndFit(vsizer)


class SendMail(wx.Dialog):
    """This dialog allows the user to send an email with sendmail."""
    def __init__(self, parent, title="SendMail"):
        super(SendMail, self).__init__(parent, -1, title, wx.DefaultPosition, wx.Size(400, 420))

        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        st1 = wx.StaticText(panel, -1, 'From: ')
        st2 = wx.StaticText(panel, -1, 'To: ')
        st3 = wx.StaticText(panel, -1, 'Subject: ')
        self.sender = wx.TextCtrl(panel, -1, size=(180, -1))
        self.mailto = wx.TextCtrl(panel, -1, size=(180, -1))
        self.subject = wx.TextCtrl(panel, -1, size=(180, -1))
        self.body = wx.TextCtrl(panel, -1, style=wx.TE_MULTILINE)
        button_send = wx.Button(panel, 1, 'Send')
        hbox1.Add(st1, 0, wx.LEFT, 10)
        hbox1.Add(self.sender, 0, wx.LEFT, 20)
        hbox2.Add(st2, 0, wx.LEFT, 10)
        hbox2.Add(self.mailto, 0, wx.LEFT, 35)
        hbox3.Add(st3, 0, wx.LEFT, 10)
        hbox3.Add(self.subject, 0)
        vbox.Add(hbox1, 0, wx.TOP, 10)
        vbox.Add(hbox2, 0, wx.TOP, 10)
        vbox.Add(hbox3, 0, wx.TOP, 10)
        vbox.Add(self.body, 1, wx.EXPAND | wx.TOP | wx.RIGHT | wx.LEFT, 15)
        vbox.Add(button_send, 0, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, 20)
        self.Bind(wx.EVT_BUTTON, self.onSend, id=1)
        panel.SetSizer(vbox)
        self.Centre()

    def setSender(self, s):
        self.sender.SetValue(s)

    def setMailto(self, s):
        self.mailto.SetValue(s)

    def setSubject(self, s):
        self.subject.SetValue(s)

    def setBody(self, text):
        self.body.SetValue(text)

    def onSend(self, event):
        sender = self.sender.GetValue()
        mailto = self.mailto.GetValue()
        subject = self.subject.GetValue()
        body = self.body.GetValue()

        header = 'From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n' % (sender, mailto, subject)
        
        try:
            retcode = sendmail(subject, body, mailto, sender=sender)

            if retcode:
                dialog = wx.MessageDialog(self, 'Email was not sent', 'Failure', wx.OK | wx.ICON_INFORMATION)
            else:
                dialog = wx.MessageDialog(self, 'Email was successfully sent', 'Success', wx.OK | wx.ICON_INFORMATION)

            dialog.ShowModal()
            dialog.Destroy()

        except:
            showErrorMessage(self)


def _user_at_host():
    from socket import gethostname
    return os.getlogin() + "@" + gethostname()


def sendmail(subject, text, mailto, sender=None):
    """
    Sends an e-mail with unix sendmail. 

    Args:
        subject:
            String with the subject of the mail.
        text:
            String with the body of the mail.
        mailto:
            String or list of string with the recipients.
        sender:
            string with the sender address.
            If sender is None, username@hostname is used.

    Returns:
        exit status
    """
    # Body of the message.
    sender = _user_at_host() if sender is None else sender
    if is_string(mailto): mailto = [mailto]

    from email.mime.text import MIMEText

    mail = MIMEText(text)
    mail["Subject"] = subject
    mail["From"] = sender
    mail["To"] = ", ".join(mailto)

    msg = mail.as_string()

    # sendmail works much better than the python interface.
    # Note that sendmail is available only on Unix-like OS.
    from subprocess import Popen, PIPE

    sendmail = which("sendmail")
    if sendmail is None: return -1

    p = Popen([sendmail, "-t"], stdin=PIPE, stderr=PIPE)

    outdata, errdata = p.communicate(msg)
    return len(errdata)


class MyApp(wx.App):
    def OnInit(self):
        dialog = SendMail(None)
        dialog.ShowModal()
        dialog.Destroy()
        return True


if __name__  == "__main__":
    MyApp(0).MainLoop()
