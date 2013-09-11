#!/usr/bin/python

import wx
import smtplib

class MailReport(wx.Dialog):
    """Simple frame that sends mails with bug reports."""
    
    def __init__(self, *args, **kw):
        super(MailReport, self).__init__(*args, **kw)
        
        self.InitUI()
        
    def InitUI(self):
        
        pnl = wx.Panel(self)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)

        st1 = wx.StaticText(pnl, label='From')
        st2 = wx.StaticText(pnl, label='To ')
        st3 = wx.StaticText(pnl, label='Subject')

        self.tc1 = wx.TextCtrl(pnl, size=(180, -1))
        self.tc2 = wx.TextCtrl(pnl, size=(180, -1))
        self.tc3 = wx.TextCtrl(pnl, size=(180, -1))

        self.tc = wx.TextCtrl(pnl, style=wx.TE_MULTILINE)
        button_send = wx.Button(pnl, label='Send')

        hbox1.Add(st1, flag=wx.LEFT, border=10)
        hbox1.Add(self.tc1, flag=wx.LEFT, border=35)
        hbox2.Add(st2, flag=wx.LEFT, border=10)
        hbox2.Add(self.tc2, flag=wx.LEFT, border=50)
        hbox3.Add(st3, flag=wx.LEFT, border=10)
        hbox3.Add(self.tc3, flag=wx.LEFT, border=20)
        vbox.Add(hbox1, flag=wx.TOP, border=10)
        vbox.Add(hbox2, flag=wx.TOP, border=10)
        vbox.Add(hbox3, flag=wx.TOP, border=10)
        vbox.Add(self.tc, proportion=1, flag=wx.EXPAND | wx.TOP |  wx.RIGHT | wx.LEFT, border=15)
        vbox.Add(button_send, flag=wx.ALIGN_CENTER | wx.TOP |  wx.BOTTOM, border=20)

        self.Bind(wx.EVT_BUTTON, self.OnSend, button_send)
        pnl.SetSizer(vbox)

        self.SetSize((400, 420))
        self.SetTitle('Tom')
        self.Centre()
        self.ShowModal()
        self.Destroy()

    def OnSend(self, e):
        
        sender = self.tc1.GetValue()
        recipient = self.tc2.GetValue()
        subject = self.tc3.GetValue()
        text = self.tc.GetValue()
        header = "From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n" %  (sender, recipient, subject)
        message = header + text

        try:
            server = smtplib.SMTP('mail.chello.sk')
            server.sendmail(sender, recipient, message)
            server.quit()
            dlg = wx.MessageDialog(self, 'Email was successfully sent', 'Success', 
                wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal()
            dlg.Destroy()

        except smtplib.SMTPException as error:
            dlg = wx.MessageDialog(self, 'Failed to send email',  'Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()


if __name__ == '__main__':
    ex = wx.App()
    MailReport(None)
    ex.MainLoop()    
