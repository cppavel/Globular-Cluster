using System;
using System.Windows.Forms;

namespace RefactoredProject
{
    public partial class AdvancedLogs : Form
    {
        public GetAdvandcedLogsQueryCallback getadvancedlogsquery;
        public AdvancedLogs()
        {
            InitializeComponent();
            this.Icon = RefactoredProject.Properties.Resources.information_icon;
        }
        delegate void ChangeTextCallback(string text, Label label);//ok
        public void ChangeText(string text, Label label)
        {
            if (label.InvokeRequired)
            {
                this.Invoke(new ChangeTextCallback(ChangeText), new object[] { text, label });
            }
            else
            {
                label.Text = text;
            }

        }//ok

        private void RefreshButton_Click(object sender, EventArgs e)
        {
           getadvancedlogsquery();
        }
    }
}
