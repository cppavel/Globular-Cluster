using System;
using System.Windows.Forms;

namespace RefactoredProject
{
    public partial class SpeedTrackBar : Form
    {
        public SpeedTrackBar(int current_speed)
        {
            InitializeComponent();
            trackBar1.Minimum = 10;
            trackBar1.Maximum = 121;
            trackBar1.Value = current_speed;
            this.Icon = RefactoredProject.Properties.Resources.settings_icon;
        }

        public event SpeedChanged speedchanged;

        private void trackBar1_Scroll(object sender, EventArgs e)
        {
            speedchanged(this, new SpeedEventArgs(this.trackBar1.Value));
        }
    }
}
