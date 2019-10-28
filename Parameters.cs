using System;
using System.Windows.Forms;

namespace RefactoredProject
{
    public partial class Parameters : Form
    {

        public event InitialTextChanged numberofbodies;
        public event InitialTextChanged dimension;
        public event InitialTextChanged numberofiterations;
        public event InitialCheckChanged unevencoordinates;
        public event InitialCheckChanged tangentialvelocities;
        public event InitialCheckChanged enableblackhole;
        public event InitialComboxBoxItemChanged comboboxitemchanged;
        public Parameters()
        {
            InitializeComponent();
            this.Icon = RefactoredProject.Properties.Resources.input_icon;
        }

        public void Default()
        {
            this.MakeVelocitiesTangentCheckbox.Checked = false;
            this.UnevenCoordinatesCheckbox.Checked = false;
            this.EnterDimensionTextbox.Text = "";
            this.EnterNumberOfBodiesTextbox.Text = "";
            this.EnterNumberOfIterationsTextbox.Text = "";
            this.CalculationModuleCombobox.SelectedIndex = 0;
        }

        private void EnterNumberOfBodiesTextbox_TextChanged(object sender, EventArgs e)
        {
            numberofbodies(this, this.EnterNumberOfBodiesTextbox.Text);
        }

        private void EnterDimensionTextbox_TextChanged(object sender, EventArgs e)
        {
            dimension(this, this.EnterDimensionTextbox.Text);

        }

        private void EnterNumberOfIterationsTextbox_TextChanged(object sender, EventArgs e)
        {
            numberofiterations(this, this.EnterNumberOfIterationsTextbox.Text);
        }

        private void UnevenCoordinatesCheckbox_CheckedChanged(object sender, EventArgs e)
        {
            unevencoordinates(this, new CheckedEventArgs(this.UnevenCoordinatesCheckbox.Checked));
        }

        private void MakeVelocitiesTangentCheckbox_CheckedChanged(object sender, EventArgs e)
        {
            tangentialvelocities(this, new CheckedEventArgs(this.MakeVelocitiesTangentCheckbox.Checked));
        }

        private void CalculationModuleCombobox_SelectedIndexChanged(object sender, EventArgs e)
        {
            comboboxitemchanged(this, CalculationModuleCombobox.SelectedIndex == 0);           
        }

        private void enable_blackhole_CheckedChanged(object sender, EventArgs e)
        {
            enableblackhole(this, new CheckedEventArgs(this.enable_blackhole.Checked));
        }
    }
}
