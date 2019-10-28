using System;
using System.Drawing;
using System.Windows.Forms;
using System.Threading;
using System.Runtime.Serialization.Formatters.Binary;
using System.IO;

namespace RefactoredProject
{
    public partial class Nbody : Form
    {

        #region Fields       
        private bool[] flags;
        private bool module = true;
        private bool enable_blackhole = false;
        private int counterofpoints;
        private int iterations;
        private double size;
        private Calculate calculate;
        private Thread th;
        private double coeff;
        private Size lastsize;
        private Point MouseCur;
        private bool lastresult;
        private bool cango;
        private string current_file;
        private bool IsSaved = true;
        private Parameters parameters;
        private SpeedTrackBar speedtrackbar;
        private bool unevencoordinates = false;
        private bool tangentialvelocities = false;
        private bool save_canceled = false;
        private ToolStripItem[] configparameters;
        private ToolStripItem[] configure;
        private ToolStripItem[] animationstart;
        private ToolStripItem[] animationstop;
        private ToolStripItem[] animationsettings;
        private ToolStripItem[] advancedlogsitem;
        private AdvancedLogs advancedlogs;

        #endregion

        #region Construction

        public Nbody()//ok
        {
            InitializeComponent();
            this.ClientSize = new System.Drawing.Size(878, 444);
            this.MaximumSize = new System.Drawing.Size(1920, 1080);
            this.MinimumSize = new System.Drawing.Size(900, 500);
            this.WindowState = FormWindowState.Minimized;
            flags = new bool[3];
            coeff = 1;
            lastresult = true;
            cango = false;
            MouseCur = new Point();
            configparameters = this.menuStrip.Items.Find("initialParametersToolStripMenuItem", true);
            configure = this.menuStrip.Items.Find("configureToolStripMenuItem", true);
            animationstart = this.menuStrip.Items.Find("startToolStripMenuItem", true);
            animationstop = this.menuStrip.Items.Find("stopToolStripMenuItem", true);
            animationsettings = this.menuStrip.Items.Find("animationSpeedToolStripMenuItem1", true);
            advancedlogsitem = this.menuStrip.Items.Find("advancedlogsToolStripMenuItem", true);
            configure[0].Enabled = false;
            animationstart[0].Enabled = false;
            animationstop[0].Enabled = false;
            animationsettings[0].Enabled = false;
            advancedlogsitem[0].Enabled = false;
            this.Icon = RefactoredProject.Properties.Resources.nbody_icon;
        }

        #endregion

        #region Methods

        private void Play()//ok
        {
            while (true)
            {
                while (lastresult)
                {
                    this.MainScreen.Invalidate();
                }
                lastresult = true;
                calculate.sphere.intgr();                             
            }
        }

        private double GetCoeff(Size Current, Size Last)//ok
        {
            double[] coeff = new double[2];
            coeff[0] = Convert.ToDouble(Current.Width) / Last.Width;
            coeff[1] = Convert.ToDouble(Current.Height) / Last.Height;
            return coeff[0] * coeff[1];
        }


        #endregion

        #region CheckedChanged

        private void UnevenCoordinatesCheckChanged(object sender, CheckedEventArgs e)//ok
        {
            this.unevencoordinates = e.check;
        }

        private void TangentialVelocitiesCheckChanged(object sender, CheckedEventArgs e)//ok
        {
            this.tangentialvelocities = e.check;
        }

        private void EnableBlackHoleCheckChanged(object sender, CheckedEventArgs e)
        {
            this.enable_blackhole = e.check;
        }

        #endregion

        #region MethodsForDelegates

        private void EnterNumberofBodiesTextChanged(object sender, string s)
        {
            flags[1] = int.TryParse(s, out counterofpoints) && counterofpoints >= 2;
            if (flags[0] & flags[1] & flags[2])
            {
                this.configure[0].Enabled = true;

            }
            else
            {
                this.configure[0].Enabled = false;
            }
        }//ok

        private void EnterDimensionTextChanged(object sender, string s)
        {
            flags[0] = double.TryParse(s, out size) && size > 0;
            if (flags[0] & flags[1] & flags[2])
            {
                this.configure[0].Enabled = true;
            }
            else
            {
                this.configure[0].Enabled = false;
            }
        }//ok

        private void EnterNumberOfIterationsTextChanged(object sender, string s)
        {
            flags[2] = int.TryParse(s, out iterations) && iterations > 0;
            if (flags[0] & flags[1] & flags[2])
            {
                this.configure[0].Enabled = true;
            }
            else
            {
                this.configure[0].Enabled = false;
            }
        }//ok

        private void IntegrateModuleComboBoxChanged(object sender, bool module)
        {
            this.module = module;
        }
        private void ChangeSpeed(object sender, SpeedEventArgs e)
        {
            if (this.calculate != null)
            {
                IsSaved = false;
                this.calculate.Speed = e.speed;
            }
        }//ok

        #endregion

        #region Visualisation

        private void MainScreen_Paint(object sender, PaintEventArgs e)//ok
        {
            if (cango)
            {
                lastresult = calculate.VisualisationStep(e.Graphics, e.ClipRectangle, coeff);
            }
        }

        private void Nbody_SizeChanged(object sender, EventArgs e)
        {
            if (cango)
            {
                coeff = GetCoeff(this.Size, this.lastsize);
            }
        }//ok
        #endregion

        #region Rotation

        private void MainScreen_MouseDown(object sender, MouseEventArgs e)
        {
            if (calculate != null)
            {
                calculate.RotateHorizontal = true;
                MouseCur.X = e.X;
                MouseCur.Y = e.Y;
            }
        }

        private void MainScreen_MouseUp(object sender, MouseEventArgs e)
        {
            if (calculate != null)
                calculate.RotateHorizontal = false;
        }

        private void MainScreen_MouseMove(object sender, MouseEventArgs e)
        {
            if (calculate != null)
            {
                this.calculate.DeltaAngleHorizontal = e.X - MouseCur.X;
                MouseCur.X = e.X;
                MouseCur.Y = e.Y;
            }

        }

        #endregion Rotation

        #region Menu

        private void enableToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.calculate.need_visualisation = true;
        }

        private void disableToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.calculate.need_visualisation = false;
        }
        private void numbersEnableToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.calculate.need_numbers = true;
        }

        private void numbersDisableToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.calculate.need_numbers = false;
        }

        private void saveToolStripMenuItem_Click(object sender, EventArgs e)//ok
        {
            if (th != null)
                th.Suspend();
            SaveState();
            if (th != null)
                th.Resume();
        }

        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (th != null)
                th.Abort();
            LoadState();
            if (calculate != null)
            {      
                calculate.Move = false;
                cango = true;
                th = new Thread(Play);
                th.IsBackground = true;
                th.Start();
                configparameters = this.menuStrip.Items.Find("initialParametersToolStripMenuItem", true);
                configure = this.menuStrip.Items.Find("configureToolStripMenuItem", true);
                animationstart = this.menuStrip.Items.Find("startToolStripMenuItem", true);
                animationstop = this.menuStrip.Items.Find("stopToolStripMenuItem", true);
                animationsettings = this.menuStrip.Items.Find("animationSpeedToolStripMenuItem1", true);
                configparameters[0].Enabled = false;
                configure[0].Enabled = false;
                advancedlogsitem[0].Enabled = true;
                animationstart[0].Enabled = true;
                animationstop[0].Enabled = true;
                animationsettings[0].Enabled = true;
            }

        }


        private void animationSpeedToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            if (speedtrackbar == null || speedtrackbar.IsDisposed)
            {
                speedtrackbar = new SpeedTrackBar(this.calculate.Speed);
                speedtrackbar.speedchanged += new SpeedChanged(this.ChangeSpeed);
            }
            else
            {
                if (speedtrackbar.WindowState == FormWindowState.Minimized)
                {
                    speedtrackbar.WindowState = FormWindowState.Normal;
                }
                speedtrackbar.Activate();

            }
            speedtrackbar.Show();
        }

        private void startToolStripMenuItem_Click(object sender, EventArgs e)
        {
            IsSaved = false;
            calculate.Move = true;
        }

        private void stopToolStripMenuItem_Click(object sender, EventArgs e)
        {
            calculate.Move = false;
        }

        private void configureToolStripMenuItem_Click(object sender, EventArgs e)
        {
            IsSaved = false;
            lastsize = this.Size;
            calculate = new Calculate(new SphereStarCluster(counterofpoints, iterations, size,module,enable_blackhole));
            calculate.sphere.flagcoordinates= this.unevencoordinates;
            calculate.sphere.makevelocitiestangent= this.tangentialvelocities;
            calculate.sphere.intgr();
            calculate.GetStarted();
            cango = true;
            th = new Thread(Play);
            th.IsBackground = true;
            th.Start();
            configparameters[0].Enabled = false;
            configure[0].Enabled = false;
            advancedlogsitem[0].Enabled = true;
            animationsettings[0].Enabled = true;
            animationstart[0].Enabled = true;
            animationstop[0].Enabled = true;
        }

        private void saveAsToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (th != null)
                th.Suspend();
            string filename;
            SaveFileDialog sf = new SaveFileDialog();
            if (sf.ShowDialog() == DialogResult.Cancel)
            {
                th.Resume();
                return;
            }
            else
            {
                filename = sf.FileName;
                current_file = filename;
            }
            BinaryFormatter bf = new BinaryFormatter();
            FileStream file = new FileStream(@filename, FileMode.Create, FileAccess.Write);
            bf.Serialize(file, calculate);
            bf.Serialize(file, lastsize);
            IsSaved = true;
            file.Close();
            if (th != null)
                th.Resume();
        }

        private void newToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (th != null)
                th.Suspend();
            DialogResult result = MessageBox.Show("Cохранить изменения?", "", MessageBoxButtons.YesNoCancel);
            if (!IsSaved && result == DialogResult.Yes)
            {
                this.saveToolStripMenuItem_Click(this, new EventArgs());
            }
            else if (result == DialogResult.Cancel)
            {
                if (th != null)
                {
                    th.Resume();
                }
                return;
            }
            if (!save_canceled)
            {
                th = null;
                IsSaved = true;
                calculate = null;
                configparameters[0].Enabled = true;
                configure[0].Enabled = false;
                animationstart[0].Enabled = false;
                animationstop[0].Enabled = false;
                animationsettings[0].Enabled = false;
                advancedlogsitem[0].Enabled = false;
                cango = false;
                MainScreen.Invalidate();
            }
            else
            {
                save_canceled = false;
            }

        }

        private void advancedlogsToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (advancedlogs == null || advancedlogs.IsDisposed)
            {
                advancedlogs = new AdvancedLogs();
                advancedlogs.getadvancedlogsquery += new GetAdvandcedLogsQueryCallback(GetAdvancedLogsQuery);
            }
            else
            {
                if (advancedlogs.WindowState == FormWindowState.Minimized)
                {
                    advancedlogs.WindowState = FormWindowState.Normal;
                }
                advancedlogs.Activate();
            }
            advancedlogs.Show();
        }

        private void initialParametersToolStripMenuItem_Click(object sender, EventArgs e)
        {
            flags[0] = flags[1] = flags[2] = false;
            if (parameters == null || parameters.IsDisposed)
            {
                parameters = new Parameters();
                parameters.dimension += new InitialTextChanged(this.EnterDimensionTextChanged);
                parameters.numberofbodies += new InitialTextChanged(this.EnterNumberofBodiesTextChanged);
                parameters.numberofiterations += new InitialTextChanged(this.EnterNumberOfIterationsTextChanged);
                parameters.tangentialvelocities += new InitialCheckChanged(this.TangentialVelocitiesCheckChanged);
                parameters.unevencoordinates += new InitialCheckChanged(this.UnevenCoordinatesCheckChanged);
                parameters.comboboxitemchanged += new InitialComboxBoxItemChanged(this.IntegrateModuleComboBoxChanged);
                parameters.enableblackhole += new InitialCheckChanged(this.EnableBlackHoleCheckChanged);
            }
            else
            {
                parameters.Default();
                if (parameters.WindowState == FormWindowState.Minimized)
                {
                    parameters.WindowState = FormWindowState.Normal;
                }
                parameters.Activate();
            }
            parameters.Show();
        }

        public void GetAdvancedLogsQuery()
        {
            if (calculate != null)
            {
                advancedlogs.ChangeText("Кинетическая энергия системы " + Convert.ToString(calculate.sphere.TotalKineticEnergy()) + " Дж", advancedlogs.KineticEnergy);
                advancedlogs.ChangeText("Потенциальная энергия системы " + Convert.ToString(calculate.sphere.TotalPotentialEnergy()) + " Дж", advancedlogs.PotentialEnergy);
                advancedlogs.ChangeText("Время эволюции системы " + Convert.ToString(calculate.totaltime / 3.1536e7) + " Лет", advancedlogs.Time);
                advancedlogs.ChangeText("Полная энергия системы " + Convert.ToString(calculate.sphere.TotalKineticEnergy() + calculate.sphere.TotalPotentialEnergy() + this.calculate.sphere.mergeenergy) + " Дж", advancedlogs.TotalEnergy);
                advancedlogs.ChangeText("Количество звезд в системе " + Convert.ToString(calculate.sphere.counterofpoints-calculate.sphere.evaporated.Count), advancedlogs.NumberofBodiesLabel);
                string evaporated = "Испарившиеся звезды: ";
                for (int i = 0; i < this.calculate.sphere.evaporated.Count; i++)
                {
                    evaporated += Convert.ToString(this.calculate.sphere.evaporated[i]) + ", ";
                }
                advancedlogs.ChangeText(evaporated, advancedlogs.EvaporatedStartsLabel);
                advancedlogs.ChangeText("Минимальное расстояние между звездами (" + this.calculate.sphere.minimal_indexi + "," + this.calculate.sphere.minimal_indexj + "): " + this.calculate.sphere.min_dst + "\n" +
                    "Максимальная скорость (в паре) у звезды " + this.calculate.sphere.max_velocity + " под номером " + this.calculate.sphere.max_velocity_index + "\n" +
                   "Шаг по времени " + 1e-2*(this.calculate.sphere.min_dst / this.calculate.sphere.max_velocity) / 3.1536e7 +" лет", this.advancedlogs.dt_parameters);
            }

        }

        #endregion Menu

        #region FileManagement

        private void getCoorsFileToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.calculate.sphere.Write();
        }

        private void SaveState()
        {
            string filename;
            if (current_file == null)
            {
                SaveFileDialog sf = new SaveFileDialog();
                if (sf.ShowDialog() == DialogResult.Cancel)
                {
                    save_canceled = true;
                    return;
                }
                else
                {
                    filename = sf.FileName;
                    current_file = filename;
                }
            }
            else
            {
                filename = current_file;
            }
            BinaryFormatter bf = new BinaryFormatter();
            FileStream file = new FileStream(@filename, FileMode.Create, FileAccess.Write);
            bf.Serialize(file, calculate);
            bf.Serialize(file, lastsize);
            IsSaved = true;
            file.Close();

        }

        private void LoadState()
        {
            string filename;
            OpenFileDialog od = new OpenFileDialog();
            if (od.ShowDialog() == DialogResult.Cancel)
            {
                return;
            }
            else
            {
                filename = od.FileName;
            }
            BinaryFormatter bf = new BinaryFormatter();
            FileStream file = new FileStream(@filename, FileMode.Open, FileAccess.Read);
            this.calculate = (Calculate)bf.Deserialize(file);
            this.lastsize = (Size)bf.Deserialize(file);
            current_file = filename;
            IsSaved = true;
            file.Close();

        }

        private void Nbody_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (th != null)
                th.Suspend();
            DialogResult result = MessageBox.Show("Cохранить изменения?", "", MessageBoxButtons.YesNoCancel);
            if (!IsSaved && result == DialogResult.Yes)
            {
                this.saveToolStripMenuItem_Click(this, new EventArgs());
            }
            else if (result == DialogResult.Cancel)
            {
                e.Cancel = true;
                if (th != null)
                {
                    th.Resume();
                }
                return;
            }
        }








        #endregion FileManagement      

        private void toolStripTextBox1_TextChanged(object sender, EventArgs e)
        {
            int a=-1;
            int.TryParse((this.menuStrip.Items.Find("toolStripTextBox1", true)[0] as ToolStripTextBox).Text, out a);
            if(a>=0&&a<this.calculate.sphere.counterofpoints)
            {
                this.calculate.to_be_highlighted = a;
            }
        }

        private void boostToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.calculate.sphere.boost = true;
        }
    }
    public delegate void GetAdvandcedLogsQueryCallback();
    public delegate void InitialComboxBoxItemChanged(object sender, bool module);
    public delegate void InitialTextChanged(object sender, string text);
    public delegate void InitialCheckChanged(object sender, CheckedEventArgs e);
    public delegate void SpeedChanged(object sender, SpeedEventArgs e);

    public class SpeedEventArgs : EventArgs
    {
        public int speed;
        public SpeedEventArgs(int speed)
        {
            this.speed = speed;
        }
    }

    public class CheckedEventArgs : EventArgs
    {
        public bool check;
        public CheckedEventArgs(bool check)
        {
            this.check = check;
        }
    }
}
