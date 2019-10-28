using System;
using System.Drawing;
using System.IO;
namespace RefactoredProject
{
    [Serializable]
    class Calculate
    {

        #region Fields

        [NonSerialized]
        private Graphics gr;
        private double coeff;
        [NonSerialized]
        private Rectangle rect;
        private double currenttime;
        private double deltatime;
        private double[][] currentpositions;
        public double totaltime;
        private double sumanglehorizontal;
        private double deltaanglehorizontal;
        private bool move;
        private bool rotatehorizontal;
        public bool need_numbers = true;
        public bool need_visualisation = true;
        private int speed;
        public SphereStarCluster sphere;
        public int to_be_highlighted = -1;

        #endregion

        #region AccessAndConstruction

        public Calculate(SphereStarCluster sphere)
        {
            this.sphere = sphere;
            speed = 50;
            totaltime = 0;
            currenttime = 0;
            deltatime = 3.1536e7 * speed;
            currentpositions = new double[sphere.counterofpoints][];
            sumanglehorizontal = 0;
            rotatehorizontal = false;
            move = false;
            for (int i = 0; i < sphere.counterofpoints; i++)
            {
                currentpositions[i] = new double[3];
            }
        }

        public void GetStarted()
        {
            for (int i = 0; i < sphere.counterofpoints; i++)
            {
                currentpositions[i][0] = sphere.datastorage[i * 3 + 0];
                currentpositions[i][1] = sphere.datastorage[i * 3 + 1];
                currentpositions[i][2] = sphere.datastorage[i * 3 + 2];
            }
        }

        public int Speed
        {
            get
            {
                return speed;
            }

            set
            {
                if (value > 0 && value <= 121)
                {
                    speed = value;
                    deltatime = 3.1536e7 * speed;
                }

            }
        }

        public bool Move
        {
            get
            {
                return move;
            }
            set
            {
                move = value;
            }
        }

        public double DeltaAngleHorizontal
        {
            get
            {
                return deltaanglehorizontal;
            }
            set
            {
                deltaanglehorizontal = value * 1e-3;
            }
        }

        public bool RotateHorizontal
        {
            get
            { return rotatehorizontal; }
            set
            { rotatehorizontal = value; }
        }

        #endregion

        #region Mathematics

        private double[] FromDecartToCylinderHorizontal(double[] decart)
        {
            double[] cylinder = new double[3];
            cylinder[0] = Math.Sqrt(decart[0] * decart[0] + decart[1] * decart[1]);
            cylinder[1] = Math.Atan2(decart[1], decart[0]);
            cylinder[2] = decart[2];
            return cylinder;
        }

        private void AllFromDecartToCylinderHorizontal(double[][] decart)
        {
            for (int i = 0; i < decart.Length; i++)
            {
                decart[i] = FromDecartToCylinderHorizontal(decart[i]);
            }
        }

        private double[] FromCylinderToDecartHorizontal(double[] cylinder)
        {
            double[] decart = new double[3];
            decart[0] = cylinder[0] * Math.Cos(cylinder[1]);
            decart[1] = cylinder[0] * Math.Sin(cylinder[1]);
            decart[2] = cylinder[2];
            return decart;
        }

        private void AllFromCylinderToDecartHorizontal(double[][] cylinder)
        {
            for (int i = 0; i < cylinder.Length; i++)
            {
                cylinder[i] = FromCylinderToDecartHorizontal(cylinder[i]);
            }
        }

        private void AddAngleHorizontal()
        {
            for (int i = 0; i < sphere.counterofpoints; i++)
            {
                currentpositions[i][1] += deltaanglehorizontal;
            }
            sumanglehorizontal += deltaanglehorizontal;
        }

        private void AddAngleHorizontal(double angle)
        {
            for (int i = 0; i < sphere.counterofpoints; i++)
            {
                currentpositions[i][1] += angle;
            }
        }

        #endregion

        #region Update

        private void GetNextIteration()
        {
            int step = PlaceBorders();
            if (step != -1)
            {
                for (int i = 0; i < sphere.counterofpoints; i++)
                {
                    Interpolate(i, step);
                }
            }
            UpdateTime();
        }

        private void Interpolate(int i, int step)
        {
            int current = step * sphere.counterofpoints + i;
            int next = (step + 1) * sphere.counterofpoints + i;
            double k = (currenttime - sphere.time[step]) / (sphere.time[step + 1] - sphere.time[step]);
            currentpositions[i][0] = (sphere.datastorage[next * 3 + 0] - sphere.datastorage[current * 3 + 0]) * k + sphere.datastorage[current * 3 + 0];
            currentpositions[i][1] = (sphere.datastorage[next * 3 + 1] - sphere.datastorage[current * 3 + 1]) * k + sphere.datastorage[current * 3 + 1];
            currentpositions[i][2] = (sphere.datastorage[next * 3 + 2] - sphere.datastorage[current * 3 + 2]) * k + sphere.datastorage[current * 3 + 2];

        }

        private int PlaceBorders()
        {
            for (int i = 0; i < sphere.time.Length - 1; i++)
            {
                if (currenttime >= sphere.time[i] && currenttime <= sphere.time[i + 1])
                {
                    return i;
                }
            }
            return -1;
        }

        private void UpdateTime()
        {
            if (currenttime + deltatime > sphere.time[sphere.time.Length - 1])
            {
                totaltime += sphere.time[sphere.time.Length - 1] - currenttime;
                currenttime += sphere.time[sphere.time.Length - 1] - currenttime;
            }
            else
            {
                currenttime += deltatime;
                totaltime += deltatime;
            }
        }

        #endregion

        #region Visualisation

        public bool VisualisationStep(Graphics gr, Rectangle rect, double coeff)
        {
            this.gr = gr;
            this.rect = rect;
            this.coeff = coeff;
            if (currenttime < sphere.time[sphere.time.Length - 1])
            {
                if (move)
                {
                    GetNextIteration();
                    AllFromDecartToCylinderHorizontal(currentpositions);
                    AddAngleHorizontal(sumanglehorizontal);
                }
                else
                {
                    AllFromDecartToCylinderHorizontal(currentpositions);
                }
                if (rotatehorizontal)
                {
                    AddAngleHorizontal();
                }
                AllFromCylinderToDecartHorizontal(currentpositions);
                Visualise(rect, coeff, gr);
            }
            else
            {
                currenttime = 0;
                return false;
            }
            return true;

        }

        private float GetScreenY(double y, Rectangle rect)
        {
            double coeff = y / (sphere.radius * 5);
            if (rect.Width >= rect.Height)
            {

                return (float)coeff * rect.Height + rect.Width / 2;
            }
            else
            {

                return (float)coeff * rect.Width + rect.Width / 2;
            }

        }

        private float GetScreenZ(double z, Rectangle rect)
        {
            double coeff = z / (sphere.radius * 5);
            if (rect.Width >= rect.Height)
            {
                return rect.Height / 2 - (float)coeff * rect.Height;
            }
            else
            {
                return rect.Height / 2 - (float)coeff * rect.Width / 2;
            }
        }

        private void Visualise(Rectangle rect, double coeff, Graphics gr)
        {
            if (need_visualisation)
            {
                for (int i = 0; i < this.sphere.counterofpoints; i++)
                {
                    float size = (4 + (float)currentpositions[i][0] / (float)(0.65 * sphere.radius)) * (float)coeff;
                    float screeny = GetScreenY(currentpositions[i][1], rect);
                    float screenz = GetScreenZ(currentpositions[i][2] + currentpositions[i][0] / 5, rect);//check
                    if (screeny <= rect.Size.Width && screeny >= 0 && screenz >= 0 && screenz < rect.Size.Height && size > 0)
                    {
                        if (i != to_be_highlighted)
                        {
                            gr.FillEllipse(new SolidBrush(Color.White), new RectangleF(new PointF(screeny, screenz), new SizeF(size, size)));
                        }
                        else
                        {
                            gr.FillEllipse(new SolidBrush(Color.Green), new RectangleF(new PointF(screeny, screenz), new SizeF(size, size)));
                        }
                        if (need_numbers)
                        {
                            gr.DrawString(Convert.ToString(i), new Font("Arial", 13), new SolidBrush(Color.Red), new PointF(screeny + size * (float)Math.Sqrt(2) / 2, screenz + size * (float)Math.Sqrt(2) / 2));
                        }
                    }
                }
            }
        }
        #endregion
    }
}


    

