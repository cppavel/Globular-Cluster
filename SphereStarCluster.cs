using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.IO;
using System.Text;
namespace RefactoredProject
{
    [Serializable]
    class SphereStarCluster
    {
        [DllImport("RandNum.DLL", EntryPoint = "dRand", CallingConvention = CallingConvention.Cdecl)]
        public extern static double dRand(double left, double right);

        [DllImport("integrate.DLL", EntryPoint = "integrate_", CallingConvention = CallingConvention.Cdecl)]
        public extern static void intgr(ref int cp, ref double x, ref double y, ref double z, ref double vx, ref double vy, ref double vz, ref double m, ref double md, ref int iter, ref double datastorage,
           ref double time, ref double startenergy, ref bool merge, ref double mergeenergy);

        [DllImport("nbodydll64.DLL", EntryPoint = "cudaIntegrate", CallingConvention = CallingConvention.Cdecl)]
        public extern static void cudaIntegrate(ref int cp, double[] x, double[] y, double[] z, double[] vx, double[] vy, double[] vz, double[] m, int iter, double[] datastorage, double[] time,double startenergy,ref bool merge, ref double mergeenergy);

        private double[] vx;
        public List<int> evaporated = new List<int>();
        private double[] vy;
        private double[] vz;
        private double[] x;
        private double[] y;
        private double[] z;
        private double[,] matrixofdistances;
        private double[] masses;
        public double[] datastorage;
        public double[] time;
        public double[] impulse_momentum;
        public double radius;
        private double startenergy;
        public double mergeenergy = 0;
        public int iterations;
        public int counterofpoints;
        public int minimal_indexi = 0;
        public int minimal_indexj = 0;
        public int max_velocity_index = 0;
        public bool flagcoordinates = false;
        private bool firstcall = true;
        public bool makevelocitiestangent = false;
        public bool fortran = true;
        public bool boost = false;
        public bool enable_blackhole = false;
        private bool merge = false;       
        private double massesSigma = 0;
        private double totaltime = 0;
        public double min_dst = 0;
        public double max_velocity = 0;
        

        public SphereStarCluster(int counterofpoints, int iterations, double radius,bool module,bool enable_blackhole)
        {
            this.counterofpoints = counterofpoints;
            this.iterations = iterations;
            this.radius = radius;
            this.fortran = module;
            this.enable_blackhole = enable_blackhole;
            vx = new double[counterofpoints];
            vy = new double[counterofpoints];
            vz = new double[counterofpoints];
            x = new double[counterofpoints];
            y = new double[counterofpoints];
            z = new double[counterofpoints];
            masses = new double[counterofpoints];
            matrixofdistances = new double[counterofpoints, counterofpoints];
            datastorage = new double[(iterations + 1) * counterofpoints * 3];
            time = new double[iterations + 1];
            impulse_momentum = new double[4];

            StreamWriter file_ev = new StreamWriter(@"evaporation.txt");
            file_ev.Write("");
            file_ev.Close();

            StreamWriter file_m = new StreamWriter(@"merge.txt");
            file_m.Write("");
            file_m.Close();

            StreamWriter file_ev_m = new StreamWriter(@"evaporated_masses.txt");
            file_ev_m.Write("");
            file_ev_m.Close();
        }

        private void RandomizeCoordinates()
        {
            double x, y, z;
            double distancefromcenter;
            int createdpoints = 0;
            while (createdpoints < counterofpoints)
            {
                x = radius * (1 - 2 * SphereStarCluster.dRand(0, 1));
                y = radius * (1 - 2 * SphereStarCluster.dRand(0, 1));
                z = radius * (1 - 2 * SphereStarCluster.dRand(0, 1));
                distancefromcenter = Math.Pow(x, 2) + Math.Pow(y, 2) + Math.Pow(z, 2);
                if (Math.Pow(radius, 2) >= distancefromcenter)
                {
                    this.x[createdpoints] = x;
                    this.y[createdpoints] = y;
                    this.z[createdpoints] = z;
                    createdpoints++;
                }
            }
        }

        private void UnevenCoordinates()
        {
            double x, y, z;
            double random;
            double func;
            double distancefromcenter;
            int createdpoints = 0;
            while (createdpoints < counterofpoints)
            {
                x = radius * (1 - 2 * SphereStarCluster.dRand(0, 1));
                y = radius * (1 - 2 * SphereStarCluster.dRand(0, 1));
                z = radius * (1 - 2 * SphereStarCluster.dRand(0, 1));
                distancefromcenter = Math.Pow(x, 2) + Math.Pow(y, 2) + Math.Pow(z, 2);
                random = SphereStarCluster.dRand(0, 1);
                func = Math.Sqrt(distancefromcenter) / radius;
                if (Math.Pow(radius, 2) >= distancefromcenter && random > func)
                {
                    this.x[createdpoints] = x;
                    this.y[createdpoints] = y;
                    this.z[createdpoints] = z;
                    createdpoints++;
                }
            }
        }

        private void RandomizeMasses()
        {
            for (int i = 0; i < counterofpoints; i++)
            {
                masses[i] = SphereStarCluster.dRand(1, 2) * 2e30;
                massesSigma += masses[i];
            }
        }

        private void RandomizeVelocities()
        {
            for (int i = 0; i < counterofpoints; i++)
            {
                vx[i] = SphereStarCluster.dRand(-1, 1) * 1e2;
                vy[i] = SphereStarCluster.dRand(-1, 1) * 1e2;
                vz[i] = SphereStarCluster.dRand(-1, 1) * 1e2;
            }
        }

        private void CalculateDistances()
        {
            for (int i = 0; i < counterofpoints; i++)
            {
                for (int j = i + 1; j < counterofpoints; j++)
                {
                    matrixofdistances[i, j] = Math.Sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]));
                    matrixofdistances[j, i] = matrixofdistances[i, j];
                }
            }
        }


        private void ChangeEnergy()
        {
            double coeff = Math.Sqrt(Math.Abs(TotalPotentialEnergy() / TotalKineticEnergy()) * 0.5);
            for (int i = 0; i < counterofpoints; i++)
            {
                vx[i] *= coeff;
                vy[i] *= coeff;
                vz[i] *= coeff;
            }
        }

        private void ChangeSpeedTangent()
        {
            double[] centermass = GetCenterMass();
            double[] velocity_tangent;
            double[] current_velocity;
            double[] current_coordinate;
            double[] velocity_parallel;
            double coeff = 1;
            for (int i = 0; i < counterofpoints; i++)
            {
                current_coordinate = new double[3] { x[i], y[i], z[i] };
                current_velocity = new double[3] { vx[i], vy[i], vz[i] };
                coeff = ScalarProduct(current_velocity, Residual(current_coordinate, centermass)) / ScalarProduct(Residual(current_coordinate, centermass), Residual(current_coordinate, centermass));
                velocity_parallel = Multiplication(Residual(current_coordinate, centermass), coeff);
                velocity_tangent = Residual(current_velocity, velocity_parallel);
                coeff = Math.Sqrt(ScalarProduct(current_velocity, current_velocity) / ScalarProduct(velocity_tangent, velocity_tangent));
                velocity_tangent = Multiplication(velocity_tangent, coeff);
                vx[i] = velocity_tangent[0];
                vy[i] = velocity_tangent[1];
                vz[i] = velocity_tangent[2];
            }



        }

        public void Firststep()
        {
            if (flagcoordinates)
            {
                UnevenCoordinates();
            }
            else
            {
                RandomizeCoordinates();
            }
            CalculateDistances();
            RandomizeMasses();
            RandomizeVelocities();
            MakeSystemStatic();
            if (makevelocitiestangent)
            {
                ChangeSpeedTangent();
            }
            ChangeEnergy();
            startenergy = TotalKineticEnergy() + TotalPotentialEnergy();
            ConservationImpulse();
            ConservationEnergy();
            startenergy = TotalKineticEnergy() + TotalPotentialEnergy();    
            if(enable_blackhole)
            {
                x[0] = 0;
                y[0] = 0;
                z[0] = 0;
                vx[0] = 0;
                vy[0] = 0;
                vz[0] = 0;
                masses[0] *= 5000;
            }
        }


        public void Write()
        {
            StreamWriter file = new StreamWriter("coors.txt");
            for (int i = 0; i < counterofpoints; i++)
            {
                file.WriteLine("{0} {1} {2}", x[i], y[i], z[i]);              
            }    
            file.Close();
        }

        private double MinimalDistance()
        {
            double min = matrixofdistances[0, 1];
            for (int i = 0; i < counterofpoints; i++)
            {
                for (int j = i + 1; j < counterofpoints; j++)
                {
                    if (matrixofdistances[j, i] < min)
                    {
                        min = matrixofdistances[j, i];
                        minimal_indexi = i;
                        minimal_indexj = j;
                    }
                }
            }
            return min;
        }

        private double GetVelocity(int i)
        {
            return Math.Sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
        }

        private double MaxVelocity()
        {
            double max = Math.Sqrt(vx[0] * vx[0] + vy[0] * vy[0] + vz[0] * vz[0]);
            for (int i = 0; i < counterofpoints; i++)
            {
                if (Math.Sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]) > max)
                {
                    max = Math.Sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
                    max_velocity_index = i;
                }

            }
            return max;
        }

        public double TotalPotentialEnergy()
        {
            double sum = 0;
            for (int i = 0; i < counterofpoints; i++)
            {
                for (int j = 0; j < counterofpoints; j++)
                {
                    if (i != j)
                    {
                        sum += -6.67408e-11 * masses[i] * masses[j] / matrixofdistances[i, j];
                    }
                }
            }
            return sum / 2;

        }

        public double TotalKineticEnergy()
        {
            double sum = 0;
            for (int i = 0; i < counterofpoints; i++)
            {
                sum += masses[i] * GetVelocity(i) * GetVelocity(i) / 2;
            }
            return sum;
        }

        public double [] CalculateImpulseMomentum()
        {
            double x_cor = 0;
            double y_cor = 0;
            double z_cor = 0;
            double L = 0;
            for(int i =0;i<counterofpoints;i++)
            {
                x_cor += masses[i] * (y[i] * vz[i] - z[i] * vy[i]);
                y_cor += -masses[i] * (x[i] * vz[i] - z[i] * vx[i]);
                z_cor += masses[i] * (x[i] * vy[i] - y[i] * vx[i]);
            }
            L = Math.Sqrt(x_cor * x_cor + y_cor * y_cor + z_cor * z_cor);
            return new double[4] { x_cor, y_cor, z_cor ,L};
        }

        private double[] GetCenterMass()
        {
            double x_centermass = 0;
            double y_centermass = 0;
            double z_centermass = 0;
            for (int i = 0; i < counterofpoints; i++)
            {               
                 x_centermass += masses[i] * x[i];
                 y_centermass += masses[i] * y[i];
                 z_centermass += masses[i] * z[i];               
            }
            x_centermass = x_centermass / massesSigma;
            y_centermass = y_centermass / massesSigma;
            z_centermass = z_centermass / massesSigma;
            return new double[3] { x_centermass, y_centermass, z_centermass };
        }

        private double[] VelocityCenterMasses()
        {
            double[] vector = new double[3];
            double sum = 0;
            for (int i = 0; i < counterofpoints; i++)
            {
                sum += masses[i];
                vector[0] += masses[i] * vx[i];
                vector[1] += masses[i] * vy[i];
                vector[2] += masses[i] * vz[i];
            }
            vector[0] /= sum;
            vector[1] /= sum;
            vector[2] /= sum;
            return vector;
        }

        public double[] TotalImpulse()
        {
            double[] impulses = new double[4];
            for (int i = 0; i < counterofpoints; i++)
            {
                impulses[0] += masses[i] * vx[i];
                impulses[1] += masses[i] * vy[i];
                impulses[2] += masses[i] * vz[i];
            }
            impulses[3] = Math.Sqrt(Math.Pow(impulses[0], 2) + Math.Pow(impulses[1], 2) + Math.Pow(impulses[2], 2));
            return impulses;
        }

        private void MakeSystemStatic()
        {
            double[] vectorV = VelocityCenterMasses();
            double[] vectorC = GetCenterMass();
            for (int i = 0; i < counterofpoints; i++)
            {
                vx[i] -= vectorV[0];
                vy[i] -= vectorV[1];
                vz[i] -= vectorV[2];
                x[i] -= vectorC[0];
                y[i] -= vectorC[1];
                z[i] -= vectorC[2];
            }
        }


        private void ConservationEnergy()
        {
            double delta = startenergy - (TotalKineticEnergy() + TotalPotentialEnergy());
            double K = TotalKineticEnergy();
            double P = TotalPotentialEnergy();
            double coeff;
            if (1 + delta / TotalKineticEnergy() > 0)
            {
                coeff = Math.Sqrt(1 + delta / TotalKineticEnergy());
            }
            else
            {
                coeff = 1;
            }
            for (int i = 0; i < counterofpoints; i++)
            {
                vx[i] *= coeff;
                vy[i] *= coeff;
                vz[i] *= coeff;

            }
        }

        private void ConservationImpulse()
        {
            double[] VectorV = VelocityCenterMasses();
            for (int i = 0; i < counterofpoints; i++)
            {
                vx[i] -= VectorV[0];
                vy[i] -= VectorV[1];
                vz[i] -= VectorV[2];
            }
        }


        private double[] Multiplication(double[] vector, double value)
        {
            double[] ret = new double[vector.Length];
            for (int i = 0; i <= vector.Length - 1; i++)
            {
                ret[i] = value * vector[i];
            }
            return ret;
        }

        private double ScalarProduct(double[] vector1, double[] vector2)
        {
            return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
        }

        private double[] Residual(double[] vector1, double[] vector2)
        {
            return new double[3] { vector1[0] - vector2[0], vector1[1] - vector2[1], vector1[2] - vector2[2] };
        }


        private void Evaporation()
        {
            double distance = 0;
            massesSigma = 0;
            double x_centermass = 0;
            double y_centermass = 0;
            double z_centermass = 0;
            for(int i=0;i<counterofpoints;i++)
            {
                if(!evaporated.Contains(i))
                {
                    massesSigma += masses[i];
                    x_centermass = masses[i] * x[i];
                    y_centermass = masses[i] * y[i];
                    z_centermass = masses[i] * z[i];

                }
            }
            x_centermass /= massesSigma;
            y_centermass /= massesSigma;
            z_centermass /= massesSigma;

            for (int i = 0; i < this.counterofpoints; i++)
            {
                if (!evaporated.Contains(i))
                {               
                    distance = Math.Sqrt((x[i]-x_centermass) * (x[i]-x_centermass) + (y[i]-y_centermass) * (y[i]-y_centermass) + (z[i]-z_centermass)* (z[i]-z_centermass));
                    if (distance > 3 * radius)
                    {
                        if ((masses[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]) / 2 + -6.67408E-11 * masses[i] * (massesSigma-masses[i]) / distance) > 0)
                        {
                            evaporated.Add(i);

                            FileStream fs = new FileStream(@"evaporated_masses.txt", FileMode.Append, FileAccess.Write);
                            string s = Convert.ToString(totaltime) + " " + Convert.ToString(masses[i]) + "\n";
                            fs.Write(Encoding.Default.GetBytes(s), 0, s.Length);
                            fs.Close();
                        }
                    }
                }
            }
        }


        public void intgr()
        {
            if (firstcall)
            {
                Firststep();
                Write();
                firstcall = false;
            }
            if (merge)
            {
                double[] vx_ = new double[counterofpoints];
                double[] vy_ = new double[counterofpoints];
                double[] vz_ = new double[counterofpoints];
                double[] x_ = new double[counterofpoints];
                double[] y_ = new double[counterofpoints];
                double[] z_ = new double[counterofpoints];
                double[] masses_ = new double[counterofpoints];
                double[,] matrixofdistances_ = new double[counterofpoints, counterofpoints];
                double[] datastorage_ = new double[(iterations + 1) * counterofpoints * 3];
                for (int i = 0; i < counterofpoints; i++)
                {
                    x_[i] = x[i];
                    y_[i] = y[i];
                    z_[i] = z[i];
                    vx_[i] = vx[i];
                    vy_[i] = vy[i];
                    vz_[i] = vz[i];
                    masses_[i] = masses[i];
                }
                matrixofdistances = matrixofdistances_;
                x = x_;
                y = y_;
                z = z_;
                vx = vx_;
                vy = vy_;
                vz = vz_;
                masses = masses_;
                CalculateDistances();
                startenergy = TotalPotentialEnergy() + TotalKineticEnergy() + mergeenergy;
                merge = false;
                datastorage = datastorage_;

                FileStream fs = new FileStream(@"merge.txt", FileMode.Append, FileAccess.Write);
                string s = Convert.ToString(totaltime) +" " + "merge" + "\n";
                fs.Write(Encoding.Default.GetBytes(s), 0, s.Length);
                fs.Close();
            }
            if (fortran)
            {
                SphereStarCluster.intgr(ref counterofpoints, ref x[0], ref y[0], ref z[0], ref vx[0], ref vy[0], ref vz[0],
                    ref masses[0], ref matrixofdistances[0, 0], ref iterations, ref datastorage[0], ref time[0], ref startenergy, ref merge, ref this.mergeenergy);
                MakeSystemStatic();
            }
            else
            {
                SphereStarCluster.cudaIntegrate(ref counterofpoints, x, y, z, vx, vy, vz, masses, iterations, datastorage, time, startenergy, ref merge, ref mergeenergy);
                CalculateDistances();
                MakeSystemStatic();

            }
            impulse_momentum = CalculateImpulseMomentum();
            totaltime += time[time.Length - 1];
            int tempcount = evaporated.Count;
            Evaporation();
            if (tempcount != evaporated.Count)
            {
                FileStream fs = new FileStream(@"evaporation.txt", FileMode.Append, FileAccess.Write);
                string s = Convert.ToString(totaltime) + " " + Convert.ToString(counterofpoints - this.evaporated.Count) + "\n";
                fs.Write(Encoding.Default.GetBytes(s), 0, s.Length);
                fs.Close();
            }
            min_dst = MinimalDistance();
            if(GetVelocity(minimal_indexi)>GetVelocity(minimal_indexj))
            {
                max_velocity = GetVelocity(minimal_indexi);
                max_velocity_index = minimal_indexi;
            }
            else
            {
                max_velocity = GetVelocity(minimal_indexj);
                max_velocity_index = minimal_indexj;
            }
            if(boost)
            {
                vx[minimal_indexi] = -vx[minimal_indexi]; 
                vy[minimal_indexi] = -vy[minimal_indexi]; 
                vz[minimal_indexi] = -vz[minimal_indexi];
                boost = false;
            }
        }

    }
}
