using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;



namespace ParticleSimulation
{
    public partial class Form1 : Form
    {
        int x, y;
        double cameraDist, rotation;

        static double timeElapsed;
        static int timeStep;

        List<Tuple<int, int, Color>> particlePlots;

        public static int maxX, minX, maxY, minY, maxZ, minZ;

        List<int> boxX;
        List<int> boxY;
        List<int> boxSize;
        List<int> boxZ;
        List<int> xcom;
        List<int> ycom;

        bool firstGen;
        int n;
        Sim sim;
        Bitmap bm;
        Graphics gbm;
        Stopwatch sw;
        

        public Form1()
        {
            Form1.timeElapsed = 0;
            Form1.timeStep = 0;
            InitializeComponent();
            Size = new System.Drawing.Size(1000, 1000);
            x = 200;
            y = 100;
            cameraDist = 10000;
            rotation = 0;

            particlePlots = new List<Tuple<int, int, Color>>();

            xcom = new List<int>();
            ycom = new List<int>();
            

            boxX = new List<int>();
            boxY = new List<int>();
            boxZ = new List<int>();
            boxSize = new List<int>();

            /*
            // ------  RUN TEST ------ //
            List<int> parameters = new List<int> { 50, 2000, 70 };
            Thread th = new Thread(new ParameterizedThreadStart(runTest), 10000000);
            th.Start(parameters);

            th.Join();

            // ------  RUN TEST ------ //
            */
            int nparticles = 250;
            firstGen = true;
            
            sim = new Sim();
            sim.init(nparticles, this);

            bm = new Bitmap(1000, 1000);
            gbm = Graphics.FromImage(bm);

            n = 0;
            sw = new Stopwatch();
            sw.Start();
        }

        private void runTest (object parameters)
        {
            List<int> p = (List<int>)parameters;
            int start = p[0];
            int end = p[1];
            int interval = p[2];
            
            sw = new Stopwatch();

            List<int> np = new List<int>();
            List<int> times = new List<int>();

            for (int i = start; i <= end; i += interval) {
                sw.Reset();
                sw.Start();
                sim = new Sim();
                sim.init(i, this);
                sim.run(150);
                sw.Stop();
                np.Add(i);
                times.Add((int)sw.ElapsedMilliseconds);
                Console.WriteLine("Finished " + i + " particles in " + sw.ElapsedMilliseconds + " milliseconds.");
            }

            Console.WriteLine("Particles:");
            np.ForEach(o => Console.WriteLine(o));
            Console.WriteLine("Times:");
            times.ForEach(o => Console.WriteLine(o));
        }

        private void Form1_Paint_1(object sender, PaintEventArgs e)
        {
            DoubleBuffered = true;
            bool drawLine = false;
            Rectangle ellipseBounds = new Rectangle(x, y, 3, 3); //like in your code sample

            int steps = 10;
            sim.run(steps);
            
            //gbm.Clear(Color.White);

            using (SolidBrush brush = new SolidBrush(Color.Black))
            {

                Pen pen = new Pen(brush);
                if (!drawLine)
                {
                    gbm.Clear(Color.White);
                    for (int i = 0; i < particlePlots.Count; i++)
                    {
                        ellipseBounds.X = particlePlots[i].Item1;
                        ellipseBounds.Y = particlePlots[i].Item2;
                        brush.Color = particlePlots[i].Item3;
                        gbm.FillEllipse(brush, ellipseBounds); //for example, do it before drawing lines.
                    }
                }
                else
                {
                    /*
                    for (int i = 0; i < xpos.Count; i++)
                    {
                        Point p1, p2;
                        Console.WriteLine("bitch;");
                        if (firstGen)
                        {
                            p1 = new Point(xpos[i], ypos[i]);
                        }
                        else
                        {
                            p1 = new Point(lastxpos[i], lastypos[i]);
                        }

                        p2 = new Point(xpos[i], ypos[i]);
                        gbm.DrawLine(pen, p1, p2);

                        lastxpos[i] = xpos[i];
                        lastypos[i] = ypos[i];
                    }
                    */
                }


                particlePlots.Clear();
                firstGen = false;
            }

            using (Brush brush = new SolidBrush(Color.Blue))
            {
                for (int i = 0; i < boxX.Count; i++)
                {
                    Pen pen = new Pen(brush);

                    double p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y, p6x, p6y, p7x, p7y, p8x, p8y;
                    int x = boxX[i], y = boxY[i], z = boxZ[i], s = boxSize[i];

                    
                    projectPoint(x, y, z, rotation, out p1x, out p1y);
                    projectPoint(x + s, y, z, rotation, out p2x, out p2y);
                    projectPoint(x + s, y + s, z, rotation, out p3x, out p3y);
                    projectPoint(x, y + s, z, rotation, out p4x, out p4y);
                    projectPoint(x, y, z + s, rotation, out p5x, out p5y);
                    projectPoint(x + s, y, z + s, rotation, out p6x, out p6y);
                    projectPoint(x + s, y + s, z + s, rotation, out p7x, out p7y);
                    projectPoint(x, y + s, z + s, rotation, out p8x, out p8y);

                    gbm.DrawLine(pen, (float)p1x, (float)p1y, (float)p2x, (float)p2y);
                    gbm.DrawLine(pen, (float)p2x, (float)p2y, (float)p3x, (float)p3y);
                    gbm.DrawLine(pen, (float)p3x, (float)p3y, (float)p4x, (float)p4y);
                    gbm.DrawLine(pen, (float)p4x, (float)p4y, (float)p1x, (float)p1y);
                    
                    gbm.DrawLine(pen, (float)p5x, (float)p5y, (float)p6x, (float)p6y);
                    gbm.DrawLine(pen, (float)p6x, (float)p6y, (float)p7x, (float)p7y);
                    gbm.DrawLine(pen, (float)p7x, (float)p7y, (float)p8x, (float)p8y);
                    gbm.DrawLine(pen, (float)p8x, (float)p8y, (float)p5x, (float)p5y);
                    
                    gbm.DrawLine(pen, (float)p1x, (float)p1y, (float)p5x, (float)p5y);
                    gbm.DrawLine(pen, (float)p2x, (float)p2y, (float)p6x, (float)p6y);
                    gbm.DrawLine(pen, (float)p3x, (float)p3y, (float)p7x, (float)p7y);
                    gbm.DrawLine(pen, (float)p4x, (float)p4y, (float)p8x, (float)p8y);
                    
                }
                boxX.Clear();
                boxY.Clear();
                boxZ.Clear();
                boxSize.Clear();

                for (int i = 0; i < xcom.Count; i++)
                {
                    Pen pen = new Pen(Color.Red);
                    gbm.DrawEllipse(pen, xcom[i], ycom[i], 2, 2);
                }
                xcom.Clear();
                ycom.Clear();
            }

            e.Graphics.DrawImage(bm, 0, 0);
            
            Form1.timeElapsed += sw.ElapsedMilliseconds;
            Form1.timeStep += 1;
            rotation += 0.005;
        }

        public void drawParticle (int x, int y, int z, Color color)
        {
            //add projected particle to drawing buffer
            double xp, yp;
            projectPoint((double)x, (double)y, (double)z, rotation, out xp, out yp);
            
            Tuple<int, int, Color> plot = new Tuple<int, int, Color>((int)xp, (int)yp, color);
            particlePlots.Add(plot);
        }

        public void drawCOM (int x, int y, int z)
        {
            double xp, yp;
            projectPoint((double)x, (double)y, (double)z, rotation, out xp, out yp);
            xcom.Add((int)xp);
            ycom.Add((int)yp);
        }

        public void drawBox(int x, int y, int z, int size)
        {
            //add points to buffer
            boxX.Add(x);
            boxY.Add(y);
            boxZ.Add(z);
            boxSize.Add(size);
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            x += 1;
            Invalidate();
        }

        private void projectPoint(double x, double y, double z, double rotation, out double xp, out double yp)
        {
            double screenOrigoDist = cameraDist / 2;
            double shiftX = 500;
            double shiftY = 300;

            x -= Sim.simSize / 2.0;
            y -= Sim.simSize / 2.0;
            z -= Sim.simSize / 2.0;

            double xr = Math.Cos(rotation) * x - Math.Sin(rotation) * z;
            double zr = Math.Sin(rotation) * x + Math.Cos(rotation) * z;

            x += Sim.simSize / 2.0;
            y += Sim.simSize / 2.0;
            z += Sim.simSize / 2.0;

            xp = ((zr + screenOrigoDist) / (zr + cameraDist + screenOrigoDist)) * xr + shiftX;
            yp = ((zr + screenOrigoDist) / (zr + cameraDist + screenOrigoDist)) * y + shiftY;
        }
    }
}