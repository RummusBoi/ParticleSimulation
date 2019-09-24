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
        double rotation;
        List<Tuple<int, int, Color>> particlePlots;

        public static int maxX, minX, maxY, minY, maxZ, minZ;

        List<(int X, int Y, int Z, int Size)> boxCoords;
        List<(int X, int Y)> com;
        

        Sim sim;
        Bitmap bm;
        Graphics gbm;
        

        public Form1()
        {
            InitializeComponent();
            Size = new System.Drawing.Size(SimConstants.IMG_WIDTH, SimConstants.IMG_HEIGHT);
            rotation = 0;

            particlePlots = new List<Tuple<int, int, Color>>();

            com = new List<(int, int)>();

            boxCoords = new List<(int, int, int, int)>();
            
            sim = new Sim();
            sim.init(this, new BarnesHutForceSolver());

            bm = new Bitmap(SimConstants.IMG_WIDTH, SimConstants.IMG_HEIGHT);
            gbm = Graphics.FromImage(bm);
        }

        private void Form1_Paint_1(object sender, PaintEventArgs e)
        {
            DoubleBuffered = true;
            Rectangle ellipseBounds = new Rectangle(0, 0, 1, 1);

            
            sim.run(SimConstants.STEPS_PER_FRAME);

            using (SolidBrush brush = new SolidBrush(Color.Black))
            {
                Pen pen = new Pen(brush);
                gbm.Clear(Color.White);
                for (int i = 0; i < particlePlots.Count; i++)
                {
                    ellipseBounds.X = particlePlots[i].Item1;
                    ellipseBounds.Y = particlePlots[i].Item2;
                    brush.Color = particlePlots[i].Item3;

                    gbm.FillRectangle(brush, ellipseBounds);
                }

                particlePlots.Clear();
            }

            using (Brush brush = new SolidBrush(Color.Blue))
            {
                for (int i = 0; i < boxCoords.Count; i++)
                {
                    Pen pen = new Pen(brush);

                    double p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y, p5x, p5y, p6x, p6y, p7x, p7y, p8x, p8y;
                    int x = boxCoords[i].X, y = boxCoords[i].Y, z = boxCoords[i].Z, s = boxCoords[i].Size;

                    
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

                boxCoords.Clear();

                for (int i = 0; i < com.Count; i++)
                {
                    Pen pen = new Pen(Color.Red);
                    gbm.DrawEllipse(pen, com[i].X, com[i].Y, 2, 2);
                }
                com.Clear();
            }

            e.Graphics.DrawImage(bm, 0, 0);

            rotation += SimConstants.ROTATION_SPEED;
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
            projectPoint(x, y, z, rotation, out xp, out yp);
            com.Add(((int)xp, (int)yp));
        }

        public void drawBox(int x, int y, int z, int size)
        {
            //add points to buffer
            boxCoords.Add((x, y, z, size));
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            Invalidate();
        }

        private void projectPoint(double x, double y, double z, double rotation, out double xp, out double yp)
        {
            double cameraDist = 7000;
            double camToScreenDist = 1000;

            double cameraX = SimConstants.SIMSIZE / 2;
            double cameraY = SimConstants.SIMSIZE / 2;
            double cameraZ = -cameraDist;
            
            x -= SimConstants.SIMSIZE / 2;
            y -= SimConstants.SIMSIZE / 2;
            z -= SimConstants.SIMSIZE / 2;
            
            double tmpx = Math.Cos(rotation) * x + Math.Sin(rotation) * z;
            z = -Math.Sin(rotation) * x + Math.Cos(rotation) * z;
            x = tmpx;
            
            x += SimConstants.SIMSIZE / 2;
            y += SimConstants.SIMSIZE / 2;
            z += SimConstants.SIMSIZE / 2;
            
            //update relative position of particle:
            x -= cameraX;
            y -= cameraY;
            z -= cameraZ;


            xp = (camToScreenDist) / z * x + SimConstants.IMG_WIDTH / 2;
            yp = (camToScreenDist) / z * y + SimConstants.IMG_HEIGHT / 2;
        }
    }
}