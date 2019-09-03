using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Drawing;
using System.Timers;
using System.Threading;


namespace ParticleSimulation
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Form1 form1 = new Form1();
            Application.Run(form1);
        }
    }
    


    public class Sim
    {
        List<Particle> particles;
        List<Particle> particles2;
        double dt;
        Random rnd;
        Node root;
        int n;
        public static int simSize;
        public static Form1 form1;

        double barnesHutDelta = 0;

        public Sim() { return; }

        public void init(int n, Form1 form)
        {
            Sim.form1 = form;
            particles = new List<Particle>();

            this.n = n;

            Sim.simSize = 2000;
            dt = 0.0001;
            rnd = new Random();

            for (int i = 0; i < n; i++)
            {
                particles.Add(new Particle(300, 300, 300, 600.0, rnd, i));
            }
            

            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    Particle.calcForces(particles[i], particles[j]);
                }
                particles[i].acc[0] = 0;
                particles[i].acc[1] = 0;
                particles[i].acc[2] = 0;
            }
            
            particles2 = new List<Particle>(particles.Count);

            foreach (Particle particle in particles)
            {
                Particle clone = new Particle(0, 0, 0, 0, rnd, particle.ID);
                clone.pos = new List<double> { particle.pos[0], particle.pos[1], particle.pos[2] };
                clone.vel = new List<double> { particle.vel[0], particle.vel[1], particle.vel[2] };
                clone.acc = new List<double> { particle.acc[0], particle.acc[1], particle.acc[2] };
                clone.mass = particle.mass;
                particles2.Add(clone);
            }
            
            root = constructTree(Sim.simSize, 2);
        }

        public void threadRun (object parameters)
        {

            List<int> p = (List<int>)parameters;
            int start = p[0];
            int end = p[1];
            int id = p[2];

            while(true)
            {
                if (threadsrun[id] == true)
                {
                    for (int i = start; i < end; i++)
                    {
                        calcForcesBarnesHut(particles[i], root);
                    }
                    threadsrun[id] = false;
                }
            }

            
        }

        public void run(int timesteps)
        {
            bool compare = false;
            bool barnes = true;
            bool threading = true;

            //compare is not yet enabled for threading

            bool plot = true;
            for (int t = 0; t < timesteps; t++)
            {

                if (threading)
                {
                    if (barnes)
                    {
                        List<int> integerList = Enumerable.Range(0, n).ToList();
                        root = constructTree(Sim.simSize, 1);
                        root.calcCOMRecursive();
                        Parallel.ForEach(integerList, i =>
                        {
                            calcForcesBarnesHut(particles[i], root);
                        });
                        Parallel.ForEach(integerList, i =>
                        {
                            particles[i].integrate(dt);
                        });
                    } else
                    {
                        List<int> integerList = Enumerable.Range(0, n).ToList();
                        Parallel.ForEach(integerList, i => {
                            for (int j = i + 1; j < n; j++)
                            {
                                Particle.calcForces(particles[i], particles[j]);
                            }
                        });
                        Parallel.ForEach(integerList, i =>
                        {
                            particles[i].integrate(dt);
                        });
                    }


                    if (t == timesteps - 1 && plot)

                    {
                        for (int i = 0; i < n; i++)
                        {
                            //root.plotCOMRecursive();
                            Sim.form1.drawBox((int)root.origin[0], (int)root.origin[1], (int)root.origin[2], Sim.simSize);
                            form1.drawParticle((int)particles[i].pos[0], (int)particles[i].pos[1], (int)particles[i].pos[2], Color.Black);
                        }
                    }
                }

                else
                {

                    root = constructTree(Sim.simSize, 1);
                    root.calcCOMRecursive();
                    for (int i = 0; i < n; i++)
                    {
                        if (barnes)
                        {
                            calcForcesBarnesHut(particles[i], root);
                            if (compare)
                            {
                                for (int j = i + 1; j < n; j++)
                                {
                                    Particle.calcForces(particles2[i], particles2[j]);
                                }
                            }

                            if (t == timesteps - 1)
                            {
                                //root.plotCOMRecursive();
                                Sim.form1.drawBox((int)root.origin[0], (int)root.origin[1], (int)root.origin[2], Sim.simSize);
                                form1.drawParticle((int)particles[i].pos[0], (int)particles[i].pos[1], (int)particles[i].pos[2], Color.Black);
                                if (compare)
                                {
                                    form1.drawParticle((int)particles2[i].pos[0], (int)particles2[i].pos[1], (int)particles2[i].pos[2], Color.Blue);
                                }
                            }
                        }
                        else
                        {
                            for (int j = i + 1; j < n; j++)
                            {
                                Particle.calcForces(particles[i], particles[j]);
                            }
                            if (t == timesteps - 1)
                            {
                                Sim.form1.drawBox(0, 0, 0, Sim.simSize);
                                form1.drawParticle((int)particles[i].pos[0], (int)particles[i].pos[1], (int)particles[i].pos[2], Color.Black);
                            }
                        }

                    }
                    for (int i = 0; i < n; i++)
                    {
                        particles[i].integrate(dt);
                    }
                    if (compare)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            particles2[i].integrate(dt);
                        }
                    }
                }
            }
        }

        public Node constructTree (double bounds, int particlesPerBox)
        {
            Node.ndivisions = 0;
            Node root = new Node(new List<double> { 0, 0, 0, }, bounds);
            foreach(Particle p in particles)
            {
                root.recursiveAddChild(p, particlesPerBox);
            }

            return root;
        }

        public void calcForcesBarnesHut (Particle p, Node root)
        {
            forceRecursive(p, root);
        }

        private void forceRecursive (Particle p, Node n)
        {
            double dx = n.com[0] - p.pos[0];
            double dy = n.com[1] - p.pos[1];
            double dz = n.com[2] - p.pos[2];

            double dist = Math.Pow(dx*dx + dy*dy + dz*dz, 0.5) + 1;
            

            bool cond = n.size / dist < barnesHutDelta;
            //if cond is true, the node is sufficiently far away

            if (cond)
            {
                double F = 6.67 * Math.Pow(10, -11) * n.mass * p.mass / dist / dist;
                double fx = F * dx / dist;
                double fy = F * dy / dist;
                double fz = F * dz / dist;

                p.acc[0] += fx / p.mass;
                p.acc[1] += fy / p.mass;
                p.acc[2] += fz / p.mass;
            } else {
                foreach (Node child in n.children)
                {
                    forceRecursive(p, child);
                }
            }
            
            if(n.children.Count == 0 && cond == false)
            {
                foreach (Particle p1 in n.particles)
                {

                    double pdx = p1.pos[0] - p.pos[0];
                    double pdy = p1.pos[1] - p.pos[1];
                    double pdz = p1.pos[2] - p.pos[2];

                    double distCorrection = 1;

                    double pdist = Math.Pow(pdx*pdx + pdy*pdy + pdz*pdz, 0.5) + distCorrection;

                    if(pdist == distCorrection) { break; }

                    double factor = 6.67 * Math.Pow(10, -11) * p1.mass * p.mass;

                    double F = factor / pdist / pdist - factor / (dist*dist*dist*dist*dist* dist * dist);
                    double fx = F * pdx / pdist;
                    double fy = F * pdy / pdist;
                    double fz = F * pdz / pdist;

                    p.acc[0] += fx / p.mass;
                    p.acc[1] += fy / p.mass;
                    p.acc[2] += fz / p.mass;
                }
            }
        }
    }

    public class Particle
    {
        public List<double> pos { get; set; }
        public List<double> vel;
        public List<double> acc;
        public int ID;

        public double mass { get; set; }

        public Particle(int xlim, int ylim, int zlim, double vellim, Random rnd, int ID)
        {
            this.ID = ID;
            pos = new List<double> { rnd.NextDouble() * xlim + Sim.simSize/2.0 - xlim/2.0, rnd.NextDouble() * ylim + Sim.simSize / 2.0 - ylim/2.0, rnd.NextDouble() * zlim + Sim.simSize / 2.0 - zlim / 2.0};
            
            vel = new List<double> { rnd.NextDouble() * vellim*2 - vellim, rnd.NextDouble() * vellim*2 - vellim, rnd.NextDouble() * vellim*2 - vellim };
            acc = new List<double> { 0, 0, 0 };
            mass = Math.Pow(10, 18);
        }

        public void calcFirstVel(double dt)
        {
            for (int i = 0; i < 3; i++)
            {
                vel[i] = vel[i] + acc[i] * dt / 2;
            }
        }

        public void integrate(double dt)
        {
            for (int i = 0; i < 3; i++)
            {
                vel[i] = vel[i] + acc[i] * dt;
                pos[i] = pos[i] + vel[i] * dt;
                
                acc[i] = 0;
            }
        }

        public static void calcForces(Particle p1, Particle p2)
        {

            double G = 6.67 * Math.Pow(10, -11);
            double dx = p2.pos[0] - p1.pos[0];
            double dy = p2.pos[1] - p1.pos[1];
            double dz = p2.pos[2] - p1.pos[2];


            double dist = Math.Pow(dx*dx + dy*dy + dz*dz, 0.5) + 1;
            double F = G * p1.mass * p2.mass / dist / dist;
            double fx = F * dx / dist;
            double fy = F * dy / dist;
            double fz = F * dz / dist;


            p1.acc[0] = p1.acc[0] + fx / p1.mass;
            p1.acc[1] = p1.acc[1] + fy / p1.mass;
            p1.acc[2] = p1.acc[2] + fz / p1.mass;

            p2.acc[0] = p2.acc[0] - fx / p2.mass;
            p2.acc[1] = p2.acc[1] - fy / p2.mass;
            p2.acc[2] = p2.acc[2] - fz / p2.mass;


        }
    }

    public class Node
    {
        public List<Node> children;
        public double size;
        public int nparticles;
        public static int ndivisions;
        public List<Particle> particles; //the particles that this region points to, ie. NOT all the particles contained within it, only the ones it points to
                                    //in the node tree.
        
            //origin is the upper left corner, size is the size of each size in the square.
        public List<double> com;
        public List<double> origin;
        public double mass;
        public Node(List<double> origin, double size) {
            com = new List<double> { 0.0, 0.0, 0.0 };
            this.size = size; mass = 0;
            this.origin = new List<double> { origin[0], origin[1], origin[2] };
            particles = new List<Particle>();
            children = new List<Node>();

            //Sim.form1.drawBox((int)origin[0], (int)origin[1], (int)origin[2], (int)this.size);
            
        }
        public void addChild(Node child) { children.Add(child);  }
        public Node getChild(int index) { return children[index]; }

        public void printInformation (int tabs)
        {
            string printString = "";
            for (int i = 0; i < tabs; i++) { printString += "   "; }
            printString += "Origin at " + origin[0] + " : " + origin[1] + " : " + origin[2] + "; size = " + size;
        

            Console.WriteLine(printString);

            if(particles.Count == 0)
            {
                
                if (children.Count != 0)
                {
                    foreach (Node child in children) {
                        child.printInformation(tabs + 1);
                    }
                }
            } else
            {
                foreach (Particle p in particles)
                {
                    string print = "";
                    for (int i = 0; i < tabs; i++) { print += "   "; }
                    print += "-";
                    print += "Particle at " + p.pos[0] + " : " + p.pos[1] + " : " + p.pos[2];
                    Console.WriteLine(print);
                }
            }
        }

        public void plotCOMRecursive ()
        {
            if (particles.Count == 0 && !(com[0]==0 && com[1] == 0 && com[2] == 0)) {
                Sim.form1.drawParticle((int)com[0], (int)com[1], (int)com[2], Color.Red);
            }
            foreach (Node child in children)
            {
                child.plotCOMRecursive();
            }
        }

        public void calcCOMRecursive ()
        {
            foreach (Node child in children)
            {
                child.calcCOMRecursive();
            }
            
            if(children.Count + particles.Count > 0)
            {
                calcCOM();
            } else
            {
                com[0] = 0;
                com[1] = 0;
                com[2] = 0;
                mass = 0;
            }
        }

        public void calcCOM ()
        {
            com[0] = 0;
            com[1] = 0;
            com[2] = 0;

            mass = 0;

            foreach (Node node in children)
            {
                mass += node.mass;
                com[0] += node.com[0] * node.mass;
                com[1] += node.com[1] * node.mass;
                com[2] += node.com[2] * node.mass;
            }

            foreach(Particle p in particles)              
            {
                mass += p.mass;
                com[0] += p.pos[0] * p.mass;
                com[1] += p.pos[1] * p.mass;
                com[2] += p.pos[2] * p.mass;
            }

            com[0] /= mass;
            com[1] /= mass;
            com[2] /= mass;
        }

        public void recursiveAddChild(Particle p, int particlesPerBox)
        {
            particles.Add(p);

            if (particles.Count > particlesPerBox && children.Count == 0)
            {
                
                children.Add(new Node(new List<double> { origin[0], origin[1], origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0], origin[1], origin[2] + size / 2 }, size / 2));
                children.Add(new Node(new List<double> { origin[0], origin[1] + size / 2, origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0], origin[1] + size / 2, origin[2] + size / 2 }, size / 2));
                
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1], origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1], origin[2] + size / 2 }, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1] + size / 2, origin[2]}, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1] + size / 2, origin[2] + size / 2 }, size / 2));
                
                Node.ndivisions += 8;

                foreach (Particle p1 in particles)
                {
                    int index = 0;
                    if(p1.pos[0] > origin[0] + size/2.0) { index += 4; }
                    if(p1.pos[1] > origin[1] + size/2.0) { index += 2; }
                    if(p1.pos[2] > origin[2] + size/2.0) { index += 1; }
                    children[index].recursiveAddChild(p1, particlesPerBox);
                }
                particles.Clear();

            } else if (children.Count != 0)
            {
                int index = 0;
                if (p.pos[0] > origin[0] + size / 2.0) { index += 4; }
                if (p.pos[1] > origin[1] + size / 2.0) { index += 2; }
                if (p.pos[2] > origin[2] + size / 2.0) { index += 1; }
                children[index].recursiveAddChild(p, particlesPerBox);
                particles.Clear();
            }
        }
    }
}