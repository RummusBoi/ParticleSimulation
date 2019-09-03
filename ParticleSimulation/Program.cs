using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;



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
        Random rnd;
        Node root;
        public static int simSize;
        public static Form1 form1;

        public Sim() { return; }

        public void init(Form1 form)
        {
            Sim.form1 = form;
            particles = new List<Particle>();
            
            rnd = new Random();

            for (int i = 0; i < SimConstants.N_PARTICLES; i++)
            {
                particles.Add(new Particle(300, 300, 300, 600.0, rnd));
            }
            
            root = constructTree(SimConstants.SIMSIZE, 1);
        }

        public void run(int timesteps)
        {
            for (int t = 0; t < timesteps; t++)
            {
                List<int> integerList = Enumerable.Range(0, SimConstants.N_PARTICLES).ToList();
                root = constructTree(Sim.simSize, 1);
                root.calcCOMRecursive();
                Parallel.ForEach(integerList, i =>
                {
                    calcForcesBarnesHut(particles[i], root);
                });
                Parallel.ForEach(integerList, i =>
                {
                    particles[i].integrate();
                });

                if (t == timesteps - 1)
                {
                    for (int i = 0; i < SimConstants.N_PARTICLES; i++)
                    {
                        Sim.form1.drawBox((int)root.origin[0], (int)root.origin[1], (int)root.origin[2], Sim.simSize);
                        form1.drawParticle((int)particles[i].pos[0], (int)particles[i].pos[1], (int)particles[i].pos[2], Color.Black);
                    }
                }
            }

            
        }

        public Node constructTree (double bounds, int particlesPerBox)
        {
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
            

            bool cond = n.size / dist < SimConstants.BARNES_HUT_DELTA;
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
        public List<double> pos;
        public List<double> vel;
        public List<double> acc;

        public double mass;

        public Particle(int xlim, int ylim, int zlim, double vellim, Random rnd)
        {
            pos = new List<double> { rnd.NextDouble() * xlim + Sim.simSize/2.0 - xlim/2.0, rnd.NextDouble() * ylim + Sim.simSize / 2.0 - ylim/2.0, rnd.NextDouble() * zlim + Sim.simSize / 2.0 - zlim / 2.0};
            
            vel = new List<double> { rnd.NextDouble() * vellim*2 - vellim, rnd.NextDouble() * vellim*2 - vellim, rnd.NextDouble() * vellim*2 - vellim };
            acc = new List<double> { 0, 0, 0 };
            mass = Math.Pow(10, 18);
        }

        public void integrate()
        {
            for (int i = 0; i < 3; i++)
            {
                vel[i] = vel[i] + acc[i] * SimConstants.DT;
                pos[i] = pos[i] + vel[i] * SimConstants.DT;
                
                acc[i] = 0;
            }
        }
    }

    public class Node
    {
        public List<Node> children;
        public double size;
        public int nparticles;
        public List<Particle> particles; 
        //the particles that this region points to, ie. NOT all the particles contained within it, only the ones it points to
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
            
        }

        public void addChild(Node child) { children.Add(child);  }

        public Node getChild(int index) { return children[index]; }

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