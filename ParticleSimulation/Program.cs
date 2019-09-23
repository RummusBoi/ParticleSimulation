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
    

    /**
     * Sim class. Holds a list of particles, and solves their next position step by step using the barnes-hut algorith.
     * Integrator used is Leapfrog integrator.
     * 
     */
    public class Sim
    {
        //Field variables
        List<Particle> particles;

        Random rnd; //rnd is used to generate starting position and velocity of all the particles.
        Node root; //root node of the octree.
        public static Form1 form1;

        //Instead of a constructor, the 'init' function is used.
        public Sim() { return; }

        
       
        /*
         * @param form The windows form that the particles will be drawn upon.
         */
        public void init(Form1 form)
        {
            form1 = form;
            particles = new List<Particle>();
            
            rnd = new Random();

            //generates all particles with random initial parameters, using the rnd variable:
            for (int i = 0; i < SimConstants.N_PARTICLES; i++)
            {
                particles.Add(new Particle(300, 300, 300, 600.0, rnd));
            }

            //generates the first generation octree
            root = constructTree(SimConstants.SIMSIZE, SimConstants.PARTICLES_PER_BOX);
        }

        /**
         * @param timesteps The amount of timesteps that the simulation is to be run before returning.
         *
         * Giving the starting position and velocity of all particles, the simulation runs 'timesteps' amount of steps into the future
         * before returning. On the final iteration, that is, when t = timesteps - 1, all the particles are plotted on the Windows Form.
         */
        public void run(int timesteps)
        {
            //runs through as many iterations as given as input
            for (int t = 0; t < timesteps; t++)
            {
                //constructs octree and calculates center of mass for all nodes
                root = constructTree(SimConstants.SIMSIZE, SimConstants.PARTICLES_PER_BOX);
                root.calcCOMRecursive();

                //uses multiple threads to approximate the force upon each particle
                List<int> integerList = Enumerable.Range(0, SimConstants.N_PARTICLES).ToList();
                Parallel.ForEach(integerList, i =>
                {
                    calcForcesBarnesHut(particles[i], root);
                });

                //uses multiple threads to integrate every particle
                Parallel.ForEach(integerList, i =>
                {
                    particles[i].integrate();
                });
                
                //if t is the final timestep, plot each particle and the root's size to the Windows Form
                if (t == timesteps - 1)
                {
                    for (int i = 0; i < SimConstants.N_PARTICLES; i++)
                    {
                        form1.drawBox((int)root.origin[0], (int)root.origin[1], (int)root.origin[2], SimConstants.SIMSIZE);
                        form1.drawParticle((int)particles[i].pos[0], (int)particles[i].pos[1], (int)particles[i].pos[2], Color.Black);
                    }
                }
            }

            
        }

        /*
         * @param bounds The size of the simulation.
         * @param particlesPerBox The maximum amount of particles in each node.
         * @return The root node of the octree.
         *
         * Constructs the entire octree using recursion. Calls 'recursiveAddChild' on the root node for each particle
         * in the simulation. Each time a particle is added, it is added in such a way that there is no box with more
         * particles than the set amount.
         */
         
        public Node constructTree (double bounds, int particlesPerBox)
        {
            Node root = new Node(new List<double> { 0, 0, 0, }, bounds);
            foreach(Particle p in particles)
            {
                root.recursiveAddChild(p, particlesPerBox);
            }

            return root;
        }

        /*
         * @param p The particle to calculate the force upon.
         * @param root The root of the octree.
         *
         * After this function is called, 'forceRecursive' is called.
         */

        public void calcForcesBarnesHut (Particle p, Node root)
        {
            
            forceRecursive(p, root);
        }

        /*
         *    --------- COMMENTS NEEDED ------------
         */
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

                    double F = factor / pdist / pdist;
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

    /*
     * Simple class representing a particle. Containts the position, velocity and acceleration of the particle, the constructor for
     * the particle as well as the integrator function.
     */
    public class Particle
    {
        public List<double> pos;
        public List<double> vel;
        public List<double> acc;

        public double mass;

        /*
         * @param xlim The upper x coordinate limit for the particle.
         * @param ylim The upper y coordinate limit for the particle.
         * @param zlim The upper z coordinate limit for the particle.
         * @param vellim The upper velocity limit for the particle.
         * @param rnd The random variable that generates a random distribution of particles.
         *
         * xlim, ylim, and zlim constitute the coordinates that the particle can have. That is, in terms of the x-coordinate, the particle
         * can have an x-coordinate anywhere in the range [-0.5 * xlim, 0.5 * xlim].
         */
        public Particle(int xlim, int ylim, int zlim, double vellim, Random rnd)
        {
            //generates random position and velocity given constraints.
            pos = new List<double> { rnd.NextDouble() * xlim + SimConstants.SIMSIZE/2.0 - xlim/2.0, rnd.NextDouble() * ylim + SimConstants.SIMSIZE / 2.0 - ylim/2.0, rnd.NextDouble() * zlim + SimConstants.SIMSIZE / 2.0 - zlim / 2.0};
            
            vel = new List<double> { rnd.NextDouble() * vellim*2 - vellim, rnd.NextDouble() * vellim*2 - vellim, rnd.NextDouble() * vellim*2 - vellim };
            acc = new List<double> { 0, 0, 0 };

            mass = SimConstants.PARTICLE_MASS;
        }

        //Leapfrog integrator function
        public void integrate()
        {
            for (int i = 0; i < 3; i++)
            {
                vel[i] = vel[i] + acc[i] * SimConstants.DT;
                pos[i] = pos[i] + vel[i] * SimConstants.DT;
                
                acc[i] = 0;
            }

            correctOutOfBounds();
        }

        private void correctOutOfBounds()
        {
            if (pos[0] > SimConstants.SIMSIZE || pos[0] < 0)
            {
                vel[0] = -vel[0];
            }

            if (pos[1] > SimConstants.SIMSIZE || pos[1] < 0)
            {
                vel[1] = -vel[1];
            }

            if (pos[2] > SimConstants.SIMSIZE || pos[2] < 0)
            {
                vel[2] = -vel[2];
            }
        }
    }

    /*
     * Node class that represents each octant in the octree. Each Node has 8 of children, each of which subdivide
     * the parent's spatial extent into 8 smaller regions. Each Node also has a list of particles. This list is only relevant if
     * the Node is a leaf Node with no Node children. In this case, the particle list containts all the particles within the
     * spatial region of the Node.
     */

    public class Node
    {
        //Lists of children and particles.
        public List<Node> children;
        public List<Particle> particles;

        
        public double size; //The size of this node. Since each node is a cube, this is simply the length of each side of the cube.
        public int nparticles; //The amount of particles contained in this node.

        public List<double> com; //Coordinates to the center of mass of the node.
        public List<double> origin; //Coordinates to the origin of the node.
        public double mass; //The sum of all the particles inside the node.

        /*
         * @param origin The origin of this node. The origin of a node is the corner closest to the origin. For example, the spatial extent
         * of the root node of the octree goes from (0, 0, 0) to (simsize, simsize, simsize), where simsize = SimConstants.SIMSIZE. In this
         * case, the origin would be (0, 0, 0), since this is the corner closest to the origin.
         */
        public Node(List<double> origin, double size) {
            //Sets field variables of the node
            com = new List<double> { 0.0, 0.0, 0.0 };
            this.size = size;
            mass = 0;
            this.origin = new List<double> { origin[0], origin[1], origin[2] };
            particles = new List<Particle>();
            children = new List<Node>();
            
        }

        public void addChild(Node child) { children.Add(child);  }

        public Node getChild(int index) { return children[index]; }

        /*
         * Recursive function for calculating the center of mass of this node. Fundamentally it works by finding the COM of all
         * children of the node, and then calculating a weighted average of these to compute the COM of this node.
         */
        public void calcCOMRecursive ()
        {
            //Calculates the COM of all the children of the node.
            foreach (Node child in children)
            {
                child.calcCOMRecursive();
            }

            //Once children COMs have been calculated, the COM of this node is calculated.
            
            if(children.Count + particles.Count > 0) 
            {
                calcCOM(); //calculates COM of this node. 
            } else //in the case where this node has no children and no particles within it, the COM is simply set to (0, 0, 0) to avoid
                    //arithmetic errors.
            {
                com[0] = 0;
                com[1] = 0;
                com[2] = 0;
                mass = 0;
            }
        }

        /*
         * Calculates the COM of this node. This function assumes that the COMs of all children of the nodes have already been
         * calculated. 
         */
        public void calcCOM ()
        {
            //Initiates the COM to (0, 0, 0) and total mass to 0.
            com[0] = 0;
            com[1] = 0;
            com[2] = 0;

            mass = 0;

            /*
             * Here we both loop through all children and all particles, however in practice only one of these lists will be iterated through.
             * In the case where the node is not a leaf node, it will have 8 children and 0 particles. In the case where the node is a
             * leaf node, the node will have 0 children and a few particles. We could have made a check here to see whether the node is
             * a leaf or not, but my implementation works just as well.
             */

            //Computes the COM of all the children COMs. 
            foreach (Node node in children)
            {
                mass += node.mass;
                com[0] += node.com[0] * node.mass;
                com[1] += node.com[1] * node.mass;
                com[2] += node.com[2] * node.mass;
            }

            //Computes the COM of all the children particles.
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

        /*
         * @param p The particle to be added to the node.
         * @param particlesPerBox The maximum amount of particles in every node.
         *
         * Recursive function that adds the particle p to the octree. If there are already enough particles in this node, the node is
         * split into 8 subregions, and all the particles of this node gets added to the correct child.
         */
        public void recursiveAddChild(Particle p, int particlesPerBox)
        {
            //Starts off by adding the particle p to this node.
            particles.Add(p);

            //Checks if the condition of maximum particles per box still holds.
            if (particles.Count > particlesPerBox && children.Count == 0)
            {
                //If it doesn't hold, and the node does not yet have any children, then 8 children are added.
                children.Add(new Node(new List<double> { origin[0], origin[1], origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0], origin[1], origin[2] + size / 2 }, size / 2));
                children.Add(new Node(new List<double> { origin[0], origin[1] + size / 2, origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0], origin[1] + size / 2, origin[2] + size / 2 }, size / 2));
                
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1], origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1], origin[2] + size / 2 }, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1] + size / 2, origin[2]}, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1] + size / 2, origin[2] + size / 2 }, size / 2));
                
                //Goes through each particle of the node and calls recursiveAddChild on the child that contains the particle.
                foreach (Particle p1 in particles)
                {
                    int index = 0;
                    if(p1.pos[0] > origin[0] + size/2.0) { index += 4; }
                    if(p1.pos[1] > origin[1] + size/2.0) { index += 2; }
                    if(p1.pos[2] > origin[2] + size/2.0) { index += 1; }
                    children[index].recursiveAddChild(p1, particlesPerBox);
                }
                //Since all the particles are now distributed among the children of the node, the particle list of this node is cleared.
                particles.Clear();

            } else if (children.Count != 0) //If the condition does hold, but the 8 children have already been added, then recursiveAddChild
                                            //is simply called on the child that contains the particle.
            {
                int index = 0;
                if (p.pos[0] > origin[0] + size / 2.0) { index += 4; }
                if (p.pos[1] > origin[1] + size / 2.0) { index += 2; }
                if (p.pos[2] > origin[2] + size / 2.0) { index += 1; }
                children[index].recursiveAddChild(p, particlesPerBox);

                //Clears the particle list.
                particles.Clear();
            }

            //If the condition still holds, that is there aren't too many particles in this node, then the particle is simply added
            //to this node, and no further action is needed.
        }
    }
}