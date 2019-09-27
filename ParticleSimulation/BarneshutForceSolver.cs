using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ParticleSimulation
{    
    class BarnesHutForceSolver : ForceSolver
    {
        public BarnesHutForceSolver() { }

        /*
         * @param particles the list of particles that should be integrated.
         * @param timestepSize the size of each timestep.
         * @param bounds the size of the simulation
         */

        public void integrateAll(List<Particle> particles, double timestepSize, double bounds)
        {
            List<int> integerList = Enumerable.Range(0, particles.Count).ToList();
            Parallel.ForEach(integerList, i =>
            {
                integrate(particles[i], timestepSize, bounds);
            });          
        }

		/*
         * @param p the particle to be integrated.
         * @param timestepSize the size of each timestep.
         * @param bounds the size of the simulation.
         */

		private void integrate(Particle p, double timestepSize, double bounds)
        {
            for (int i = 0; i < 3; i++)
            {
                p.vel[i] = p.vel[i] + p.acc[i] * timestepSize;
                p.pos[i] = p.pos[i] + p.vel[i] * timestepSize;

                p.acc[i] = 0;
            }

            correctOutOfBounds(p, bounds);
        }

        /*
         * @param p the particle whose position is to be checked.
         * @param bounds the size of the simulation. 
         * Function that determines what happens when a particle goes out of bounds.
         * In this case, the particles simply bounce off the edges of the box.
         */

        private void correctOutOfBounds(Particle p, double bounds)
        {
            if (p.pos[0] > SimConstants.SIMSIZE || p.pos[0] < 0)
            {
                p.vel[0] = -p.vel[0]*1;
            }

            if (p.pos[1] > SimConstants.SIMSIZE || p.pos[1] < 0)
            {
                p.vel[1] = -p.vel[1] * 1;
            }

            if (p.pos[2] > SimConstants.SIMSIZE || p.pos[2] < 0)
            {
                p.vel[2] = -p.vel[2] * 1;
            }
        }

        /*
         * @param p The particle to calculate the force upon.
         * @param root The root of the octree.
         *
         * After this function is called, 'forceRecursive' is called.
         */

        public void calcForces(List<Particle> particles, double bounds)
        {
            Node root = constructTree(bounds, SimConstants.PARTICLES_PER_BOX, particles);

            List<int> integerList = Enumerable.Range(0, SimConstants.N_PARTICLES).ToList();
            Parallel.ForEach(integerList, i =>
            {
                forceRecursive(particles[i], root);
            });
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


        private Node constructTree(double bounds, int particlesPerBox, List<Particle> particles)
        {
            Node root = new Node(new List<double> { 0, 0, 0, }, bounds);
            foreach (Particle p in particles)
            {
                root.recursiveAddChild(p, particlesPerBox);
            }
            root.calcCOMRecursive();
            return root;
        }

        /*
         * @param p the force from all other particles in the node n (and all children nodes) on p will be calculated, and the acceleration
         *          of p will be set accordingly.
         * @param n the root node of the tree to be traversed.
         *
         * This is the recursive function that calculates the force of the entire tree, n, on the particle p. On average, only log(n) nodes
         * will be traversed when the function is called, which means that entire timestep can be calculated in O(n log n) time, since there
         * are n particles.
         */
        private void forceRecursive(Particle p, Node n)
        {
            double dx = n.com[0] - p.pos[0];
            double dy = n.com[1] - p.pos[1];
            double dz = n.com[2] - p.pos[2];

            double dist = Math.Pow(dx * dx + dy * dy + dz * dz, 0.5) + 1;


            bool cond = n.size / dist < SimConstants.BARNES_HUT_DELTA;

            //if cond is true, the node is sufficiently far away
            if (cond)
            {
                //calculate the force between p and n.
                double F = 6.67 * Math.Pow(10, -11) * n.mass * p.mass / dist / dist;
                double fx = F * dx / dist;
                double fy = F * dy / dist;
                double fz = F * dz / dist;

                //update the acceleration of p
                p.acc[0] += fx / p.mass;
                p.acc[1] += fy / p.mass;
                p.acc[2] += fz / p.mass;
            }
            else
            {
                foreach (Node child in n.children)
                {
                    forceRecursive(p, child);
                }
            }

            if (n.children.Count == 0 && cond == false)
            {
                foreach (Particle p1 in n.particles)
                {

                    double pdx = p1.pos[0] - p.pos[0];
                    double pdy = p1.pos[1] - p.pos[1];
                    double pdz = p1.pos[2] - p.pos[2];

                    double distCorrection = 1;

                    double pdist = Math.Pow(pdx * pdx + pdy * pdy + pdz * pdz, 0.5) + distCorrection;

                    if (pdist == distCorrection) { break; }

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
        public Node(List<double> origin, double size)
        {
            //Sets field variables of the node
            com = new List<double> { 0.0, 0.0, 0.0 };
            this.size = size;
            mass = 0;
            this.origin = new List<double> { origin[0], origin[1], origin[2] };
            particles = new List<Particle>();
            children = new List<Node>();

        }

        public void addChild(Node child) { children.Add(child); }

        public Node getChild(int index) { return children[index]; }

        /*
         * Recursive function for calculating the center of mass of this node. Fundamentally it works by finding the COM of all
         * children of the node, and then calculating a weighted average of these to compute the COM of this node.
         */
        public void calcCOMRecursive()
        {
            //Calculates the COM of all the children of the node.
            foreach (Node child in children)
            {
                child.calcCOMRecursive();
            }

            //Once children COMs have been calculated, the COM of this node is calculated.

            if (children.Count + particles.Count > 0)
            {
                calcCOM(); //calculates COM of this node. 
            }
            else //in the case where this node has no children and no particles within it, the COM is simply set to (0, 0, 0) to avoid
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
        public void calcCOM()
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
            foreach (Particle p in particles)
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
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1] + size / 2, origin[2] }, size / 2));
                children.Add(new Node(new List<double> { origin[0] + size / 2, origin[1] + size / 2, origin[2] + size / 2 }, size / 2));

                //Goes through each particle of the node and calls recursiveAddChild on the child that contains the particle.
                foreach (Particle p1 in particles)
                {
                    int index = 0;
                    if (p1.pos[0] > origin[0] + size / 2.0) { index += 4; }
                    if (p1.pos[1] > origin[1] + size / 2.0) { index += 2; }
                    if (p1.pos[2] > origin[2] + size / 2.0) { index += 1; }
                    children[index].recursiveAddChild(p1, particlesPerBox);
                }
                //Since all the particles are now distributed among the children of the node, the particle list of this node is cleared.
                particles.Clear();

            }
            else if (children.Count != 0) //If the condition does hold, but the 8 children have already been added, then recursiveAddChild
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
