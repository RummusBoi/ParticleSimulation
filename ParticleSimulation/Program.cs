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
        ForceSolver forceSolver;

        Random rnd; //rnd is used to generate starting position and velocity of all the particles.
        Node root; //root node of the octree.
        public static Form1 form1;

        //Instead of a constructor, the 'init' function is used.
        public Sim() { return; }

        
       
        /*
         * @param form The windows form that the particles will be drawn upon.
         */
        public void init(Form1 form, ForceSolver forceSolver)
        {
            form1 = form;
            particles = new List<Particle>();
            this.forceSolver = forceSolver;
            
            rnd = new Random();

            //generates all particles with random initial parameters, using the rnd variable:
            for (int i = 0; i < SimConstants.N_PARTICLES; i++)
            {
                particles.Add(new Particle(300, 300, 300, 600.0, rnd));
            }
        }

        /**
         * @param timesteps The amount of timesteps that the simulation is to be run before returning.
         *
         * Giving the starting position and velocity of all particles, the simulation runs 'timesteps' amount of steps into the future
         * before returning. On the final iteration, that is, when t = timesteps - 1, all the particles are plotted on the Windows Form.
         */
        public void run(int timesteps)
        {

            for(int t = 0; t < timesteps; t++)
            {
                forceSolver.calcForces(particles, SimConstants.SIMSIZE);

                forceSolver.integrateAll(particles, SimConstants.DT, SimConstants.SIMSIZE);
                

                if (t == timesteps - 1)
                {
                    for (int i = 0; i < SimConstants.N_PARTICLES; i++)
                    {
                        form1.drawBox(0, 0, 0, SimConstants.SIMSIZE);
                        form1.drawParticle((int)particles[i].pos[0], (int)particles[i].pos[1], (int)particles[i].pos[2], Color.Black);
                    }
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
    }
}