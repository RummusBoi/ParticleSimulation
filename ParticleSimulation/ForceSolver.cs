using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ParticleSimulation
{
    public interface ForceSolver
    {
        void integrateAll(List<Particle> particles, double timestepSize, double bounds);
        void calcForces(List<Particle> particles, double bounds);
    }
}
