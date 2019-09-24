using System;

public class SimConstants
{

    public static double DT = 0.001;
    public static int N_PARTICLES = 100;
    public static double BARNES_HUT_DELTA = 0.5;
    public static int SIMSIZE = 5000;
    public static int PARTICLES_PER_BOX = 1;
    public static int STEPS_PER_FRAME = 10;
    public static double ROTATION_SPEED = 0.005;
    public static int IMG_WIDTH = 1000;
    public static int IMG_HEIGHT = 1000;
    public static double PARTICLE_MASS = Math.Pow(10, 18);


    public SimConstants() { return;  }
}