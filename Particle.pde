class Particle {
  // static parameters
  public float mass; // [kg]

  // dynamic parameters
  public float x, y, z;
  public float x_next, y_next, z_next;
  public float vx, vy, vz;
  public float vx_next, vy_next, vz_next;
  public float ax, ay, az;
  public float rho; // [kg/m^3]
  public float pressure;

  public Particle(float mass, float x, float y, float z) {
    this.mass = mass;
    this.x = x;
    this.y = y;
    this.z = z;
  }
}