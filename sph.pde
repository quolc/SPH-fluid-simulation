import java.util.*;

// simulation constants
int M = 10;
int N = M * M * M;
float dt = 0.005; // [s]

// physical constants
float volume = 1.2; // [L]
float part_mass = volume * 1.0 / N; // [kg]
float h = 0.05; // kernel width
float k = 0.80; // gass constant
float mu = 0.60; // viscosity constant
float g = 9.8; // gravitation
float sigma = 0.005; // tension
float rho_0 = 0;

// 0: rectangle pool, cube drop
// 1: rectangle pool, sphere drop
int exp_type = 0;

// boundary condition
float bin_size = 0.125;

// initial condition
float cube_size = pow(volume, 1.0/3) * 0.1;
float sphere_rad = pow(3/(4*PI)*volume/1000, 1.0/3);

boolean particle_on = true;   // enable particle rendering
boolean marching_on = false;  // enable surface rendering
boolean output_obj = false;   // enable writing .obj files

Particle[] ps = new Particle[N];
ArrayList<Triangle> tri_buf = new ArrayList<Triangle>();

// marching-cube related variables
boolean norm_interpolation = true; // rendering with normal interpolation (gouraud shading)
int GRID_N = 30;
float GRID_SIZE = bin_size*2; // [-GRID_SIZE/2, GRID_SIZE/2] * [-GRID_SIZE/2, GRID_SIZE/2] * [0, GRID_SIZE]
float threshold = 0.6;

float[][][] field = new float[GRID_N][GRID_N][GRID_N];
PVector[][][] normals = new PVector[GRID_N][GRID_N][GRID_N];

void setup() {
  size(800, 800, P3D);
  frameRate(30);

  sphereDetail(6);

  // initialize particle initial states
  switch (exp_type) {
    case 0:
      float random_range = 0.00;
      for (int i=0; i<M; i++) {
        for (int j=0; j<M; j++) {
          for (int k=0; k<M; k++) {
            int ind = i*M*M + j*M + k;
            ps[ind] = new Particle(part_mass, 
              (i - 0.5*(M-1) + random(-random_range, random_range)) * (cube_size / M),
              (j - 0.5*(M-1) + random(-random_range, random_range)) * (cube_size / M), 
              (k + random(-random_range, random_range)) * (cube_size / M) + 0.12);
          }
        }
      }
      break;
    case 1:
      for (int i=0; i<N; i++) {
        while (true) {
          float x = random(-sphere_rad, sphere_rad);
          float y = random(-sphere_rad, sphere_rad);
          float z = random(-sphere_rad, sphere_rad);
          if (sq(x)+sq(y)+sq(z) < sq(sphere_rad)) {
            ps[i] = new Particle(part_mass, x, y, z + sphere_rad*2);
            break;
          }
        }
      }
      break;
  }
  
  for (int i=0; i<GRID_N; i++) {
    for (int j=0; j<GRID_N; j++) {
      for (int k=0; k<GRID_N; k++) {
        field[i][j][k] = 0;
      }
    }
  }
}

// example of dynamic field
void computeField() {
  // use color field as scalar field
  for (int i=0; i<GRID_N; i++) {
    for (int j=0; j<GRID_N; j++) {
      for (int k=0; k<GRID_N; k++) {
        field[i][j][k] = 0;
      }
    }
  }
  
  for (Particle p : ps) {
    float i_f = (p.x+GRID_SIZE/2) / GRID_SIZE * GRID_N;
    float j_f = (p.y+GRID_SIZE/2) / GRID_SIZE * GRID_N;
    float k_f = p.z / GRID_SIZE * GRID_N;
    int ne = 2;
    for (int i=max(int(i_f)-ne, 0); i<=int(i_f)+ne && i<GRID_N; i++) {
      for (int j=max(int(j_f)-ne, 0); j<=int(j_f)+ne && j<GRID_N; j++) {
        for (int k=max(int(k_f)-ne, 0); k<=int(k_f)+ne && k<GRID_N; k++) {
          float x = -GRID_SIZE/2 + GRID_SIZE / (GRID_N-1) * i;
          float y = -GRID_SIZE/2 + GRID_SIZE / (GRID_N-1) * j;
          float z = GRID_SIZE / (GRID_N-1) * k;
          field[i][j][k] += p.mass / p.rho * wSpiky(
            new Particle(0, x, y, z), p
          );
        }
      }
    }
  }
}

// calculates gradient of scalar field for normal interpolation
void computeNormals() {
  for (int i=0; i<GRID_N; i++) {
    for (int j=0; j<GRID_N; j++) {
      for (int k=0; k<GRID_N; k++) {
        float dx = GRID_SIZE / (GRID_N-1);
        float dudx = 0, dudy = 0, dudz = 0;
        if (i == 0) dudx = (field[i+1][j][k] - field[i][j][k]) / dx;
        else if (i == GRID_N-1) dudx = (field[i][j][k] - field[i-1][j][k]) / dx;
        else dudx = (field[i+1][j][k] - field[i-1][j][k]) / (2*dx);
        if (j == 0) dudy = (field[i][j+1][k] - field[i][j][k]) / dx;
        else if (j == GRID_N-1) dudy = (field[i][j][k] - field[i][j-1][k]) / dx;
        else dudy = (field[i][j+1][k] - field[i][j-1][k]) / (2*dx);
        if (k == 0) dudz = (field[i][j][k+1] - field[i][j][k]) / dx;
        else if (k == GRID_N-1) dudz = (field[i][j][k] - field[i][j][k-1]) / dx;
        else dudz = (field[i][j][k+1] - field[i][j][k-1]) / (2*dx);
        normals[i][j][k] = new PVector(-dudx, -dudy, -dudz);
        normals[i][j][k].normalize();
      }
    }
  }
}

void draw() {
  long ta = millis();
  background(255);

  // do simulation
  for (int t=0; t<1; t++) {
    simulate();
  }
  
  // compute visualization
  computeField(); 
  computeNormals();

  // rendering
  background(255);  

  // camera setup
  float fov = PI/3.0;
  perspective(fov, float(width)/float(height), 
    0.01, 1000);

  float cam_r = 0.5;
  float cam_theta = PI/6;
  float cam_phi = TWO_PI / 800 * frameCount ;
  camera(cam_r * sin(cam_theta) * cos(cam_phi), 
    cam_r * sin(cam_theta) * sin(cam_phi), 
    cam_r * cos(cam_theta), 
    0, 0, M/2 * 0.009, 
    0, 0, -1);

  directionalLight(161, 161, 161, -0.2, 0.3, -0.5);
  ambientLight(91, 110, 127);

  // x-y-z axes
  /*
  stroke(255, 0, 0);
  line(0, 0, 0, 100, 0, 0);
  stroke(0, 255, 0);
  line(0, 0, 0, 0, 100, 0);
  stroke(0, 0, 255);
  line(0, 0, 0, 0, 0, 100);
//  */

  // ground
  noStroke();
  switch (exp_type) {
    case 0:
    case 1:
      ambient(127,127,127);
      rect(-bin_size, -bin_size, bin_size*2, bin_size*2);
      pushMatrix();
      rotateX(PI/2);
      translate(0, 0, bin_size);
      rect(-bin_size, 0, bin_size*2, bin_size/2);
      translate(0, 0, -2*bin_size);
      rect(-bin_size, 0, bin_size*2, bin_size/2);
      popMatrix();
      pushMatrix();
      rotateY(PI/2);
      translate(0, 0, bin_size);
      rect(-bin_size/2, -bin_size, bin_size/2, bin_size*2);
      translate(0, 0, -2*bin_size);
      rect(-bin_size/2, -bin_size, bin_size/2, bin_size*2);
      popMatrix();
      break;
    case 2:
      ambient(127,127,127);
      pushMatrix();
      rotateX(PI/2);
      translate(0, 0, cube_size/2);
      rect(-cube_size/2, 0, cube_size, cube_size);
      translate(0, 0, -cube_size);
      rect(-cube_size/2, 0, cube_size, cube_size);
      popMatrix();
      pushMatrix();
      rotateY(PI/2);
      translate(0, 0, cube_size/2);
      rect(-cube_size, -cube_size/2, cube_size, cube_size);
      translate(0, 0, -cube_size);
      rect(-cube_size, -cube_size/2, cube_size, cube_size);
      popMatrix();
      break;
  }

  // particles
  if (particle_on) {
    ambient(127,127,127);
    for (Particle p : ps) {
      pushMatrix();
      translate(p.x, p.y, p.z);
      sphere(pow(p.mass / p.rho, 1.0/3));
      popMatrix();
    }
  }

  // marching-cube
  if (marching_on) {
    tri_buf.clear();
    stroke(0,0,0);
    strokeWeight(0.5);
    noStroke();
    ambient(127,200,255);
    specular(255, 255, 255);
    beginShape(TRIANGLES);
    for (int i=0; i<GRID_N-1; i++) {
      for (int j=0; j<GRID_N-1; j++) {
        for (int k=0; k<GRID_N-1; k++) {
          Triangle[] tris = Polygonise(i, j, k);
          tri_buf.addAll(Arrays.asList(tris));
          for (Triangle tri : tris) {
            for (int l=0; l<3; l++) {
              PVector p = ijk2xyz(tri.ps[l]);
              if (norm_interpolation) {
                normal(tri.ns[l].x, tri.ns[l].y, tri.ns[l].z);
              }
              vertex(p.x, p.y, p.z);
            }
          }
        }
      }
    }
    endShape();
    if (output_obj) {
      outputObj();
    }
  }

  // collapsed time
  long tb = millis();
  if (frameCount % 30 == 0) {
    println("computation time: ", tb - ta);
  }
  
  save(frameCount + ".png");
}

// smoothing kernel functions
public float w(Particle pa, Particle pb) {
  float r = sqrt(sq(pa.x - pb.x) + sq(pa.y - pb.y) + sq(pa.z - pb.z));

  if (r > h) return 0;
  return (float)315/(64*PI*pow(h, 9)) * pow((h*h - r*r), 3);
}
public float wSpiky(Particle pa, Particle pb) {
  float r = sqrt(sq(pa.x - pb.x) + sq(pa.y - pb.y) + sq(pa.z - pb.z));

  if (r > h) return 0;
  return (float)15/(PI*pow(h, 6)) * pow((h - r), 3);
}

public PVector gradw(Particle pa, Particle pb) {
  float r = sqrt(sq(pa.x - pb.x) + sq(pa.y - pb.y) + sq(pa.z - pb.z));

  if (r > h) return new PVector(0, 0, 0);
  float coef = -6 * (float)315/(64*PI*pow(h, 9)) * pow((h*h - r*r), 2);
  return new PVector(
    coef * (pa.x - pb.x), 
    coef * (pa.y - pb.y), 
    coef * (pa.z - pb.z));
}
public PVector gradwSpiky(Particle pa, Particle pb) {
  float r = sqrt(sq(pa.x - pb.x) + sq(pa.y - pb.y) + sq(pa.z - pb.z));

  if (r > h) return new PVector(0, 0, 0);
  float coef = -3 * (float)15/(PI*pow(h, 6))*sq(h-r)/r;
  return new PVector(
    coef * (pa.x - pb.x), 
    coef * (pa.y - pb.y), 
    coef * (pa.z - pb.z));
}

public float laplw(Particle pa, Particle pb) {
  float r2 = sq(pa.x - pb.x) + sq(pa.y - pb.y) + sq(pa.z - pb.z);
  if (r2 > h*h) return 0;
  
  return (float)315/(64*PI*pow(h, 9)) * (h*h - r2) * (42 * r2 - 18*h*h);
}
public float laplwVisc(Particle pa, Particle pb) {
  float r = sqrt(sq(pa.x - pb.x) + sq(pa.y - pb.y) + sq(pa.z - pb.z));

  if (r > h) return 0;
  return 45/(PI*pow(h, 6))*(h-r);
}

void simulate() {
  // a. calculate mass densities of each particle
  rho_0 = 0;
  for (Particle p : ps) {
    p.rho = 0;
    for (Particle q : ps) {
      p.rho += q.mass * w(p, q);
    }
    rho_0 += p.rho / N;
  }

  // b. calculate pressure at each particle
  for (Particle p : ps) {
    p.pressure = k * (p.rho - rho_0);
  }

  // c. calculate acceleration of particles
  for (Particle p : ps) {
    float fx=0, fy=0, fz=0;

    // 1. pressure term
    float fpx=0, fpy=0, fpz=0;
    for (Particle q : ps) {
      if (p == q) continue;
      PVector gw = gradwSpiky(p, q);
      fpx += -q.mass * (p.pressure + q.pressure) / (2 * p.rho) * gw.x;
      fpy += -q.mass * (p.pressure + q.pressure) / (2 * p.rho) * gw.y;
      fpz += -q.mass * (p.pressure + q.pressure) / (2 * p.rho) * gw.z;
    }
    fx += fpx;
    fy += fpy;
    fz += fpz;

    // 2. viscosity term
    float fvx=0, fvy=0, fvz=0;
    for (Particle q : ps) {
      float lw = laplwVisc(p, q);
      fvx += mu * q.mass * (q.vx - p.vx) / q.rho * lw;
      fvy += mu * q.mass * (q.vy - p.vy) / q.rho * lw;
      fvz += mu * q.mass * (q.vz - p.vz) / q.rho * lw;
    }
    fx += fvx;
    fy += fvy;
    fz += fvz;

    // 3. surface tension
    float ftx=0, fty=0, ftz=0;
    float nx=0, ny=0, nz=0; // color field normal
    for (Particle q : ps) {
      PVector gw = gradw(p, q);
      nx += q.mass / q.rho * gw.x;
      ny += q.mass / q.rho * gw.y;
      nz += q.mass / q.rho * gw.z;
    }
    float nabs = sqrt(nx*nx + ny*ny + nz*nz);
    float nthreshold = 0.1;
    if (nabs > nthreshold) {
      float lapl = 0;
      for (Particle q : ps) {
        float lw = laplw(p, q);
        lapl += q.mass / q.rho * lw;
      }
      ftx = -sigma * lapl * nx / nabs;
      fty = -sigma * lapl * ny / nabs;
      ftz = -sigma * lapl * nz / nabs;
      fx += ftx;
      fy += fty;
      fz += ftz;
    }
//     if (frameCount % 10 == 0) println(nabs);

    // 4. external force
    fz -= g * p.rho;

    p.ax = fx / p.rho;
    p.ay = fy / p.rho;
    p.az = fz / p.rho;
  }

  // d. time integration (ToDo: Runge-Kutta)
  for (Particle p : ps) {
    p.vx_next = p.vx + p.ax * dt;
    p.vy_next = p.vy + p.ay * dt;
    p.vz_next = p.vz + p.az * dt;
    p.x_next = p.x + p.vx_next * dt;
    p.y_next = p.y + p.vy_next * dt;
    p.z_next = p.z + p.vz_next * dt;

    // repulsion
    switch (exp_type) {
      case 0:
      case 1:
        if (p.x_next < -bin_size) {
          p.vx_next = -p.vx_next;
          p.x_next = -bin_size + (-bin_size - p.x_next);
        }
        if (p.x_next > bin_size) {
          p.vx_next = -p.vx_next;
          p.x_next = bin_size - (p.x_next - bin_size);
        }
        if (p.y_next < -bin_size) {
          p.vy_next = -p.vy_next;
          p.y_next = -bin_size + (-bin_size - p.y_next);
        }
        if (p.y_next > bin_size) {
          p.vy_next = -p.vy_next;
          p.y_next = bin_size - (p.y_next - bin_size);
        }
        break;
      case 2:
        if (p.x_next < -cube_size/2) {
          p.vx_next = -p.vx_next;
          p.x_next = -cube_size/2 + (-cube_size/2 - p.x_next);
        }
        if (p.x_next > cube_size/2*3) {
          p.vx_next = -p.vx_next;
          p.x_next = cube_size/2 - (p.x_next - cube_size/2);
        }
        if (p.y_next < -cube_size/2) {
          p.vy_next = -p.vy_next;
          p.y_next = -cube_size/2 + (-cube_size/2 - p.y_next);
        }
        if (p.y_next > cube_size/2) {
          p.vy_next = -p.vy_next;
          p.y_next = cube_size/2 - (p.y_next - cube_size/2);
        }
        break;
    }
    if (p.z_next < 0) {
      p.vz_next = -p.vz_next;
      p.z_next = -p.z_next;
    }

    // apply
    p.vx = p.vx_next;
    p.vy = p.vy_next;
    p.vz = p.vz_next;
    p.x = p.x_next;
    p.y = p.y_next;
    p.z = p.z_next;
  }
}

void outputObj() {
  ArrayList<String> strs = new ArrayList<String>();
  
  // v entry
  for (int i=0; i<tri_buf.size(); i++) {
    Triangle tri = tri_buf.get(i);
    for (int j=0; j<3; j++) {
      strs.add(String.format("v %f %f %f", ijk2xyz(tri.ps[j]).y, ijk2xyz(tri.ps[j]).z, ijk2xyz(tri.ps[j]).x));
    }
  }
  
  // n entry
  for (int i=0; i<tri_buf.size(); i++) {
    Triangle tri = tri_buf.get(i);
    for (int j=0; j<3; j++) {
      strs.add(String.format("vn %f %f %f", tri.ns[j].y, tri.ns[j].z, tri.ns[j].x));
    }
  }
  
  // f entry
  for (int i=0; i<tri_buf.size(); i++) {
    Triangle tri = tri_buf.get(i);
    strs.add(String.format("f %d//%d %d//%d %d//%d", i*3+1, i*3+1, i*3+2, i*3+2, i*3+3, i*3+3));
  }
  
  saveStrings(frameCount+".obj", strs.toArray(new String[]{}));
}