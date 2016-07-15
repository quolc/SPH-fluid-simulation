/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
  0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
PVector ijk2xyz(PVector v) {
  return new PVector(
    -GRID_SIZE/2 + GRID_SIZE * v.x / (GRID_N-1),
    -GRID_SIZE/2 + GRID_SIZE * v.y / (GRID_N-1),
    GRID_SIZE * v.z / (GRID_N-1)
  );
}

Triangle[] Polygonise(int x, int y, int z)
{
   int cubeindex;
   PVector[] vertlist = new PVector[12];
   PVector[] normlist = new PVector[12];

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   int[] x_o = {0,1,1,0,0,1,1,0};
   int[] y_o = {0,0,1,1,0,0,1,1};
   int[] z_o = {0,0,0,0,1,1,1,1};
   cubeindex = 0;
   for (int i=0; i<8; i++) {
     float f = field[x + x_o[i]][y + y_o[i]][z + z_o[i]];
     if (f < threshold) {
       cubeindex |= (1 << i);
     }
   }

   /* Cube is entirely in/out of the surface */
   if (edgeTable[cubeindex] == 0)
      return new Triangle[] {};

   /* Find the vertices where the surface intersects the cube */
   int[][] edge_vert_table = {
     {0, 1}, {1, 2}, {2, 3}, {3, 0},
     {4, 5}, {5, 6}, {6, 7}, {7, 4},
     {0, 4}, {1, 5}, {2, 6}, {3, 7}
   };
   for (int i=0; i<12; i++) {
     if ((edgeTable[cubeindex] & (1 << i)) > 0) {
       int v1 = edge_vert_table[i][0];
       int v2 = edge_vert_table[i][1];
       vertlist[i] = vertexInterp(
         x + x_o[v1], y + y_o[v1], z + z_o[v1],
         x + x_o[v2], y + y_o[v2], z + z_o[v2]
       );
       normlist[i] = normInterp(
         x + x_o[v1], y + y_o[v1], z + z_o[v1],
         x + x_o[v2], y + y_o[v2], z + z_o[v2]
       );
     }
   }

   /* Create the triangle */
   ArrayList<Triangle> tris = new ArrayList<Triangle>();
   for (int i=0; triTable[cubeindex][i] != -1; i += 3) {
     Triangle tri = new Triangle(
       vertlist[triTable[cubeindex][i]],
       vertlist[triTable[cubeindex][i+1]],
       vertlist[triTable[cubeindex][i+2]],
       normlist[triTable[cubeindex][i]],
       normlist[triTable[cubeindex][i+1]],
       normlist[triTable[cubeindex][i+2]]
     );
     tris.add(tri);
   }

   return tris.toArray(new Triangle[]{});
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
PVector vertexInterp(int i1, int j1, int k1, int i2, int j2, int k2)
{
  float mu;
  float valp1 = field[i1][j1][k1];
  float valp2 = field[i2][j2][k2];
 
  if (abs(threshold-valp1) < 0.00001)
    return(new PVector(i1, j1, k1));
  if (abs(threshold-valp2) < 0.00001)
    return(new PVector(i2, j2, k2));
  if (abs(valp1-valp2) < 0.00001)
    return(new PVector(i1, j1, k1));
  mu = (threshold - valp1) / (valp2 - valp1);
  
 return new PVector(
   i1 + mu * (i2 - i1),
   j1 + mu * (j2 - j1),
   k1 + mu * (k2 - k1)
 );
}
PVector normInterp(int i1, int j1, int k1, int i2, int j2, int k2)
{
  float mu;
  float valp1 = field[i1][j1][k1];
  float valp2 = field[i2][j2][k2];
  
  PVector norm1 = normals[i1][j1][k1];
  PVector norm2 = normals[i2][j2][k2];
   
  if (abs(threshold-valp1) < 0.00001)
    return norm1;
  if (abs(threshold-valp2) < 0.00001)
    return norm2;
  if (abs(valp1-valp2) < 0.00001)
    return norm1;
  mu = (threshold - valp1) / (valp2 - valp1);
  
  return new PVector(
    norm1.x + mu * (norm2.x - norm1.x),
    norm1.y + mu * (norm2.y - norm1.y),
    norm1.z + mu * (norm2.z - norm1.z)
  );
}