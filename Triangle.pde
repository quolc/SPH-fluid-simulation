class Triangle {
  public PVector[] ps;
  public PVector[] ns;
  
  public Triangle(PVector p1, PVector p2, PVector p3,
                  PVector n1, PVector n2, PVector n3) {
    this.ps = new PVector[3];
    this.ns = new PVector[3];
    this.ps[0] = p1;
    this.ps[1] = p2;
    this.ps[2] = p3;
    this.ns[0] = n1;
    this.ns[1] = n2;
    this.ns[2] = n3;
  }
}