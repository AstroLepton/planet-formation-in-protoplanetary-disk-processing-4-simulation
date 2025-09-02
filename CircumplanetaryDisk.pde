// CircumplanetaryDisk class - Disk around massive planets

class CircumplanetaryDisk {
  Particle planet;
  ArrayList<Particle> diskParticles;
  float innerRadius, outerRadius;
  float angularMomentum;
  
  CircumplanetaryDisk(Particle p) {
    planet = p;
    diskParticles = new ArrayList<Particle>();
    innerRadius = planet.radius * 3;
    outerRadius = pow(planet.mass / (3 * centralStar.mass), 1.0/3.0) * 
                  PVector.dist(planet.pos, centralStar.pos) * 0.4; // Hill sphere fraction
    angularMomentum = 0;
  }
  
  void update() {
    // Remove particles that have been accreted
    for (int i = diskParticles.size()-1; i >= 0; i--) {
      Particle dp = diskParticles.get(i);
      float distToPlanet = PVector.dist(dp.pos, planet.pos);
      
      if (distToPlanet < planet.radius * 1.5) {
        // Accrete onto planet
        planet.mass += dp.mass;
        planet.radius = pow(planet.mass * 0.75 / PI, 1.0/3.0) * 2;
        diskParticles.remove(i);
      } else if (distToPlanet > outerRadius) {
        // Particle escapes disk
        diskParticles.remove(i);
      }
    }
    
    // Update disk particle orbits
    for (Particle dp : diskParticles) {
      updateDiskParticleOrbit(dp);
    }
  }
  
  void updateDiskParticleOrbit(Particle dp) {
    PVector toPlanet = PVector.sub(planet.pos, dp.pos);
    float distance = toPlanet.mag();
    
    // Gravitational force from planet
    float planetForce = G * planet.mass * dp.mass / (distance * distance);
    toPlanet.normalize();
    dp.applyForce(PVector.mult(toPlanet, planetForce));
    
    // Orbital decay due to gas drag (simplified)
    dp.vel.mult(0.999);
  }
  
  void display() {
    // Draw disk structure
    stroke(red(planet.col), green(planet.col), blue(planet.col), 60);
    strokeWeight(2);
    noFill();
    ellipse(planet.pos.x, planet.pos.y, outerRadius*2, outerRadius*2);
    
    stroke(red(planet.col), green(planet.col), blue(planet.col), 100);
    strokeWeight(1);
    ellipse(planet.pos.x, planet.pos.y, innerRadius*2, innerRadius*2);
    
    // Draw disk particles
    for (Particle dp : diskParticles) {
      fill(red(planet.col), green(planet.col), blue(planet.col), 150);
      noStroke();
      ellipse(dp.pos.x, dp.pos.y, 3, 3);
    }
  }
  
  void addParticle(Particle p) {
    diskParticles.add(p);
  }
}
