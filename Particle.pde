// Particle class - Main physics body in the simulation

class Particle {
  PVector pos, vel, acc;
  float mass, radius, temperature, density;
  float stickinessCoeff, porosityFactor;
  color col;
  ArrayList<PVector> trail;
  boolean isInHillSphere;
  float eccentricity, semiMajorAxis;
  String compositionType; // "rocky", "icy", "mixed"
  CircumplanetaryDisk circumDisk;
  float lastCollisionTime;
  float lastDistanceFromStar;
  
  Particle(float x, float y, float m) {
    pos = new PVector(x, y);
    vel = new PVector(random(-1, 1), random(-1, 1));
    acc = new PVector(0, 0);
    mass = m;
    radius = pow(mass * 0.75 / PI, 1.0/3.0) * 2;
    
    // Safe temperature calculation - check if centralStar exists
    if (centralStar != null) {
      temperature = calculateTemperature(pos);
    } else {
      temperature = 300; // Default temperature if star not yet initialized
    }
    
    density = mass / (4.0/3.0 * PI * pow(radius, 3));
    col = getParticleColor();
    trail = new ArrayList<PVector>();
    isInHillSphere = false;
    
    // Safe orbital elements calculation
    if (centralStar != null) {
      calculateOrbitalElements();
      determineComposition();
    } else {
      // Set default values if star not yet initialized
      eccentricity = 0.1;
      semiMajorAxis = 100;
      compositionType = "mixed";
    }
    
    porosityFactor = 1.0;
    stickinessCoeff = calculateStickiness();
    lastCollisionTime = 0;
    lastDistanceFromStar = 0;
    
    // Create circumplanetary disk for massive bodies
    if (mass > 20) {
      circumDisk = new CircumplanetaryDisk(this);
    }
  }
  
  void determineComposition() {
    if (centralStar == null) {
      compositionType = "mixed";
      return;
    }
    
    float distFromStar = PVector.dist(pos, centralStar.pos);
    if (distFromStar < 180) {
      compositionType = "rocky";
      density *= 1.5; // Higher density for rocky materials
    } else if (distFromStar > 350) {
      compositionType = "icy";
      density *= 0.7; // Lower density for icy materials
      stickinessCoeff *= 1.3; // Ice sticks better at low temperatures
    } else {
      compositionType = "mixed";
    }
  }
  
  float calculateStickiness() {
    float baseStickiness = 0.5;
    if (temperature < 273) { // Below water freezing point
      baseStickiness *= 1.5; // Ice is stickier
    }
    if (temperature < 150) { // Very cold, other volatiles condense
      baseStickiness *= 2.0;
    }
    
    // Safe access to global coagulationEfficiency
    try {
      return baseStickiness * coagulationEfficiency;
    } catch (Exception e) {
      return baseStickiness * 0.8; // Default value if global not accessible
    }
  }
  
  void calculateOrbitalElements() {
    if (centralStar == null) {
      return;
    }
    
    PVector r = PVector.sub(pos, centralStar.pos);
    PVector v = vel.copy();
    
    float mu = G * centralStar.mass;
    semiMajorAxis = 1.0 / (2.0/r.mag() - v.magSq()/mu);
    
    // Simplified eccentricity calculation using cross product
    PVector h = cross(r, v);
    float h_mag = h.mag();
    eccentricity = sqrt(1 + 2 * (v.magSq()/2 - mu/r.mag()) * h_mag*h_mag / (mu*mu));
    eccentricity = constrain(eccentricity, 0, 0.9);
  }
  
  color getParticleColor() {
    // Color based on composition and temperature
    float tempFactor = map(temperature, 100, 1000, 0, 1);
    float massFactor = map(mass, 0, 50, 0, 1);
    
    if (compositionType != null && compositionType.equals("rocky")) {
      return color(200 + mass*2, 100 + mass*2, 50); // Reddish-brown for rocky
    } else if (compositionType != null && compositionType.equals("icy")) {
      return color(150, 200 + mass, 255 - mass); // Bluish for icy
    } else {
      return color(180 + mass, 180 + mass, 150 + mass*2); // Grayish for mixed
    }
  }
  
  void update() {
    vel.add(acc);
    vel.limit(8); // Slightly higher velocity limit
    pos.add(vel);
    acc.mult(0);
    
    // Update temperature based on position - safe access
    if (centralStar != null) {
      temperature = calculateTemperature(pos);
    }
    col = getParticleColor();
    
    // Update circumplanetary disk
    if (circumDisk != null) {
      circumDisk.update();
    }
    
    // Update stickiness based on current temperature
    stickinessCoeff = calculateStickiness();
    
    // Store trail
    try {
      if (showTrails && frameCount % 3 == 0) {
        trail.add(pos.copy());
        if (trail.size() > 60) trail.remove(0);
      }
    } catch (Exception e) {
      // Default trail behavior if global variable not accessible
      trail.add(pos.copy());
      if (trail.size() > 60) trail.remove(0);
    }
    
    // Realistic boundary conditions (no wrap-around) - safe access
    if (centralStar != null) {
      float distFromCenter = PVector.dist(pos, centralStar.pos);
      try {
        if (distFromCenter > diskOuterRadius * 1.2) {
          // Particle escapes the system
          vel.mult(0.95);
        }
      } catch (Exception e) {
        // Use default boundary if global variable not accessible
        if (distFromCenter > 500) {
          vel.mult(0.95);
        }
      }
    }
  }
  
  void display() {
    // Draw circumplanetary disk first
    try {
      if (circumDisk != null && showCircumplanetaryDisks) {
        circumDisk.display();
      }
    } catch (Exception e) {
      // Default circumplanetary disk drawing if global variable not accessible
      if (circumDisk != null) {
        circumDisk.display();
      }
    }
    
    // Draw Hill sphere if particle is massive enough
    try {
      if (mass > 15 && showForces && centralStar != null) {
        float hillRadius = pow(mass / (3 * centralStar.mass), 1.0/3.0) * 
                          PVector.dist(pos, centralStar.pos) * 50; // Scaled for visibility
        stroke(255, 50);
        strokeWeight(1);
        noFill();
        ellipse(pos.x, pos.y, hillRadius*2, hillRadius*2);
      }
    } catch (Exception e) {
      // Default Hill sphere drawing if global variable not accessible
      if (mass > 15 && centralStar != null) {
        float hillRadius = pow(mass / (3 * centralStar.mass), 1.0/3.0) * 
                          PVector.dist(pos, centralStar.pos) * 50; // Scaled for visibility
        stroke(255, 50);
        strokeWeight(1);
        noFill();
        ellipse(pos.x, pos.y, hillRadius*2, hillRadius*2);
      }
    }
    
    // Draw trail with orbital path
    try {
      if (showTrails) {
        stroke(red(col), green(col), blue(col), 120);
        strokeWeight(1 + mass/20);
        noFill();
        beginShape();
        for (PVector p : trail) {
          vertex(p.x, p.y);
        }
        endShape();
      }
    } catch (Exception e) {
      // Default trail drawing if global variable not accessible
      if (trail.size() > 0) {
        stroke(red(col), green(col), blue(col), 120);
        strokeWeight(1 + mass/20);
        noFill();
        beginShape();
        for (PVector p : trail) {
          vertex(p.x, p.y);
        }
        endShape();
      }
    }
    
    // Draw particle with realistic size
    fill(col);
    stroke(255, 100);
    strokeWeight(0.5);
    ellipse(pos.x, pos.y, radius*2, radius*2);
    
    // Draw velocity vector
    try {
      if (showForces) {
        stroke(0, 255, 0, 150);
        strokeWeight(1);
        PVector velVis = PVector.mult(vel, 10);
        line(pos.x, pos.y, pos.x + velVis.x, pos.y + velVis.y);
        
        // Draw acceleration vector
        stroke(255, 0, 0, 150);
        PVector accVis = PVector.mult(acc, 1000);
        line(pos.x, pos.y, pos.x + accVis.x, pos.y + accVis.y);
      }
    } catch (Exception e) {
      // Default force vector drawing if global variable not accessible
      stroke(0, 255, 0, 150);
      strokeWeight(1);
      PVector velVis = PVector.mult(vel, 10);
      line(pos.x, pos.y, pos.x + velVis.x, pos.y + velVis.y);
      
      // Draw acceleration vector
      stroke(255, 0, 0, 150);
      PVector accVis = PVector.mult(acc, 1000);
      line(pos.x, pos.y, pos.x + accVis.x, pos.y + accVis.y);
    }
  }
  
  void applyForce(PVector force) {
    PVector f = PVector.div(force, mass);
    acc.add(f);
  }
}
