// Input.pde - User input handling and controls

void keyPressed() {
  if (key == ' ') {
    showTrails = !showTrails;
  } else if (key == 'f' || key == 'F') {
    showForces = !showForces;
  } else if (key == 'g' || key == 'G') {
    showGas = !showGas;
  } else if (key == 'm' || key == 'M') {
    showMagnetic = !showMagnetic;
  } else if (key == 't' || key == 'T') {
    showTemperature = !showTemperature;
  } else if (key == 'c' || key == 'C') {
    showCoagulation = !showCoagulation;
  } else if (key == 'a' || key == 'A') {
    showAccretionStreams = !showAccretionStreams;
  } else if (key == 'w' || key == 'W') {
    showSpiralWaves = !showSpiralWaves;
  } else if (key == 'd' || key == 'D') {
    showCircumplanetaryDisks = !showCircumplanetaryDisks;
  } else if (key == 'r' || key == 'R') {
    resetSimulation();
  } else if (key == '1') {
    coagulationEfficiency = max(0.1, coagulationEfficiency - 0.1);
    println("Coagulation efficiency: " + coagulationEfficiency);
  } else if (key == '2') {
    coagulationEfficiency = 0.8; // Default
    println("Coagulation efficiency reset to default: " + coagulationEfficiency);
  } else if (key == '3') {
    coagulationEfficiency = min(1.5, coagulationEfficiency + 0.1);
    println("Coagulation efficiency: " + coagulationEfficiency);
  } else if (key == '4') {
    gasAccretionRate = max(0.0001, gasAccretionRate - 0.0005);
    println("Gas accretion rate: " + gasAccretionRate);
  } else if (key == '5') {
    gasAccretionRate = 0.001; // Default
    println("Gas accretion rate reset to default: " + gasAccretionRate);
  } else if (key == '6') {
    gasAccretionRate = min(0.01, gasAccretionRate + 0.0005);
    println("Gas accretion rate: " + gasAccretionRate);
  } else if (key == '7') {
    diskPlanetCouplingStrength = max(0.01, diskPlanetCouplingStrength - 0.02);
    println("Disk-planet coupling: " + diskPlanetCouplingStrength);
  } else if (key == '8') {
    diskPlanetCouplingStrength = 0.1; // Default
    println("Disk-planet coupling reset to default: " + diskPlanetCouplingStrength);
  } else if (key == '9') {
    diskPlanetCouplingStrength = min(0.5, diskPlanetCouplingStrength + 0.02);
    println("Disk-planet coupling: " + diskPlanetCouplingStrength);
  } else if (key == 'p' || key == 'P') {
    // Pause/unpause
    if (isLooping()) {
      noLoop();
      println("Simulation paused");
    } else {
      loop();
      println("Simulation resumed");
    }
  } else if (key == 's' || key == 'S') {
    // Save current state
    saveSimulationState();
  } else if (key == 'x' || key == 'X') {
    // Clear all trails
    for (Particle p : particles) {
      p.trail.clear();
    }
    for (Particle p : planetesimals) {
      p.trail.clear();
    }
    println("All trails cleared");
  } else if (key == 'z' || key == 'Z') {
    // Add random planetesimal
    addRandomPlanetesimal();
  }
}

void mousePressed() {
  // Enhanced particle creation with composition based on location
  float mass = random(5, 20);
  Particle newParticle = new Particle(mouseX, mouseY, mass);
  
  // Calculate realistic orbital velocity
  PVector toCenter = PVector.sub(centralStar.pos, newParticle.pos);
  float distance = toCenter.mag();
  
  if (distance > 60) {
    // Keplerian orbital velocity
    float orbitalSpeed = sqrt(G * centralStar.mass / distance);
    PVector vel = new PVector(-toCenter.y, toCenter.x);
    vel.normalize();
    vel.mult(orbitalSpeed * random(0.8, 1.2));
    
    // Add orbital eccentricity
    vel.add(PVector.random2D().mult(random(0.5, 2.0)));
    
    newParticle.vel = vel;
  }
  
  // Determine where to add the particle based on mass
  if (newParticle.mass < planetesimalFormationThreshold) {
    dust.add(newParticle);
    println("Added dust particle: Mass=" + nf(newParticle.mass, 1, 2));
  } else if (newParticle.mass < 15) {
    planetesimals.add(newParticle);
    println("Added planetesimal: Mass=" + nf(newParticle.mass, 1, 2) + 
            ", Composition=" + newParticle.compositionType);
  } else {
    particles.add(newParticle);
    println("Added planet: Mass=" + nf(newParticle.mass, 1, 2) + 
            ", Composition=" + newParticle.compositionType);
  }
}

void mouseDragged() {
  // Enhanced dust trail creation with realistic properties
  if (frameCount % 2 == 0) {
    float dustMass = random(0.3, 1.5);
    Particle dustParticle = new Particle(mouseX + random(-8, 8), mouseY + random(-8, 8), dustMass);
    
    // Give dust velocity based on mouse movement and local orbital motion
    PVector mouseVel = new PVector(mouseX - pmouseX, mouseY - pmouseY);
    dustParticle.vel = PVector.mult(mouseVel, 0.2);
    
    // Add some orbital component
    PVector toCenter = PVector.sub(centralStar.pos, dustParticle.pos);
    float distance = toCenter.mag();
    if (distance > 50) {
      float orbitalSpeed = sqrt(G * centralStar.mass / distance) * 0.3;
      PVector orbitalVel = new PVector(-toCenter.y, toCenter.x);
      orbitalVel.normalize();
      orbitalVel.mult(orbitalSpeed);
      dustParticle.vel.add(orbitalVel);
    }
    
    dust.add(dustParticle);
  }
}

void addRandomPlanetesimal() {
  float angle = random(TWO_PI);
  float distance = random(150, 400);
  float x = centralStar.pos.x + cos(angle) * distance;
  float y = centralStar.pos.y + sin(angle) * distance;
  
  float orbitalSpeed = sqrt(G * centralStar.mass / distance);
  PVector vel = new PVector(-sin(angle), cos(angle));
  vel.mult(orbitalSpeed * random(0.9, 1.1));
  
  Particle planetesimal = new Particle(x, y, random(planetesimalFormationThreshold, 12));
  planetesimal.vel = vel;
  planetesimals.add(planetesimal);
  
  println("Added planetesimal: Mass=" + nf(planetesimal.mass, 1, 2) + 
          ", Distance=" + nf(distance, 1, 0));
}

void resetSimulation() {
  particles.clear();
  dust.clear();
  planetesimals.clear();
  gasCells.clear();
  accretionStreams.clear();
  spiralWaves.clear();
  
  // Recreate everything with enhanced physics
  createGasDisk();
  createRealisticDustCloud(500);
  createSeedPlanetesimals(8);
  createInitialPlanets(3);
  
  println("=== SIMULATION RESET ===");
  println("Enhanced planet formation physics active:");
  println("- Dust coagulation and growth");
  println("- Planetesimal dynamics");
  println("- Gas accretion flows");
  println("- Disk-planet interactions");
  println("- Circumplanetary disk formation");
  println("- Spiral density waves");
}

void saveSimulationState() {
  // Enhanced save function including all new physics
  ArrayList<String> data = new ArrayList<String>();
  data.add("# Enhanced Planet Accretion Simulation State");
  data.add("# Format: Type,PosX,PosY,VelX,VelY,Mass,Temperature,Composition");
  
  // Save planets
  for (Particle p : particles) {
    String diskInfo = (p.circumDisk != null) ? "DISK" : "NONE";
    data.add("PLANET," + p.pos.x + "," + p.pos.y + "," + p.vel.x + "," + p.vel.y + "," + 
             p.mass + "," + p.temperature + "," + p.compositionType + "," + diskInfo);
  }
  
  // Save planetesimals
  for (Particle p : planetesimals) {
    data.add("PLANETESIMAL," + p.pos.x + "," + p.pos.y + "," + p.vel.x + "," + p.vel.y + "," + 
             p.mass + "," + p.temperature + "," + p.compositionType);
  }
  
  // Save dust (sample only, as there might be many)
  for (int i = 0; i < min(100, dust.size()); i++) {
    Particle d = dust.get(i);
    data.add("DUST," + d.pos.x + "," + d.pos.y + "," + d.vel.x + "," + d.vel.y + "," + 
             d.mass + "," + d.temperature + "," + d.compositionType);
  }
  
  // Save physics parameters
  data.add("PARAMS," + coagulationEfficiency + "," + gasAccretionRate + "," + 
           diskPlanetCouplingStrength + "," + stickingVelocity);
  
  String[] dataArray = new String[data.size()];
  for (int i = 0; i < data.size(); i++) {
    dataArray[i] = data.get(i);
  }
  saveStrings("data/enhanced_accretion_" + year() + month() + day() + hour() + minute() + ".csv", dataArray);
  println("Enhanced simulation state saved with " + data.size() + " entries!");
}
