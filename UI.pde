// UI.pde - User interface and information display

void drawUI() {
  fill(255, 220);
  textAlign(LEFT);
  textSize(12);
  
  // Main statistics
  text("=== ENHANCED ACCRETION SIMULATION ===", 10, 20);
  text("Planets: " + particles.size(), 10, 40);
  text("Planetesimals: " + planetesimals.size(), 10, 55);
  text("Dust: " + dust.size(), 10, 70);
  text("Gas Cells: " + gasCells.size(), 10, 85);
  text("Accretion Streams: " + accretionStreams.size(), 10, 100);
  text("Spiral Waves: " + spiralWaves.size(), 10, 115);
  
  // Physics toggles
  text("=== VISUALIZATIONS ===", 10, 140);
  text("Trails: " + (showTrails ? "ON" : "OFF") + " (SPACE)", 10, 155);
  text("Forces: " + (showForces ? "ON" : "OFF") + " (F)", 10, 170);
  text("Gas: " + (showGas ? "ON" : "OFF") + " (G)", 10, 185);
  text("Magnetic: " + (showMagnetic ? "ON" : "OFF") + " (M)", 10, 200);
  text("Temperature: " + (showTemperature ? "ON" : "OFF") + " (T)", 10, 215);
  text("Coagulation: " + (showCoagulation ? "ON" : "OFF") + " (C)", 10, 230);
  text("Accretion Streams: " + (showAccretionStreams ? "ON" : "OFF") + " (A)", 10, 245);
  text("Spiral Waves: " + (showSpiralWaves ? "ON" : "OFF") + " (W)", 10, 260);
  text("Circumplanetary Disks: " + (showCircumplanetaryDisks ? "ON" : "OFF") + " (D)", 10, 275);
  
  // Physics parameters
  text("=== PHYSICS PARAMETERS ===", 10, 300);
  text("Coagulation Efficiency: " + nf(coagulationEfficiency, 1, 2) + " (1-3)", 10, 315);
  text("Gas Accretion Rate: " + nf(gasAccretionRate, 1, 4) + " (4-6)", 10, 330);
  text("Disk-Planet Coupling: " + nf(diskPlanetCouplingStrength, 1, 2) + " (7-9)", 10, 345);
  text("Sticking Velocity: " + nf(stickingVelocity, 1, 1) + " m/s", 10, 360);
  
  // Largest bodies info
  text("=== LARGEST PLANETS ===", 10, 385);
  // Create a copy and sort it
  ArrayList<Particle> sortedParticles = new ArrayList<Particle>();
  for (Particle p : particles) {
    sortedParticles.add(p);
  }
  // Simple bubble sort for Processing 4 compatibility
  for (int i = 0; i < sortedParticles.size(); i++) {
    for (int j = i + 1; j < sortedParticles.size(); j++) {
      if (sortedParticles.get(i).mass < sortedParticles.get(j).mass) {
        Particle temp = sortedParticles.get(i);
        sortedParticles.set(i, sortedParticles.get(j));
        sortedParticles.set(j, temp);
      }
    }
  }
  
  for (int i = 0; i < min(3, sortedParticles.size()); i++) {
    Particle p = sortedParticles.get(i);
    float distFromStar = PVector.dist(p.pos, centralStar.pos);
    String diskInfo = (p.circumDisk != null) ? " [DISK]" : "";
    text("Planet " + (i+1) + ": M=" + nf(p.mass, 1, 1) + 
         ", R=" + nf(distFromStar, 1, 0) + 
         ", T=" + nf(p.temperature, 1, 0) + "K" +
         ", " + p.compositionType + diskInfo, 10, 405 + i*15);
  }
  
  // Planetesimal info
  if (planetesimals.size() > 0) {
    text("=== LARGEST PLANETESIMALS ===", 10, 460);
    ArrayList<Particle> sortedPlanetesimals = new ArrayList<Particle>();
    for (Particle p : planetesimals) {
      sortedPlanetesimals.add(p);
    }
    // Simple bubble sort
    for (int i = 0; i < sortedPlanetesimals.size(); i++) {
      for (int j = i + 1; j < sortedPlanetesimals.size(); j++) {
        if (sortedPlanetesimals.get(i).mass < sortedPlanetesimals.get(j).mass) {
          Particle temp = sortedPlanetesimals.get(i);
          sortedPlanetesimals.set(i, sortedPlanetesimals.get(j));
          sortedPlanetesimals.set(j, temp);
        }
      }
    }
    
    for (int i = 0; i < min(2, sortedPlanetesimals.size()); i++) {
      Particle p = sortedPlanetesimals.get(i);
      text("Planetesimal " + (i+1) + ": M=" + nf(p.mass, 1, 1) + 
           ", " + p.compositionType, 10, 480 + i*15);
    }
  }
  
  // Controls
  text("=== CONTROLS ===", 10, 520);
  text("Click: Add particle", 10, 535);
  text("Drag: Add dust trail", 10, 550);
  text("R: Reset simulation", 10, 565);
  text("P: Pause/Resume", 10, 580);
  
  // Current simulation time and performance
  text("Simulation Time: " + nf(frameCount * 0.01, 1, 2) + " Myr", 10, 605);
  text("FPS: " + nf(frameRate, 1, 1), 10, 620);
  
  // Disk properties on the right side
  if (showGas) {
    textAlign(RIGHT);
    text("=== DISK PROPERTIES ===", width - 10, 20);
    text("Inner Radius: " + diskInnerRadius, width - 10, 40);
    text("Outer Radius: " + diskOuterRadius, width - 10, 55);
    text("Gas Resolution: " + gasResolution + "x" + gasResolution, width - 10, 70);
    
    // Current mouse position gas properties
    if (mouseX > 0 && mouseY > 0) {
      GasCell nearestGas = findNearestGasCell(new PVector(mouseX, mouseY));
      if (nearestGas != null) {
        text("=== CURSOR GAS PROPERTIES ===", width - 10, 95);
        text("Density: " + nf(nearestGas.density, 1, 3), width - 10, 110);
        text("Pressure: " + nf(nearestGas.pressure, 1, 1), width - 10, 125);
        text("Temperature: " + nf(nearestGas.temperature, 1, 0) + "K", width - 10, 140);
        text("Velocity: " + nf(nearestGas.vel.mag(), 1, 2), width - 10, 155);
      }
    }
    
    // Accretion information
    if (accretionStreams.size() > 0) {
      text("=== ACTIVE ACCRETION ===", width - 10, 180);
      text("Active Streams: " + accretionStreams.size(), width - 10, 195);
      
      // Count planets with circumplanetary disks
      int disksCount = 0;
      for (Particle p : particles) {
        if (p.circumDisk != null) disksCount++;
      }
      text("Circumplanetary Disks: " + disksCount, width - 10, 210);
    }
  }
  
  textAlign(LEFT);
}
