// Creation.pde - Functions for creating initial conditions and new objects

void createGasDisk() {
  gasCells.clear();
  // Create structured gas disk with radial and azimuthal variations
  for (int i = 0; i < gasResolution; i++) {
    for (int j = 0; j < gasResolution; j++) {
      float angle = map(j, 0, gasResolution-1, 0, TWO_PI);
      float radius = map(i, 0, gasResolution-1, diskInnerRadius, diskOuterRadius);
      
      float x = centralStar.pos.x + cos(angle) * radius;
      float y = centralStar.pos.y + sin(angle) * radius;
      
      // Density profile: higher near star, with spiral structure
      float baseDensity = pow(100.0 / radius, 1.5);
      float spiralDensity = 1 + 0.3 * sin(2 * angle + radius * 0.01);
      float density = baseDensity * spiralDensity * random(0.8, 1.2);
      
      GasCell gas = new GasCell(x, y, density);
      
      // Initial orbital velocity for gas
      float orbitalSpeed = sqrt(G * centralStar.mass / radius) * 0.95;
      gas.vel = new PVector(-sin(angle), cos(angle));
      gas.vel.mult(orbitalSpeed * 0.1); // Scale down for stability
      
      gasCells.add(gas);
    }
  }
}

void createRealisticDustCloud(int numDust) {
  for (int i = 0; i < numDust; i++) {
    // Power-law radial distribution
    float radius = pow(random(0.001, 1), -0.5) * (diskOuterRadius - diskInnerRadius) + diskInnerRadius;
    float angle = random(TWO_PI);
    
    float x = centralStar.pos.x + cos(angle) * radius;
    float y = centralStar.pos.y + sin(angle) * radius;
    
    // Keplerian orbital velocity with random perturbations
    float orbitalSpeed = sqrt(G * centralStar.mass / radius);
    PVector vel = new PVector(-sin(angle), cos(angle));
    vel.mult(orbitalSpeed * random(0.7, 1.1));
    
    // Add turbulent velocity component
    vel.add(PVector.random2D().mult(random(0.2, 0.8)));
    
    Particle d = new Particle(x, y, random(0.3, 1.2));
    d.vel = vel;
    dust.add(d);
  }
}

void createSeedPlanetesimals(int numPlanetesimals) {
  for (int i = 0; i < numPlanetesimals; i++) {
    // Distribute planetesimals with orbital resonances and random perturbations
    float baseRadius = 100 + i * 40 + random(-20, 20);
    float angle = random(TWO_PI);
    
    float x = centralStar.pos.x + cos(angle) * baseRadius;
    float y = centralStar.pos.y + sin(angle) * baseRadius;
    
    // Keplerian velocity with eccentricity
    float orbitalSpeed = sqrt(G * centralStar.mass / baseRadius);
    PVector vel = new PVector(-sin(angle), cos(angle));
    vel.mult(orbitalSpeed * random(0.9, 1.1));
    
    // Add random eccentricity
    vel.add(PVector.random2D().mult(random(0.3, 1.0)));
    
    Particle p = new Particle(x, y, random(planetesimalFormationThreshold, 15));
    p.vel = vel;
    planetesimals.add(p);
  }
}

void createInitialPlanets(int numPlanets) {
  float[] planetDistances = {200, 320, 420}; // Specific orbital distances
  
  for (int i = 0; i < numPlanets && i < planetDistances.length; i++) {
    float distance = planetDistances[i];
    float angle = i * TWO_PI / numPlanets + random(-0.5, 0.5);
    
    float x = centralStar.pos.x + cos(angle) * distance;
    float y = centralStar.pos.y + sin(angle) * distance;
    
    // Keplerian orbital velocity
    float orbitalSpeed = sqrt(G * centralStar.mass / distance);
    PVector vel = new PVector(-sin(angle), cos(angle));
    vel.mult(orbitalSpeed * random(0.95, 1.05));
    
    Particle planet = new Particle(x, y, random(15, 30));
    planet.vel = vel;
    particles.add(planet);
    
    println("Created planet " + (i+1) + ": Mass=" + nf(planet.mass, 1, 1) + 
            ", Distance=" + nf(distance, 1, 0));
  }
}

void createAccretionStream(Particle planet) {
  // Find nearby gas to create accretion stream
  float searchRadius = 80;
  GasCell sourceGas = null;
  float maxDensity = 0;
  
  for (GasCell gas : gasCells) {
    float distance = PVector.dist(planet.pos, gas.pos);
    if (distance < searchRadius && gas.density > maxDensity) {
      maxDensity = gas.density;
      sourceGas = gas;
    }
  }
  
  if (sourceGas != null) {
    float flowRate = sourceGas.density * planet.mass * 0.001;
    AccretionStream stream = new AccretionStream(sourceGas.pos, planet, flowRate);
    accretionStreams.add(stream);
  }
}

void createSpiralWave(Particle planet) {
  float wavelength = PVector.dist(planet.pos, centralStar.pos) * 0.5;
  float amplitude = planet.mass * 2;
  int armNumber = 2; // Typical for planet-induced waves
  
  SpiralWave wave = new SpiralWave(centralStar.pos, amplitude, wavelength, armNumber);
  spiralWaves.add(wave);
}

void createDiskGap(Particle planet) {
  float gapRadius = pow(planet.mass / centralStar.mass, 1.0/3.0) * 
                   PVector.dist(planet.pos, centralStar.pos) * 2;
  
  // Reduce gas density around planet
  for (GasCell gas : gasCells) {
    float distanceToPlanet = PVector.dist(gas.pos, planet.pos);
    if (distanceToPlanet < gapRadius) {
      float gapFactor = distanceToPlanet / gapRadius;
      gas.density *= (0.1 + 0.9 * gapFactor * gapFactor);
    }
  }
  
  // Clear dust from gap region
  for (int i = dust.size()-1; i >= 0; i--) {
    Particle dustParticle = dust.get(i);
    float distanceToPlanet = PVector.dist(dustParticle.pos, planet.pos);
    if (distanceToPlanet < gapRadius * 0.8) {
      if (random(1) < 0.01) { // Small chance to remove dust
        dust.remove(i);
      }
    }
  }
}
