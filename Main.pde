// Enhanced Planet Accretion Simulation
// Main file - Contains global variables, setup, draw, and helper functions

ArrayList<Particle> particles;
ArrayList<Particle> dust;
ArrayList<Particle> planetesimals;
ArrayList<GasCell> gasCells;
ArrayList<AccretionStream> accretionStreams;
ArrayList<SpiralWave> spiralWaves;
Star centralStar;

// Physics constants (scaled for simulation)
float G = 0.5;
float dustSize = 2;
float minAccretionDistance = 8;
boolean showTrails = false;
boolean showForces = false;
boolean showGas = true;
boolean showMagnetic = false;
boolean showTemperature = false;
boolean showCoagulation = true;
boolean showAccretionStreams = true;
boolean showSpiralWaves = true;
boolean showCircumplanetaryDisks = true;

// New physics parameters
float gasViscosity = 0.01;
float magneticStrength = 0.3;
float dustGasCoupling = 0.05;
float temperatureGradient = 0.002;
float pressureScale = 100.0;
float hillSphereRadius = 25.0;
float stickingVelocity = 2.0;
float coagulationEfficiency = 0.8;
float planetesimalFormationThreshold = 5.0;
float gasAccretionRate = 0.001;
float diskPlanetCouplingStrength = 0.1;

// Disk parameters
float diskInnerRadius = 80;
float diskOuterRadius = 450;
int gasResolution = 40;

void setup() {
  size(1200, 800);
  
  particles = new ArrayList<Particle>();
  dust = new ArrayList<Particle>();
  planetesimals = new ArrayList<Particle>();
  gasCells = new ArrayList<GasCell>();
  accretionStreams = new ArrayList<AccretionStream>();
  spiralWaves = new ArrayList<SpiralWave>();
  
  // Create central star
  centralStar = new Star(width/2, height/2, 50);
  
  // Initialize gas disk
  createGasDisk();
  
  // Create initial dust cloud with realistic distribution
  createRealisticDustCloud(400);
  
  // Create seed planetesimals with orbital resonances
  createSeedPlanetesimals(8);
  
  // Create initial planets
  createInitialPlanets(3);
  
  println("=== ENHANCED CONTROLS ===");
  println("SPACE: Toggle trails");
  println("F: Toggle force vectors");
  println("G: Toggle gas visualization");
  println("M: Toggle magnetic field lines");
  println("T: Toggle temperature visualization");
  println("C: Toggle coagulation visualization");
  println("A: Toggle accretion streams");
  println("W: Toggle spiral waves");
  println("D: Toggle circumplanetary disks");
  println("R: Reset simulation");
  println("P: Pause/Resume");
  println("S: Save state");
  println("X: Clear trails");
  println("Z: Add random planetesimal");
  println("1-3: Adjust coagulation efficiency");
  println("4-6: Adjust gas accretion rate");
  println("7-9: Adjust disk-planet coupling");
  println("Mouse: Add particles");
}

void draw() {
  // Dynamic background based on temperature if enabled
  if (showTemperature) {
    drawTemperatureBackground();
  } else if (showTrails) {
    fill(5, 5, 15, 30);
    rect(0, 0, width, height);
  } else {
    background(5, 5, 15);
  }
  
  // Update all physics systems
  updateGasDynamics();
  updateMagneticField();
  updateGravityAndOrbits();
  updateDustGasInteractions();
  updateCoagulationProcess();
  updatePlanetesimalDynamics();
  updateGasAccretionFlows();
  updateDiskPlanetInteractions();
  checkEnhancedAccretion();
  updateParticles();
  
  // Draw everything in correct order
  if (showGas) drawGasDisk();
  if (showSpiralWaves) drawSpiralWaves();
  if (showMagnetic) drawMagneticField();
  if (showAccretionStreams) drawAccretionStreams();
  centralStar.display();
  if (showCircumplanetaryDisks) drawCircumplanetaryDisks();
  drawDust();
  drawPlanetesimals();
  drawParticles();
  if (showCoagulation) drawCoagulationEvents();
  drawUI();
}

// PVector cross product helper method
PVector cross(PVector v1, PVector v2) {
  return new PVector(0, 0, v1.x * v2.y - v1.y * v2.x);
}

// Global temperature calculation function
float calculateTemperature(PVector position) {
  float distanceFromStar = PVector.dist(position, centralStar.pos);
  return 1000 * pow(100.0 / distanceFromStar, 0.5); // Realistic temperature gradient
}
