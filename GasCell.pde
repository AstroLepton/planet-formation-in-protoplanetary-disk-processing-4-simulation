// GasCell class - Represents gas elements in the protoplanetary disk

class GasCell {
  PVector pos, vel;
  float density, pressure, temperature, viscosity;
  color gasColor;
  
  GasCell(float x, float y, float d) {
    pos = new PVector(x, y);
    vel = new PVector(0, 0);
    density = d;
    temperature = 1000 * pow(100.0 / PVector.dist(pos, centralStar.pos), 0.5);
    pressure = density * temperature * pressureScale;
    viscosity = gasViscosity;
    updateColor();
  }
  
  void updateColor() {
    float densityFactor = map(log(density + 1), 0, 3, 0, 255);
    float tempFactor = map(temperature, 100, 1000, 0, 255);
    gasColor = color(tempFactor * 0.3, densityFactor * 0.5, 255 - tempFactor * 0.2, 60);
  }
  
  void update() {
    // Simple gas dynamics - pressure gradient force
    PVector pressureGradient = calculatePressureGradient();
    vel.add(PVector.mult(pressureGradient, 0.01));
    
    // Viscous damping
    vel.mult(1 - viscosity);
    
    // Update position
    pos.add(PVector.mult(vel, 0.1));
    
    // Recalculate properties
    temperature = 1000 * pow(100.0 / PVector.dist(pos, centralStar.pos), 0.5);
    pressure = density * temperature * pressureScale;
    updateColor();
  }
  
  PVector calculatePressureGradient() {
    // Simplified pressure gradient calculation
    float distFromStar = PVector.dist(pos, centralStar.pos);
    PVector gradient = PVector.sub(centralStar.pos, pos);
    gradient.normalize();
    gradient.mult(-pressure / (distFromStar * distFromStar));
    return gradient;
  }
  
  void display() {
    fill(gasColor);
    noStroke();
    ellipse(pos.x, pos.y, 8, 8);
  }
}
