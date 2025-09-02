// SpiralWave class - Density waves in the disk caused by planets

class SpiralWave {
  float amplitude, wavelength, propagationSpeed;
  float phase, armNumber;
  PVector center;
  color waveColor;
  
  SpiralWave(PVector c, float amp, float wl, int arms) {
    center = c.copy();
    amplitude = amp;
    wavelength = wl;
    propagationSpeed = 0.02;
    armNumber = arms;
    phase = 0;
    waveColor = color(255, 255, 150, 60);
  }
  
  void update() {
    phase += propagationSpeed;
  }
  
  void display() {
    stroke(waveColor);
    strokeWeight(1);
    noFill();
    
    for (int arm = 0; arm < armNumber; arm++) {
      beginShape();
      for (float r = diskInnerRadius; r < diskOuterRadius; r += 5) {
        float baseAngle = arm * TWO_PI / armNumber;
        float spiralAngle = baseAngle + r * 0.01 + phase;
        float waveAmplitude = amplitude * sin(r / wavelength + phase);
        
        float x = center.x + (r + waveAmplitude) * cos(spiralAngle);
        float y = center.y + (r + waveAmplitude) * sin(spiralAngle);
        
        vertex(x, y);
      }
      endShape();
    }
  }
}
