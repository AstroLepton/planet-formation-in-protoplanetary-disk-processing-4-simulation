// Star class - The central star of the system

class Star {
  PVector pos;
  float mass, radius, luminosity;
  color col;
  
  Star(float x, float y, float m) {
    pos = new PVector(x, y);
    mass = m;
    radius = sqrt(mass) * 2;
    luminosity = pow(mass, 3.5); // Mass-luminosity relation
    col = color(255, 255, 200);
  }
  
  void display() {
    // Dynamic stellar wind visualization
    for (int i = 0; i < 8; i++) {
      float windIntensity = 20 + sin(frameCount * 0.1 + i) * 5;
      fill(255, 255, 150, windIntensity);
      ellipse(pos.x, pos.y, radius*2 + i*6, radius*2 + i*6);
    }
    
    // Core with realistic temperature gradient
    fill(col);
    ellipse(pos.x, pos.y, radius*2, radius*2);
    
    // Central hot spot
    fill(255, 255, 255, 200);
    ellipse(pos.x, pos.y, radius, radius);
  }
}
