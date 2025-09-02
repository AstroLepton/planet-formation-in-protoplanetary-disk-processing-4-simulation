// Rendering.pde - All drawing and visualization functions

void drawParticles() {
  for (Particle p : particles) {
    p.display();
  }
}

void drawPlanetesimals() {
  for (Particle p : planetesimals) {
    // Draw planetesimals with distinct appearance
    fill(150, 100, 80, 200);
    stroke(100, 70, 50);
    strokeWeight(1);
    ellipse(p.pos.x, p.pos.y, p.radius*2, p.radius*2);
    
    // Draw composition indicator
    if (p.compositionType.equals("rocky")) {
      fill(200, 100, 50, 100);
      ellipse(p.pos.x, p.pos.y, p.radius, p.radius);
    } else if (p.compositionType.equals("icy")) {
      fill(150, 200, 255, 100);
      ellipse(p.pos.x, p.pos.y, p.radius, p.radius);
    }
  }
}

void drawDust() {
  for (Particle d : dust) {
    // Color dust based on composition and size
    if (d.compositionType.equals("rocky")) {
      fill(150, 100, 80, 180);
    } else if (d.compositionType.equals("icy")) {
      fill(150, 200, 255, 180);
    } else {
      fill(120, 150, 200, 180);
    }
    
    noStroke();
    float displaySize = dustSize + d.mass * 0.5; // Size based on mass
    ellipse(d.pos.x, d.pos.y, displaySize, displaySize);
  }
}

void drawGasDisk() {
  for (GasCell gas : gasCells) {
    gas.display();
  }
}

void drawMagneticField() {
  // Draw simplified magnetic field lines
  stroke(255, 100, 255, 80);
  strokeWeight(1);
  
  for (int i = 0; i < 12; i++) {
    float angle = i * TWO_PI / 12;
    PVector start = PVector.add(centralStar.pos, PVector.fromAngle(angle).mult(60));
    PVector end = PVector.add(centralStar.pos, PVector.fromAngle(angle + PI).mult(300));
    
    // Curved field lines
    for (int j = 0; j < 20; j++) {
      float t = j / 19.0;
      PVector pos = PVector.lerp(start, end, t);
      
      // Add curvature
      float curvature = sin(t * PI) * 30;
      pos.add(PVector.fromAngle(angle + PI/2).mult(curvature));
      
      if (j > 0) {
        PVector prevPos = PVector.lerp(start, end, (j-1) / 19.0);
        float prevCurvature = sin((j-1) / 19.0 * PI) * 30;
        prevPos.add(PVector.fromAngle(angle + PI/2).mult(prevCurvature));
        
        line(prevPos.x, prevPos.y, pos.x, pos.y);
      }
    }
  }
}

void drawTemperatureBackground() {
  // Draw temperature gradient background efficiently
  noStroke();
  int numCircles = 50;
  for (int i = numCircles; i > 0; i--) {
    float radius = i * (width / (float)numCircles);
    float distFromStar = radius;
    float temp = 1000 * pow(100.0 / (distFromStar + 10), 0.5);
    
    float r = map(temp, 100, 1000, 0, 255);
    float g = map(temp, 100, 1000, 0, 150);
    float b = map(temp, 100, 1000, 100, 0);
    
    r = constrain(r, 0, 255);
    g = constrain(g, 0, 255);
    b = constrain(b, 0, 255);
    
    fill(r, g, b, 20);
    ellipse(centralStar.pos.x, centralStar.pos.y, radius * 2, radius * 2);
  }
}

void drawAccretionStreams() {
  for (AccretionStream stream : accretionStreams) {
    stream.display();
  }
}

void drawSpiralWaves() {
  for (SpiralWave wave : spiralWaves) {
    wave.display();
  }
}

void drawCircumplanetaryDisks() {
  for (Particle planet : particles) {
    if (planet.circumDisk != null) {
      planet.circumDisk.display();
    }
  }
}

void drawCoagulationEvents() {
  // Draw visual indicators for recent coagulation events
  for (Particle d : dust) {
    if (frameCount - d.lastCollisionTime < 30) {
      stroke(255, 255, 0, 200);
      strokeWeight(2);
      noFill();
      ellipse(d.pos.x, d.pos.y, d.radius * 4, d.radius * 4);
    }
  }
}
