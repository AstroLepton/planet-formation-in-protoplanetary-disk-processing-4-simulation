// AccretionStream class - Gas flows toward massive planets

class AccretionStream {
  PVector startPos, endPos;
  Particle targetPlanet;
  float flowRate, temperature;
  ArrayList<PVector> streamPath;
  color streamColor;
  float creationTime;
  
  AccretionStream(PVector start, Particle target, float rate) {
    startPos = start.copy();
    endPos = target.pos.copy();
    targetPlanet = target;
    flowRate = rate;
    temperature = calculateTemperature(startPos);
    streamPath = new ArrayList<PVector>();
    createStreamPath();
    streamColor = color(255, 150, 100, 180);
    creationTime = frameCount;
  }
  
  void createStreamPath() {
    streamPath.clear();
    int pathSteps = 20;
    
    for (int i = 0; i <= pathSteps; i++) {
      float t = i / (float)pathSteps;
      PVector pathPoint = PVector.lerp(startPos, endPos, t);
      
      // Add gravitational curvature
      PVector toStar = PVector.sub(centralStar.pos, pathPoint);
      float starDist = toStar.mag();
      toStar.normalize();
      toStar.mult(G * centralStar.mass / (starDist * starDist) * t * t * 100);
      pathPoint.add(toStar);
      
      streamPath.add(pathPoint);
    }
  }
  
  void update() {
    endPos = targetPlanet.pos.copy();
    createStreamPath();
    
    // Stream fades over time
    float age = frameCount - creationTime;
    float alpha = map(age, 0, 60, 180, 0);
    alpha = constrain(alpha, 0, 180);
    streamColor = color(255, 150, 100, alpha);
  }
  
  void display() {
    stroke(streamColor);
    strokeWeight(2 + flowRate * 50);
    noFill();
    
    beginShape();
    for (PVector point : streamPath) {
      vertex(point.x, point.y);
    }
    endShape();
    
    // Draw flow direction indicators
    for (int i = 0; i < streamPath.size()-1; i += 3) {
      PVector p1 = streamPath.get(i);
      PVector p2 = streamPath.get(i+1);
      PVector dir = PVector.sub(p2, p1);
      dir.normalize();
      dir.mult(5);
      
      stroke(streamColor);
      strokeWeight(1);
      line(p1.x, p1.y, p1.x + dir.x, p1.y + dir.y);
    }
  }
  
  boolean isExpired() {
    return frameCount - creationTime > 60;
  }
}
