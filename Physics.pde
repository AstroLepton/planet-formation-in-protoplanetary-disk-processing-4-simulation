// Physics.pde - All physics calculations and updates

void updateGasDynamics() {
  for (GasCell gas : gasCells) {
    gas.update();
  }
}

void updateMagneticField() {
  if (!showMagnetic) return;
  
  // Apply magnetic forces to charged particles (simplified)
  for (Particle p : particles) {
    if (p.mass > 5) { // Only larger particles affected significantly
      PVector magneticForce = calculateMagneticForce(p);
      p.applyForce(magneticForce);
    }
  }
}

PVector calculateMagneticForce(Particle p) {
  // Simplified magnetic field: dipolar field from star
  PVector r = PVector.sub(p.pos, centralStar.pos);
  float distance = r.mag();
  
  // Magnetic field strength decreases with distance
  float fieldStrength = magneticStrength * centralStar.mass / (distance * distance);
  
  // Force perpendicular to velocity (Lorentz force simplified)
  PVector magneticForce = new PVector(-p.vel.y, p.vel.x);
  magneticForce.mult(fieldStrength * 0.1);
  
  return magneticForce;
}

void updateGravityAndOrbits() {
  // Enhanced gravitational interactions with relativistic corrections for close approaches
  
  // Star attracts all particles
  for (Particle p : particles) {
    PVector starForce = calculateEnhancedGravity(centralStar.pos, centralStar.mass, p.pos, p.mass);
    p.applyForce(starForce);
  }
  
  for (Particle d : dust) {
    PVector starForce = calculateEnhancedGravity(centralStar.pos, centralStar.mass, d.pos, d.mass);
    d.applyForce(starForce);
  }
  
  for (Particle p : planetesimals) {
    PVector starForce = calculateEnhancedGravity(centralStar.pos, centralStar.mass, p.pos, p.mass);
    p.applyForce(starForce);
  }
  
  // N-body gravitational interactions between particles
  for (int i = 0; i < particles.size(); i++) {
    for (int j = i + 1; j < particles.size(); j++) {
      Particle p1 = particles.get(i);
      Particle p2 = particles.get(j);
      
      PVector force = calculateEnhancedGravity(p1.pos, p1.mass, p2.pos, p2.mass);
      p1.applyForce(force);
      p2.applyForce(PVector.mult(force, -1));
      
      // Check for Hill sphere interactions
      checkHillSphereInteraction(p1, p2);
    }
  }
}

PVector calculateEnhancedGravity(PVector pos1, float mass1, PVector pos2, float mass2) {
  PVector force = PVector.sub(pos1, pos2);
  float distance = force.mag();
  distance = constrain(distance, 3, 200); // Prevent singularities
  
  float strength = (G * mass1 * mass2) / (distance * distance);
  
  // Add relativistic correction for very close approaches
  if (distance < 20) {
    float relativisticCorrection = 1 + 3 * G * (mass1 + mass2) / (distance * 299792458 * 299792458 * 0.0001);
    strength *= relativisticCorrection;
  }
  
  force.normalize();
  force.mult(strength);
  
  return force;
}

void checkHillSphereInteraction(Particle p1, Particle p2) {
  float distance = PVector.dist(p1.pos, p2.pos);
  float hillRadius1 = pow(p1.mass / (3 * centralStar.mass), 1.0/3.0) * 
                     PVector.dist(p1.pos, centralStar.pos);
  float hillRadius2 = pow(p2.mass / (3 * centralStar.mass), 1.0/3.0) * 
                     PVector.dist(p2.pos, centralStar.pos);
  
  if (distance < max(hillRadius1, hillRadius2)) {
    p1.isInHillSphere = true;
    p2.isInHillSphere = true;
  }
}

void updateDustGasInteractions() {
  // Aerodynamic drag and dust-gas coupling
  for (Particle d : dust) {
    GasCell nearestGas = findNearestGasCell(d.pos);
    if (nearestGas != null) {
      applyAerodynamicDrag(d, nearestGas);
    }
  }
  
  // Gas accretion onto particles
  for (Particle p : particles) {
    if (p.mass > 10) { // Only massive particles can accrete gas significantly
      accreteGas(p);
    }
  }
}

GasCell findNearestGasCell(PVector pos) {
  GasCell nearest = null;
  float minDistance = Float.MAX_VALUE;
  
  for (GasCell gas : gasCells) {
    float distance = PVector.dist(pos, gas.pos);
    if (distance < minDistance) {
      minDistance = distance;
      nearest = gas;
    }
  }
  
  return minDistance < 20 ? nearest : null;
}

void applyAerodynamicDrag(Particle particle, GasCell gas) {
  // Relative velocity between particle and gas
  PVector relativeVel = PVector.sub(particle.vel, gas.vel);
  float relativeSpeed = relativeVel.mag();
  
  if (relativeSpeed > 0.1) {
    // Drag force proportional to relative velocity squared
    float dragCoeff = dustGasCoupling * gas.density * particle.radius * particle.radius;
    PVector dragForce = PVector.mult(relativeVel, -dragCoeff * relativeSpeed);
    
    particle.applyForce(dragForce);
    
    // Dust settling toward midplane (simplified)
    float settlingForce = -0.01 * (particle.pos.y - centralStar.pos.y);
    particle.applyForce(new PVector(0, settlingForce));
  }
}

void accreteGas(Particle planet) {
  float accretionRadius = planet.radius * 2;
  float accretedMass = 0;
  
  for (int i = gasCells.size()-1; i >= 0; i--) {
    GasCell gas = gasCells.get(i);
    float distance = PVector.dist(planet.pos, gas.pos);
    
    if (distance < accretionRadius) {
      // Accrete a fraction of the gas
      float accretionRate = 0.001 * planet.mass / (distance + 1);
      accretedMass += gas.density * accretionRate;
      gas.density *= (1 - accretionRate);
      
      if (gas.density < 0.01) {
        gasCells.remove(i);
      }
    }
  }
  
  if (accretedMass > 0) {
    planet.mass += accretedMass;
    planet.radius = pow(planet.mass * 0.75 / PI, 1.0/3.0) * 2;
    planet.col = planet.getParticleColor();
  }
}

void updateCoagulationProcess() {
  // Dust coagulation and growth into larger aggregates
  for (int i = dust.size()-1; i >= 0; i--) {
    for (int j = i-1; j >= 0; j--) {
      Particle d1 = dust.get(i);
      Particle d2 = dust.get(j);
      
      float distance = PVector.dist(d1.pos, d2.pos);
      float collisionRadius = (d1.radius + d2.radius) * 1.2;
      
      if (distance < collisionRadius) {
        PVector relativeVel = PVector.sub(d1.vel, d2.vel);
        float relativeSpeed = relativeVel.mag();
        
        // Coagulation probability based on sticking velocity
        float stickingProb = exp(-relativeSpeed / stickingVelocity) * 
                           min(d1.stickinessCoeff, d2.stickinessCoeff);
        
        if (random(1) < stickingProb) {
          // Successful coagulation
          Particle newAggregate = coagulateParticles(d1, d2);
          dust.remove(max(i, j));
          dust.remove(min(i, j));
          
          // Check if aggregate is large enough to become planetesimal
          if (newAggregate.mass > planetesimalFormationThreshold) {
            planetesimals.add(newAggregate);
          } else {
            dust.add(newAggregate);
          }
          break;
        } else {
          // Bouncing collision - no sticking
          bounceCollision(d1, d2);
        }
      }
    }
  }
}

Particle coagulateParticles(Particle p1, Particle p2) {
  // Create new aggregate particle
  PVector newPos = PVector.lerp(p1.pos, p2.pos, p1.mass / (p1.mass + p2.mass));
  Particle aggregate = new Particle(newPos.x, newPos.y, p1.mass + p2.mass);
  
  // Conservation of momentum
  PVector totalMomentum = PVector.add(
    PVector.mult(p1.vel, p1.mass),
    PVector.mult(p2.vel, p2.mass)
  );
  aggregate.vel = PVector.div(totalMomentum, aggregate.mass);
  
  // Porosity increases with coagulation (fluffy aggregates)
  aggregate.porosityFactor = (p1.porosityFactor + p2.porosityFactor) * 1.1;
  aggregate.radius = aggregate.radius * sqrt(aggregate.porosityFactor);
  
  // Update density accounting for porosity
  aggregate.density = aggregate.mass / (4.0/3.0 * PI * pow(aggregate.radius, 3));
  
  // Mixed composition
  if (p1.compositionType.equals(p2.compositionType)) {
    aggregate.compositionType = p1.compositionType;
  } else {
    aggregate.compositionType = "mixed";
  }
  
  aggregate.col = aggregate.getParticleColor();
  aggregate.lastCollisionTime = frameCount;
  
  return aggregate;
}

void updatePlanetesimalDynamics() {
  // Planetesimal orbital evolution and gravitational interactions
  
  // Planetesimal-planetesimal interactions
  for (int i = 0; i < planetesimals.size(); i++) {
    for (int j = i + 1; j < planetesimals.size(); j++) {
      Particle p1 = planetesimals.get(i);
      Particle p2 = planetesimals.get(j);
      
      PVector force = calculateEnhancedGravity(p1.pos, p1.mass, p2.pos, p2.mass);
      p1.applyForce(force);
      p2.applyForce(PVector.mult(force, -1));
      
      // Check for planetesimal accretion
      float distance = PVector.dist(p1.pos, p2.pos);
      if (distance < (p1.radius + p2.radius) * 0.9) {
        checkPlanetesimalAccretion(p1, p2, i, j);
      }
    }
  }
  
  // Planetesimal-planet interactions
  for (Particle planet : particles) {
    for (int i = planetesimals.size()-1; i >= 0; i--) {
      Particle planetesimal = planetesimals.get(i);
      
      PVector force = calculateEnhancedGravity(planet.pos, planet.mass, 
                                             planetesimal.pos, planetesimal.mass);
      planetesimal.applyForce(force);
      planet.applyForce(PVector.mult(force, -1));
      
      // Check for capture into Hill sphere
      float hillRadius = pow(planet.mass / (3 * centralStar.mass), 1.0/3.0) * 
                        PVector.dist(planet.pos, centralStar.pos);
      float distance = PVector.dist(planet.pos, planetesimal.pos);
      
      if (distance < hillRadius * 0.7) {
        // Planetesimal captured - becomes part of planet or circumplanetary disk
        if (distance < planet.radius * 2) {
          // Direct accretion
          planet.mass += planetesimal.mass;
          planet.radius = pow(planet.mass * 0.75 / PI, 1.0/3.0) * 2;
          planetesimals.remove(i);
        } else if (planet.circumDisk != null) {
          // Add to circumplanetary disk
          planet.circumDisk.addParticle(planetesimal);
          planetesimals.remove(i);
        }
      }
    }
  }
}

void checkPlanetesimalAccretion(Particle p1, Particle p2, int i, int j) {
  PVector relativeVel = PVector.sub(p1.vel, p2.vel);
  float impactSpeed = relativeVel.mag();
  float escapeVelocity = sqrt(2 * G * (p1.mass + p2.mass) / (p1.radius + p2.radius));
  
  if (impactSpeed < escapeVelocity * 0.8) {
    // Successful accretion
    if (p1.mass >= p2.mass) {
      mergeParticlesEnhanced(p1, p2);
      planetesimals.remove(j);
    } else {
      mergeParticlesEnhanced(p2, p1);
      planetesimals.remove(i);
    }
  } else {
    // Collision may result in fragmentation or bouncing
    if (impactSpeed > escapeVelocity * 2.0) {
      // Fragmentation - break apart into smaller pieces
      fragmentPlanetesimals(p1, p2, i, j);
    } else {
      // Bouncing collision
      bounceCollision(p1, p2);
    }
  }
}

void fragmentPlanetesimals(Particle p1, Particle p2, int i, int j) {
  // Create fragments from planetesimal collision
  PVector collisionCenter = PVector.lerp(p1.pos, p2.pos, p1.mass / (p1.mass + p2.mass));
  
  int numFragments = (int)random(4, 10);
  float totalMass = p1.mass + p2.mass;
  
  planetesimals.remove(max(i, j));
  planetesimals.remove(min(i, j));
  
  for (int k = 0; k < numFragments; k++) {
    float fragmentMass = totalMass * random(0.03, 0.3);
    PVector fragmentPos = PVector.add(collisionCenter, PVector.random2D().mult(random(8, 20)));
    
    Particle fragment = new Particle(fragmentPos.x, fragmentPos.y, fragmentMass);
    fragment.vel = PVector.add(
      PVector.lerp(p1.vel, p2.vel, 0.5),
      PVector.random2D().mult(random(2, 5))
    );
    fragment.temperature += 150; // Heating from impact
    
    if (fragmentMass < planetesimalFormationThreshold) {
      dust.add(fragment);
    } else {
      planetesimals.add(fragment);
    }
  }
}

void updateGasAccretionFlows() {
  // Create and update gas accretion streams toward massive planets
  
  // Remove expired streams
  for (int i = accretionStreams.size()-1; i >= 0; i--) {
    AccretionStream stream = accretionStreams.get(i);
    if (stream.isExpired()) {
      accretionStreams.remove(i);
    } else {
      stream.update();
    }
  }
  
  // Create new accretion streams for massive planets
  for (Particle planet : particles) {
    if (planet.mass > 15 && random(1) < 0.05) { // 5% chance per frame for massive planets
      createAccretionStream(planet);
    }
  }
  
  // Gas accretion onto planets
  for (Particle planet : particles) {
    if (planet.mass > 8) {
      float accretionRadius = planet.radius * 3;
      float totalAccretedMass = 0;
      
      for (int i = gasCells.size()-1; i >= 0; i--) {
        GasCell gas = gasCells.get(i);
        float distance = PVector.dist(planet.pos, gas.pos);
        
        if (distance < accretionRadius) {
          // Gas accretion rate depends on planet mass and local gas density
          float accretionRate = gasAccretionRate * planet.mass * gas.density / (distance + 1);
          float accretedMass = gas.density * accretionRate;
          
          totalAccretedMass += accretedMass;
          gas.density *= (1 - accretionRate);
          
          // Create accretion heating
          planet.temperature += accretedMass * 10;
          
          if (gas.density < 0.005) {
            gasCells.remove(i);
          }
        }
      }
      
      if (totalAccretedMass > 0) {
        planet.mass += totalAccretedMass;
        planet.radius = pow(planet.mass * 0.75 / PI, 1.0/3.0) * 2;
        planet.col = planet.getParticleColor();
        
        // Create circumplanetary disk if planet becomes massive enough
        if (planet.mass > 20 && planet.circumDisk == null) {
          planet.circumDisk = new CircumplanetaryDisk(planet);
        }
      }
    }
  }
}

void updateDiskPlanetInteractions() {
  // Disk-planet gravitational interactions creating gaps and spiral waves
  
  for (Particle planet : particles) {
    if (planet.mass > 12) { // Only massive planets significantly affect the disk
      
      // Create gap in disk around planet
      createDiskGap(planet);
      
      // Generate spiral density waves
      if (random(1) < 0.02) { // Occasional wave generation
        createSpiralWave(planet);
      }
      
      // Torque exchange between planet and disk
      applyDiskTorque(planet);
    }
  }
  
  // Update existing spiral waves
  for (int i = spiralWaves.size()-1; i >= 0; i--) {
    SpiralWave wave = spiralWaves.get(i);
    wave.update();
    
    // Remove old waves
    if (frameCount % 600 == 0 && random(1) < 0.3) {
      spiralWaves.remove(i);
    }
  }
}

void applyDiskTorque(Particle planet) {
  // Simplified torque from disk-planet interaction
  float planetOrbitalRadius = PVector.dist(planet.pos, centralStar.pos);
  
  // Inner disk exerts positive torque, outer disk negative torque
  float netTorque = 0;
  
  for (GasCell gas : gasCells) {
    float gasRadius = PVector.dist(gas.pos, centralStar.pos);
    float distance = PVector.dist(gas.pos, planet.pos);
    
    if (distance < 50) { // Only nearby gas contributes significantly
      float torqueContribution = gas.density * planet.mass / (distance * distance);
      
      if (gasRadius < planetOrbitalRadius) {
        netTorque += torqueContribution; // Inner disk speeds up planet
      } else {
        netTorque -= torqueContribution; // Outer disk slows down planet
      }
    }
  }
  
  // Apply torque as tangential force
  PVector toPlanet = PVector.sub(planet.pos, centralStar.pos);
  PVector tangentialForce = new PVector(-toPlanet.y, toPlanet.x);
  tangentialForce.normalize();
  tangentialForce.mult(netTorque * diskPlanetCouplingStrength * 0.0001);
  
  planet.applyForce(tangentialForce);
}

void checkEnhancedAccretion() {
  // Enhanced accretion with collision physics and fragmentation
  
  // Particle-particle collisions with impact velocity consideration
  for (int i = particles.size()-1; i >= 0; i--) {
    for (int j = i-1; j >= 0; j--) {
      Particle p1 = particles.get(i);
      Particle p2 = particles.get(j);
      
      float distance = PVector.dist(p1.pos, p2.pos);
      float collisionDistance = (p1.radius + p2.radius) * 0.9;
      
      if (distance < collisionDistance) {
        PVector relativeVel = PVector.sub(p1.vel, p2.vel);
        float impactSpeed = relativeVel.mag();
        
        // Determine collision outcome based on impact speed and masses
        float escapeVelocity = sqrt(2 * G * (p1.mass + p2.mass) / collisionDistance);
        
        if (impactSpeed < escapeVelocity * 0.5) {
          // Accretion - merge particles
          if (p1.mass >= p2.mass) {
            mergeParticlesEnhanced(p1, p2);
            particles.remove(j);
          } else {
            mergeParticlesEnhanced(p2, p1);
            particles.remove(i);
            break;
          }
        } else if (impactSpeed < escapeVelocity * 1.5) {
          // Bouncing collision
          bounceCollision(p1, p2);
        } else {
          // Catastrophic fragmentation
          fragmentParticles(p1, p2, i, j);
        }
      }
    }
  }
  
  // Enhanced dust accretion with sticking probability
  for (Particle p : particles) {
    for (int i = dust.size()-1; i >= 0; i--) {
      Particle d = dust.get(i);
      float distance = PVector.dist(p.pos, d.pos);
      
      if (distance < p.radius + 3) {
        // Sticking probability based on relative velocity
        PVector relativeVel = PVector.sub(p.vel, d.vel);
        float stickingProb = 1.0 / (1 + relativeVel.magSq() * 10);
        
        if (random(1) < stickingProb) {
          p.mass += d.mass;
          p.radius = pow(p.mass * 0.75 / PI, 1.0/3.0) * 2;
          p.col = p.getParticleColor();
          dust.remove(i);
        }
      }
    }
  }
}

void mergeParticlesEnhanced(Particle larger, Particle smaller) {
  // Enhanced merger with energy dissipation and composition mixing
  
  // Conservation of momentum
  PVector totalMomentum = PVector.add(
    PVector.mult(larger.vel, larger.mass),
    PVector.mult(smaller.vel, smaller.mass)
  );
  
  float totalMass = larger.mass + smaller.mass;
  larger.vel = PVector.div(totalMomentum, totalMass);
  
  // Energy dissipation during collision
  larger.vel.mult(0.95); // Some energy lost to heat/deformation
  
  // Update properties
  larger.mass = totalMass;
  larger.radius = pow(larger.mass * 0.75 / PI, 1.0/3.0) * 2;
  larger.density = larger.mass / (4.0/3.0 * PI * pow(larger.radius, 3));
  
  // Temperature increase from collision energy
  larger.temperature += 50;
  larger.col = larger.getParticleColor();
  
  // Position adjustment (center of mass)
  PVector centerOfMass = PVector.add(
    PVector.mult(larger.pos, larger.mass - smaller.mass),
    PVector.mult(smaller.pos, smaller.mass)
  );
  larger.pos = PVector.div(centerOfMass, larger.mass);
}

void bounceCollision(Particle p1, Particle p2) {
  // Elastic collision with some energy loss
  PVector relativePos = PVector.sub(p2.pos, p1.pos);
  relativePos.normalize();
  
  PVector relativeVel = PVector.sub(p1.vel, p2.vel);
  float velocityAlongNormal = PVector.dot(relativeVel, relativePos);
  
  if (velocityAlongNormal > 0) return; // Objects separating
  
  float restitution = 0.6; // Coefficient of restitution
  float impulse = -(1 + restitution) * velocityAlongNormal / (1/p1.mass + 1/p2.mass);
  
  PVector impulseVector = PVector.mult(relativePos, impulse);
  p1.vel.add(PVector.div(impulseVector, p1.mass));
  p2.vel.sub(PVector.div(impulseVector, p2.mass));
}

void fragmentParticles(Particle p1, Particle p2, int i, int j) {
  // Create fragments from catastrophic collision
  PVector collisionCenter = PVector.lerp(p1.pos, p2.pos, p1.mass / (p1.mass + p2.mass));
  
  int numFragments = (int)random(3, 8);
  float totalMass = p1.mass + p2.mass;
  
  particles.remove(max(i, j));
  particles.remove(min(i, j));
  
  for (int k = 0; k < numFragments; k++) {
    float fragmentMass = totalMass * random(0.05, 0.4);
    PVector fragmentPos = PVector.add(collisionCenter, PVector.random2D().mult(random(5, 15)));
    
    Particle fragment = new Particle(fragmentPos.x, fragmentPos.y, fragmentMass);
    fragment.vel = PVector.add(
      PVector.lerp(p1.vel, p2.vel, 0.5),
      PVector.random2D().mult(random(1, 4))
    );
    fragment.temperature += 200; // Heating from impact
    
    if (fragmentMass < 2) {
      dust.add(fragment);
    } else {
      particles.add(fragment);
    }
  }
}

void updateParticles() {
  for (Particle p : particles) {
    p.update();
  }
  for (Particle d : dust) {
    d.update();
  }
  for (Particle p : planetesimals) {
    p.update();
  }
}
