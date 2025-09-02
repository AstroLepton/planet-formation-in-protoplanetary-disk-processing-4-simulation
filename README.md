# Enhanced Planet Accretion Simulation

A sophisticated Processing-based simulation that models the complex physics of planet formation in protoplanetary disks, from dust coagulation to full planetary systems.

## Overview

This simulation recreates the  process of planet formation, starting from microscopic dust particles and evolving into a complete planetary system. It incorporates realistic physics including gravitational interactions, gas dynamics, magnetic fields, and complex collision outcomes.

## Physics Engine

### Core Physics Systems

#### 1. **Gravitational Dynamics**
- **N-body gravitational interactions** between all particles
- **Relativistic corrections** for close approaches
- **Hill sphere mechanics** for planetary influence zones
- **Orbital evolution** with realistic eccentricity and semi-major axis calculations

#### 2. **Gas Disk Physics**
- **Protoplanetary disk** with realistic temperature and density gradients
- **Pressure gradient forces** driving gas motion
- **Viscous evolution** of the gas disk
- **Gas accretion** onto growing planets

#### 3. **Dust Dynamics**
- **Aerodynamic drag** between dust and gas
- **Dust settling** toward the disk midplane
- **Temperature-dependent composition** (rocky, icy, mixed)
- **Porosity evolution** during coagulation

#### 4. **Collision Physics**
- **Sticking probability** based on impact velocity and material properties
- **Coagulation efficiency** for dust growth
- **Fragmentation** for high-velocity impacts
- **Energy dissipation** during collisions
- **Composition mixing** during mergers

#### 5. **Advanced Features**
- **Magnetic field interactions** with charged particles
- **Spiral density waves** from planet-disk interactions
- **Disk gaps** around massive planets
- **Circumplanetary disks** for large bodies
- **Accretion streams** toward massive planets

### Physical Parameters

- **Gravitational constant**: G = 0.5 (scaled for simulation)
- **Gas viscosity**: 0.01
- **Magnetic field strength**: 0.3
- **Dust-gas coupling**: 0.05
- **Temperature gradient**: 0.002
- **Hill sphere radius**: 25.0
- **Sticking velocity**: 2.0
- **Coagulation efficiency**: 0.8

## What You'll See

### 1. **Initial State**
- **Central star** with dynamic stellar wind visualization
- **Gas disk** extending from 80 to 450 units radius
- **Dust cloud** with 400 initial particles
- **Seed planetesimals** with orbital resonances
- **Initial planets** to kickstart the system

### 2. **Evolution Process**

#### **Phase 1: Dust Coagulation (Early)**
- Small dust particles begin sticking together
- Formation of larger aggregates
- Temperature-dependent composition zones emerge
- Rocky materials closer to star, icy materials farther out

#### **Phase 2: Planetesimal Formation (Mid)**
- Aggregates reach critical mass threshold
- Gravitational focusing increases collision rates
- Hill spheres begin to form around larger bodies
- Orbital resonances and gravitational interactions

#### **Phase 3: Planetary Growth (Advanced)**
- Massive bodies create gaps in the gas disk
- Spiral density waves form from disk-planet interactions
- Gas accretion onto growing planets
- Circumplanetary disks develop around large bodies

#### **Phase 4: System Maturation (Late)**
- Stable orbital configurations emerge
- Gas disk dissipation
- Final planetary system architecture

### 3. **Visual Phenomena**

#### **Particle Types & Colors**
- **Dust**: Small, numerous, various colors based on composition
- **Planetesimals**: Medium-sized, growing bodies
- **Planets**: Large, massive bodies with distinct colors
- **Gas cells**: Semi-transparent, temperature/density colored

#### **Dynamic Effects**
- **Trails**: Particle motion history (toggle with SPACE)
- **Force vectors**: Gravitational and other forces (toggle with F)
- **Magnetic field lines**: Field visualization (toggle with M)
- **Temperature gradients**: Heat distribution (toggle with T)
- **Accretion streams**: Gas flow toward planets (toggle with A)
- **Spiral waves**: Density wave patterns (toggle with W)

## Controls

### **Visualization Toggles**
- **SPACE**: Toggle particle trails
- **F**: Toggle force vector display
- **G**: Toggle gas disk visualization
- **M**: Toggle magnetic field lines
- **T**: Toggle temperature visualization
- **C**: Toggle coagulation events
- **A**: Toggle accretion streams
- **W**: Toggle spiral waves
- **D**: Toggle circumplanetary disks

### **Simulation Control**
- **R**: Reset simulation
- **P**: Pause/Resume
- **S**: Save current state
- **X**: Clear all trails
- **Z**: Add random planetesimal
- **Mouse**: Click to add particles

### **Physics Adjustments**
- **1-3**: Adjust coagulation efficiency
- **4-6**: Adjust gas accretion rate
- **7-9**: Adjust disk-planet coupling strength

## Expected Observations

### **Short-term (First few minutes)**
- Dust particles moving in orbital patterns
- Initial coagulation events creating larger aggregates
- Gas disk showing pressure-driven motion
- Temperature gradients visible in particle colors

### **Medium-term (5-15 minutes)**
- Planetesimals forming and growing
- Gravitational interactions creating orbital resonances
- Hill spheres becoming visible around larger bodies
- Gas accretion beginning on massive particles

### **Long-term (15+ minutes)**
- Stable planetary system emerging
- Clear gaps in gas disk around planets
- Spiral density waves propagating through disk
- Circumplanetary disks around large planets
- Accretion streams feeding growing bodies

### **Key Physics Demonstrations**
1. **Orbital mechanics** - Elliptical orbits, resonances
2. **Collision outcomes** - Sticking vs. bouncing vs. fragmentation
3. **Gravitational focusing** - Increased collision rates near massive bodies
4. **Disk-planet interactions** - Gap formation, spiral waves
5. **Composition zoning** - Temperature-dependent material properties
6. **Energy conservation** - Momentum and energy in collisions
7. **Magnetic effects** - Charged particle interactions
8. **Gas dynamics** - Pressure gradients, viscous evolution

## Scientific Accuracy

This simulation incorporates many real physical processes found in actual planet formation:

- **Realistic collision physics** with velocity-dependent outcomes
- **Temperature-dependent composition** matching solar system patterns
- **Hill sphere mechanics** for planetary influence zones
- **Gas disk evolution** with pressure and viscous forces
- **Magnetic field interactions** affecting charged particles
- **Energy conservation** in all physical processes

## Technical Details

- **Built with Processing** (Java-based graphics library)
- **Real-time physics simulation** with 60 FPS target
- **Modular architecture** with separate physics, rendering, and UI systems
- **Scalable parameters** for different simulation scenarios
- **Memory-efficient** particle management

## Educational Value

This simulation is excellent for:
- **Astronomy education** - Planet formation processes
- **Physics learning** - Gravitational mechanics, collision physics
- **Computer science** - Real-time simulation programming
- **Visual learning** - Complex physical processes made visible

## Future Enhancements

Potential improvements could include:
- **3D visualization** for more realistic disk geometry
- **Chemical reaction networks** for complex molecule formation
- **Radiative transfer** for more accurate temperature calculations
- **Multi-threading** for larger particle counts
- **Export capabilities** for data analysis

