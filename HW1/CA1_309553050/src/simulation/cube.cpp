#include "cube.h"
#include <iostream>
#include "Eigen/Dense"
#include <windows.h>
#include "../util/helper.h"
namespace simulation {
constexpr float g_cdK = 2500.0f;
constexpr float g_cdD = 50.0f;

Cube::Cube()
    : particleNumPerEdge(10),
      cubeLength(2.0),
      initialPosition(Eigen::Vector3f(0.0, 0.0, 0.0)),
      springCoefStruct(g_cdK),
      springCoefShear(g_cdK),
      springCoefBending(g_cdK),
      damperCoefStruct(g_cdD),
      damperCoefShear(g_cdD),
      damperCoefBending(g_cdD) {
    particleNumPerFace = particleNumPerEdge * particleNumPerEdge;
    initializeParticle();
    initializeSpring();
}

Cube::Cube(const Eigen::Vector3f &a_kInitPos, const float cubeLength, const int numAtEdge, const float dSpringCoef,
           const float dDamperCoef)
    : particleNumPerEdge(numAtEdge),
      cubeLength(cubeLength),
      initialPosition(a_kInitPos),
      springCoefStruct(dSpringCoef),
      springCoefShear(dSpringCoef),
      springCoefBending(dSpringCoef),
      damperCoefStruct(dDamperCoef),
      damperCoefShear(dDamperCoef),
      damperCoefBending(dDamperCoef) {
    particleNumPerFace = numAtEdge * numAtEdge;
    initializeParticle();
    initializeSpring();
}

int Cube::getParticleNum() const { return static_cast<int>(particles.size()); }

int Cube::getSpringNum() const { return static_cast<int>(springs.size()); }

int Cube::getNumAtEdge() const { return particleNumPerEdge; }

unsigned int Cube::getPointMap(const int a_ciSide, const int a_ciI, const int a_ciJ) {
    int r = -1;

    switch (a_ciSide) {
        case 1:  // [a_ciI][a_ciJ][0] bottom face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ;
            break;
        case 6:  // [a_ciI][a_ciJ][9] top face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ + particleNumPerEdge - 1;
            break;
        case 2:  // [a_ciI][0][a_ciJ] front face
            r = particleNumPerFace * a_ciI + a_ciJ;
            break;
        case 5:  // [a_ciI][9][a_ciJ] back face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * (particleNumPerEdge - 1) + a_ciJ;
            break;
        case 3:  // [0][a_ciI][a_ciJ] left face
            r = particleNumPerEdge * a_ciI + a_ciJ;
            break;
        case 4:  // [9][a_ciI][a_ciJ] ra_ciIght face
            r = particleNumPerFace * (particleNumPerEdge - 1) + particleNumPerEdge * a_ciI + a_ciJ;
            break;
    }

    return r;
}

Particle &Cube::getParticle(int particleIdx) { return particles[particleIdx]; }

std::vector<Particle> *Cube::getParticlePointer() { return &particles; }

Spring &Cube::getSpring(int springIdx) { return springs[springIdx]; }

void Cube::setSpringCoef(const float springCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        springCoefStruct = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        springCoefShear = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        springCoefBending = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::BENDING);
    }
}

void Cube::setDamperCoef(const float damperCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        damperCoefStruct = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        damperCoefShear = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        damperCoefBending = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::BENDING);
    }
}

void Cube::resetCube(const Eigen::Vector3f &offset, const float &rotate) {
    float dTheta = util::radians(rotate);  //  change angle from degree to
                                           //  radian

    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        int i = uiI / particleNumPerFace;
        int j = (uiI / particleNumPerEdge) % particleNumPerEdge;
        int k = uiI % particleNumPerEdge;
        float offset_x = (float)((i - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
        float offset_y = (float)((j - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
        float offset_z = (float)((k - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));

        Eigen::Vector3f RotateVec(offset_x, offset_y,
                                  offset_z);  //  vector from center of cube to the particle

        Eigen::AngleAxis<float> rotation(dTheta, Eigen::Vector3f(1.0f, 0.0f, 1.0f).normalized());

        RotateVec = rotation * RotateVec;

        particles[uiI].setPosition(initialPosition + offset + RotateVec);
        particles[uiI].setForce(Eigen::Vector3f::Zero());
        particles[uiI].setVelocity(Eigen::Vector3f::Zero());
    }
}

void Cube::addForceField(const Eigen::Vector3f &force) {
    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        particles[uiI].setAcceleration(force);
    }
}


void Cube::computeInternalForce() {
    // TODO

    Eigen::Vector3f SF;  // computeSpringForce
    Eigen::Vector3f DF;  // computeDamperForce

        for (int i = 0; i < springs.size(); i++) {
            
            
            SF = computeSpringForce(
                particles[springs[i].getSpringStartID()].getPosition(),
                particles[springs[i].getSpringEndID()].getPosition(),
                springs[i].getSpringCoef(),
                springs[i].getSpringRestLength()
            );
            
            DF = computeDamperForce(
                particles[springs[i].getSpringStartID()].getPosition(),
                particles[springs[i].getSpringEndID()].getPosition(),
                                    particles[springs[i].getSpringStartID()].getForce(),
                particles[springs[i].getSpringEndID()].getForce(),  
                springs[i].getDamperCoef()
            );
            

            particles[springs[i].getSpringStartID()].addForce(SF+DF);           //彈簧力+阻尼力
            particles[springs[i].getSpringEndID()].addForce(-1* (SF+DF));  //彈簧力+阻尼力的反作用力
            
        }
   
}

Eigen::Vector3f Cube::computeSpringForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const float springCoef, const float restLength) {
    // TODO
    Eigen::Vector3f SF;
    Eigen::Vector3f positionBA = positionA - positionB;
    float distance =
    sqrt(positionBA[0] * positionBA[0] + positionBA[1] * positionBA[1] + positionBA[2] * positionBA[2]);
    SF = -1 * springCoef * (distance - restLength) * positionBA / distance;
    return SF;
}

Eigen::Vector3f Cube::computeDamperForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const Eigen::Vector3f &velocityA, const Eigen::Vector3f &velocityB,
                                         const float damperCoef) {
    // TODO
    Eigen::Vector3f DF;
    Eigen::Vector3f positionBA = positionA - positionB;
    Eigen::Vector3f velocityBA = velocityA - velocityB;
    float distance =
        sqrt(positionBA[0] * positionBA[0] + positionBA[1] * positionBA[1] + positionBA[2] * positionBA[2]);
    float v_dot_x = velocityBA[0] * positionBA[0] + velocityBA[1] * positionBA[1] + velocityBA[2] * positionBA[2];
    DF = -1 * damperCoef * v_dot_x * positionBA / (distance*distance);
    return DF;
}

void Cube::initializeParticle() {
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                Particle Particle;
                float offset_x = (float)((i - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                float offset_y = (float)((j - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                float offset_z = (float)((k - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                Particle.setPosition(Eigen::Vector3f(initialPosition(0) + offset_x, initialPosition(1) + offset_y,
                                                     initialPosition(2) + offset_z));
                //std::cout<<Particle.getPosition();
                particles.push_back(Particle);
                
            }
        }
    }
}
#define particleNumPerEdge 10
#define STRUCT_SPRING_NUM  2700 // 9 * 10 * 20 + 9 * 10 * 10
#define BENDING_SPRING_NUM 2400 //  8 * 10 * 20 + 10 * 10 * 8 
#define SHEAR_SPRING_NUM 12150

void Cube::initializeSpring() {
    int iParticleID = 0; 
    int iNeighborID = 0;
    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
    float struct_spring_length = sqrt(((particles[1].getPosition() - particles[0].getPosition())[0]) *((particles[1].getPosition() - particles[0].getPosition())[0]) +
                                     ((particles[1].getPosition() - particles[0].getPosition())[1]) *((particles[1].getPosition() - particles[0].getPosition())[1]) +
                                 ((particles[1].getPosition() - particles[0].getPosition())[2]) *((particles[1].getPosition() - particles[0].getPosition())[2])) ;

    float bend_spring_length =  struct_spring_length * 2;

    float shear_spring_length =  sqrt(struct_spring_length * struct_spring_length + struct_spring_length * struct_spring_length);//平面對角長度
    float shear_spring_length3 =
        sqrt(struct_spring_length * struct_spring_length + struct_spring_length * struct_spring_length +  struct_spring_length * struct_spring_length); //立體對角長度

    enum {
        X_id = 1,
        Y_id,  //  = 2
        Z_id   //  = 3
    };
    // TODO FINISH
    particleNumPerFace = particleNumPerEdge*particleNumPerEdge;
    float springCoef = 0.05;
    float damperCoef = 0.05;
    Spring *struct_spring = new Spring[STRUCT_SPRING_NUM];  
    
    //Z axis
    int counter = 0;
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge - 1; k++){
                
                //STRUCT spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + j * particleNumPerEdge + (k + 1);
                struct_spring[counter] = Spring(iParticleID, iNeighborID, struct_spring_length, springCoef, damperCoef,
                                                simulation::Spring::SpringType::STRUCT);
                springs.push_back(struct_spring[counter]); 
                counter++;
            }
        }
    }
    // Y axis
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge-1; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                // STRUCT spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + (j+1) * particleNumPerEdge + k;
                struct_spring[counter] = Spring(iParticleID, iNeighborID, struct_spring_length, springCoef, damperCoef,
                                                simulation::Spring::SpringType::STRUCT);
                springs.push_back(struct_spring[counter]);
                counter++;
            }
        }
    }
    // X axis
    for (int i = 0; i < particleNumPerEdge-1; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                // STRUCT spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i+1) * particleNumPerFace + j * particleNumPerEdge + k;
                struct_spring[counter] = Spring(iParticleID, iNeighborID, struct_spring_length, springCoef, damperCoef,
                                                simulation::Spring::SpringType::STRUCT);
                springs.push_back(struct_spring[counter]);
                counter++;
            }
        }
    }
    //std::cout << counter << std::endl;
    Spring *bend_spring = new Spring[BENDING_SPRING_NUM];  
    
    // Z axis
    counter = 0;
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge - 2; k++) {
                // BENDING spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + j * particleNumPerEdge + (k + 2);
                bend_spring[counter] = Spring(iParticleID, iNeighborID, bend_spring_length, springCoef, damperCoef,
                                                simulation::Spring::SpringType::BENDING);
                springs.push_back(bend_spring[counter]);
                counter++;
            }
        }
    }
    // Y axis
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge - 2; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                // BENDING spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + (j + 2) * particleNumPerEdge + k;
                bend_spring[counter] = Spring(iParticleID, iNeighborID, bend_spring_length, springCoef, damperCoef,
                                                simulation::Spring::SpringType::BENDING);
                springs.push_back(bend_spring[counter]);
                counter++;
            }
        }
    }
    // X axis
    for (int i = 0; i < particleNumPerEdge - 2; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                // BENDING spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i + 2) * particleNumPerFace + j * particleNumPerEdge + k;
                bend_spring[counter] = Spring(iParticleID, iNeighborID, bend_spring_length, springCoef, damperCoef,
                                                simulation::Spring::SpringType::BENDING);
                springs.push_back(bend_spring[counter]);
                counter++;
            }
        }
    }
    
    //std::cout << counter << std::endl;
    Spring *shear_spring = new Spring[SHEAR_SPRING_NUM];  
    counter = 0;

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge-1; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + (j+1) * particleNumPerEdge + (k + 1);
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 0; i < particleNumPerEdge-1; i++) {
        for (int j = 0; j < particleNumPerEdge - 1; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i+1) * particleNumPerFace + (j + 1) * particleNumPerEdge + k;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 0; i < particleNumPerEdge - 1; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i + 1) * particleNumPerFace + j * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 0; i < particleNumPerEdge - 1; i++) {
        for (int j = 0; j < particleNumPerEdge-1; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i + 1) * particleNumPerFace + (j+1) * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length3, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    
    
    for (int i = 1; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge-1; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i-1) * particleNumPerFace + (j + 1) * particleNumPerEdge + k;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    
    for (int i = 1; i < particleNumPerEdge ; i++) {
        for (int j = 0; j < particleNumPerEdge-1 ; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i  * particleNumPerFace + (j+1) * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    
    for (int i = 1; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i-1) * particleNumPerFace + j * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    
    for (int i = 1; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge-1; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i-1) * particleNumPerFace + (j+1)  * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length3, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + (j-1) * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 0; i < particleNumPerEdge-1; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i+1)  * particleNumPerFace +j  * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 0; i < particleNumPerEdge-1; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i+1) * particleNumPerFace + (j-1) * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length3, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }

    for (int i = 1; i < particleNumPerEdge; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = i * particleNumPerFace + (j-1)  * particleNumPerEdge + (k+1);
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 1; i < particleNumPerEdge; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i-1) * particleNumPerFace + j * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 1; i < particleNumPerEdge; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge-1; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i-1) * particleNumPerFace + (j-1) * particleNumPerEdge + k+1;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length3, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge-1; i++) {
        for (int j = 0; j < particleNumPerEdge-1; j++) {
            for (int k = 1; k < particleNumPerEdge; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i + 1) * particleNumPerFace + (j + 1) * particleNumPerEdge + k;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    for (int i = 0; i < particleNumPerEdge-1; i++) {
        for (int j = 1; j < particleNumPerEdge; j++) {
            for (int k = 1; k < particleNumPerEdge; k++) {
                // SHEAR spring
                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                iNeighborID = (i + 1) * particleNumPerFace + (j - 1) * particleNumPerEdge + k;
                shear_spring[counter] = Spring(iParticleID, iNeighborID, shear_spring_length, springCoef, damperCoef,
                                               simulation::Spring::SpringType::SHEAR);
                springs.push_back(shear_spring[counter]);
                counter++;
            }
        }
    }
    //std::cout << counter << std::endl;
}

void Cube::updateSpringCoef(const float a_cdSpringCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setSpringCoef(a_cdSpringCoef);
        }
    }
}

void Cube::updateDamperCoef(const float a_cdDamperCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setDamperCoef(a_cdDamperCoef);
        }
    }
}
}  //  namespace simulation
