#include "terrain.h"
#include "iostream"
#include <stdexcept>

#include "../util/helper.h"

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Sphere:
            return std::make_unique<SphereTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        case simulation::TerrainType::TiltedPlane:
            return std::make_unique<TiltedPlaneTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }





void PlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f; //kr
    constexpr float coefFriction = 0.3f; //kf
    // TODO 
    int par_num = cube.getParticleNum();

    Eigen::Vector3f coll;
    Eigen::Vector3f fc;
    coll[0] = 0;
    coll[1] = 4.9f;
    coll[2] = 0;
    for (int i = 0; i < 1000; i += 1) {
        if (cube.getParticle(i).getPosition()[1] < position[1] + eEPSILON) {  //¸I¼²§P©w
            cube.getParticle(i).addForce(coll);
        }
    }
}

// SphereTerrain //

SphereTerrain::SphereTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType SphereTerrain::getType() { return TerrainType::Sphere; }

void SphereTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    


}

// BowlTerrain //

BowlTerrain::BowlTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO
}

// TiltedPlaneTerrain //

TiltedPlaneTerrain::TiltedPlaneTerrain() { modelMatrix = util::rotateDegree(0, 0, -45) * util::scale(60, 1, 60); }

TerrainType TiltedPlaneTerrain::getType() { return TerrainType::TiltedPlane; }

void TiltedPlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;

    float radius = 3.0f;
    float mass = 10.0f;
    // TODO
    int par_num = cube.getParticleNum();

    Eigen::Vector3f coll;
    coll[0] = 4.9f * sqrt(2);
    coll[1] = 4.9f * sqrt(2);
    coll[2] = 0;

    for (int i = 0; i < 1000; i += 1) {
        if (cube.getParticle(i).getPosition()[1] + cube.getParticle(i).getPosition()[0] <
            position[1]  + eEPSILON + 1) {  //¸I¼²§P©w
            cube.getParticle(i).addForce(coll);
        }
    }


}
}  // namespace simulation
