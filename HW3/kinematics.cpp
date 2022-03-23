#include "simulation/kinematics.h"
#include <math.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include "Eigen/Dense"
#include "acclaim/bone.h"
#include "util/helper.h"
using namespace std;
namespace kinematics {
Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    Eigen::MatrixXd Jacobian_plus;
    Eigen::Matrix3Xd J = Jacobian.topRows(3);
    Eigen::Vector3d V = target.head<3>();
    Eigen::VectorXd result;
    Jacobian_plus = J.transpose() * (J * J.transpose()).inverse();  // J+ = (JT*J)^-1*JT =>>> JT * (J * JT) ^ -1
    result = (Jacobian_plus * V).transpose();
    return result;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
                             acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;

    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;
    acclaim::Bone* current_bone = end_bone;
    size_t bone_num = 0;
  
    float D = current_bone->length;
    float d = current_bone->length;
    while (current_bone->parent->idx != start_bone->idx && current_bone->parent->idx != 0) {
        if (current_bone->parent->idx != 0) {
            current_bone = current_bone->parent;
            D += current_bone->length;
            d -= current_bone->length;
            bone_num++;
        }
    }

    bone_num = bone_num + 2;  // include both side  = 9
    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);  // (4,27)
    Eigen::Vector4d r;
    Eigen::Vector3d J_temp;
    Eigen::Vector3d ax;
    Eigen::Vector3d ay;
    Eigen::Vector3d az;
    Eigen::Vector3d ea;
    

    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        Eigen::Vector3d r3, J_num;
        Jacobian.setZero();
        if (desiredVector.norm() < epsilon) {
            break;
        }

        current_bone = end_bone;
        for (int i = 0; i < bone_num; i++) {
            r = target_pos - (current_bone->end_position);
            r3 = r.head<3>();  // 3ºûªºr
            ax = current_bone->rotation * Eigen::Vector4d(1, 0, 0, 0).normalized().head<3>();
            ay = current_bone->rotation * Eigen::Vector4d(0, 1, 0, 0).normalized().head<3>();
            az = current_bone->rotation * Eigen::Vector4d(0, 0, 1, 0).normalized().head<3>();
            if (current_bone->dofrx == 1) {
                Jacobian(0, 3*i) = ax.cross(r3)[0];
                Jacobian(1, 3*i) = ax.cross(r3)[1];
                Jacobian(2, 3*i) = ax.cross(r3)[2];
            } 
            if (current_bone->dofry == 1) {
                Jacobian(0, 3 * i + 1) = ay.cross(r3)[0];
                Jacobian(1, 3 * i + 1) = ay.cross(r3)[1];
                Jacobian(2, 3 * i + 1) = ay.cross(r3)[2];
            }
            if (current_bone->dofrz == 1) {
                Jacobian(0, 3 * i + 2) = az.cross(r3)[0];
                Jacobian(1, 3 * i + 2) = az.cross(r3)[1];
                Jacobian(2, 3 * i + 2) = az.cross(r3)[2];
            } 
            current_bone = current_bone->parent;
        }
        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        current_bone = end_bone;
        int counter = 0;
        for (int i = 0; i < bone_num; i++) {
            posture.bone_rotations[current_bone->idx][0] =
                posture.bone_rotations[current_bone->idx][0] + deltatheta[counter];
            counter++;
            posture.bone_rotations[current_bone->idx][1] =
                posture.bone_rotations[current_bone->idx][1] + deltatheta[counter];
            counter++;
            posture.bone_rotations[current_bone->idx][2] =
                posture.bone_rotations[current_bone->idx][2] + deltatheta[counter];
            counter++;
            current_bone = current_bone->parent;
        }
    }

    if ((target_pos - start_bone->start_position).norm() > d && (target_pos - start_bone->start_position).norm() < D) {
        return true;
    }
    else {
        return false;
    }
}
}  // namespace kinematics
