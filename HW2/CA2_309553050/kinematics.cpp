#include "simulation/kinematics.h"
#include <iostream>
#include "Eigen/Dense"
#include <cmath>
#include "acclaim/bone.h"
#include "util/helper.h"
using namespace std;

namespace kinematics {

void connectbone(const acclaim::Posture& posture, acclaim::Bone* bone, int i) {

        std::vector<Eigen::Vector4d> rotation = posture.bone_rotations;  
        acclaim::Bone* root = bone;
        std::vector<Eigen::Vector4d> V(31);
        Eigen::Vector3d temp;
        Eigen::Matrix<double, 3, 3> rotationMatrix;
        bone->start_position = bone->parent->end_position;
        rotationMatrix = bone->parent->rotation.linear()* bone->rot_parent_current.linear() *
                         util::rotateDegreeZYX(rotation[i]).toRotationMatrix(); 

        temp = rotationMatrix * (Eigen::Vector3d(bone->dir[0], bone->dir[1], bone->dir[2]) * bone->length);
        V[i] = Eigen::Vector4d(temp[0], temp[1], temp[2], 0);  // V2*i2
        bone->end_position = bone->start_position + V[i];
        bone->rotation.linear() = rotationMatrix;

}

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO
    
    std::vector<Eigen::Vector4d> rotation = posture.bone_rotations;        // 31 item 很多都有值
    std::vector<Eigen::Vector4d> translation = posture.bone_translations;  // 31 item 只有第一筆有值
    acclaim::Bone* root = bone;
    std::vector<Eigen::Vector4d> V(31);
    Eigen::Vector3d temp;
    Eigen::Matrix<double, 3, 3> rotationMatrix;
    int i;
    
    //root 初始值設定 0
    root->start_position = translation[0];
  
    rotationMatrix = util::rotateDegreeZYX(rotation[0]).toRotationMatrix() ;  // Rasf * Ramc
    temp = rotationMatrix * (Eigen::Vector3d(bone->dir[0], bone->dir[1], bone->dir[2]) * bone->length);
    V[0] = Eigen::Vector4d(temp[0], temp[1], temp[2], 0);
    root->end_position = translation[0] + V[0]; 
    bone->rotation.linear() = rotationMatrix;
   
    
    bone = root;
    for (int i = 1; i < 6; i++) {  // lhipjoint > ltoes 1-5
        bone = bone->child;
        connectbone(posture, bone, bone->idx);
    }
    bone = root->child->sibling;  // rlipjoint > rtoes 6-10
    for (int i = 6; i < 10; i++) {
        connectbone(posture, bone, bone->idx);
        bone = bone->child;
    }
    connectbone(posture, bone, bone->idx);
    bone = root->child->sibling->sibling;  // lowerback > head 11-16
    for (int i = 11; i < 16; i++) {
        connectbone(posture, bone, bone->idx);
        bone = bone->child;
    }
    connectbone(posture, bone, bone->idx);
    bone = root->child->sibling->sibling->child->child->child->sibling;  // lclavicle > lfingers 17-22
    for (int i = 17; i < 22; i++) {
        connectbone(posture, bone, bone->idx);
        bone = bone->child;
    }
    connectbone(posture, bone, bone->idx);
    bone =
        root->child->sibling->sibling->child->child->child->sibling->child->child->child->child->sibling;  // lthumb  23
    connectbone(posture, bone, bone->idx);
    bone = root->child->sibling->sibling->child->child->child->sibling->sibling;  // rclavicle > rfingers 24-29
    for (int i = 24; i < 29; i++) {

        connectbone(posture, bone, bone->idx);
        bone = bone->child;
    }
    connectbone(posture, bone, bone->idx);
    bone = root->child->sibling->sibling->child->child->child->sibling->sibling->child->child->child->child
               ->sibling;  // rthumb //30
    connectbone(posture, bone, bone->idx);
   
}

std::vector<acclaim::Posture> timeWarper(const std::vector<acclaim::Posture>& postures, int keyframe_old,
                                         int keyframe_new) {
    int total_frames = static_cast<int>(postures.size());
    int total_bones = static_cast<int>(postures[0].bone_rotations.size());
    float frame_1, frame_2;
    std::vector<acclaim::Posture> new_postures = postures;
    for (int i = 0; i < total_frames; ++i) {
        for (int j = 0; j < total_bones; ++j) {
            // TODO
            //將150重新切成160 (0.94:0.06)

            if (i > 0 && i < 150) {
                frame_1 = (i % 16);
                frame_1 = frame_1 / 16;
                frame_2 = (16 - i % 16);
                frame_2 = frame_2 / 16;
                new_postures[i].bone_translations[j] = postures[i - 1 + i/16].bone_translations[j] * frame_1 + postures[i+i/16].bone_translations[j] * frame_2;
                new_postures[i].bone_rotations[j] = postures[i - 1 + i / 16].bone_rotations[j] * frame_1 +
                                                    postures[i + i / 16].bone_rotations[j] * frame_2;


            } else if (i == 0) {
                new_postures[i].bone_translations[j] = postures[0].bone_translations[j];
                new_postures[i].bone_rotations[j] = postures[0].bone_rotations[j];
            } else if (i == 150){
                new_postures[i].bone_translations[j] = postures[160].bone_translations[j];
                new_postures[i].bone_rotations[j] = postures[160].bone_rotations[j];
            } else if (i > 150 && i <450) { //i>150
                frame_1 = (i % 10);
                frame_1 = frame_1 / 10;
                frame_2 = (10 - i % 10);
                frame_2 = frame_2 / 10;
                new_postures[i].bone_translations[j] = postures[i+10 - (i-150)/10].bone_translations[j] * frame_1 +
                                                       postures[i+11 -(i-150)/10].bone_translations[j] * frame_2;

                new_postures[i].bone_rotations[j] = postures[i + 10 - (i - 150) / 10].bone_rotations[j] * frame_1 +
                                                    postures[i + 11 - (i - 150) / 10].bone_rotations[j] * frame_2;
            } 
        }
    }
    return new_postures;
}
}  // namespace kinematics
