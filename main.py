"""Just for testing the computation burden of the blade_elements"""
from Bladelement import BladeAeroCalculator as BAC
import numpy as np
import random
import time

Test_bac = BAC('SimpleFlapper0Wing40BLE_X_pos.npy',\
                'SimpleFlapper0Wing40BLE_Y_pos.npy',bladeElementNo=40)
Test_bac.SetVirturalWingPlaneRelative2Wing([1,0,0,0,1,0,0,0,1])

indicator_1  = 0
Compare_loop_length = 4000


start_time = time.time()
while indicator_1 < Compare_loop_length:
    LU_OrVec = [1,0,0,0,1,0,0,0,1]
    vel_FreeFlow = np.array([0,0,0])
    LU_velVec = np.random.randn(6)

    LU_wing_Rotation_matrix = \
        np.matrix([[LU_OrVec[0], LU_OrVec[1], LU_OrVec[2]],\
                [LU_OrVec[3], LU_OrVec[4], LU_OrVec[5]],\
                [LU_OrVec[6], LU_OrVec[7], LU_OrVec[8]]])

    d_Real_motor_LU_wing_joint_sensor_value = random.uniform(1,30)

    Test_bac.RequestWingOrientation2InertiaFrame (LU_OrVec)
    Test_bac.RequestWingPlaneDirection()
    Test_bac.RequestVelocities(vel_FreeFlow, LU_velVec[0:3],\
                                            LU_velVec[3:6],d_Real_motor_LU_wing_joint_sensor_value)
    Test_bac.CalcEffectiveVelocity()
    Test_bac.CalcAoA()
    Test_bac.CopmputeAerodynamicForce()

    LU_r_shift_in_wing = np.array([Test_bac. X_pos_r, Test_bac. Y_pos_r, 0])
    LU_r_shift = np.array(np.matmul(LU_wing_Rotation_matrix, LU_r_shift_in_wing)).squeeze().tolist()
    LU_r = np.array(np.sum(Test_bac. F_r, axis=0)).squeeze()

    LU_a_shift_in_wing = np.array([Test_bac. X_pos_a, Test_bac. Y_pos_a, 0])
    LU_a_shift = np.array(np.matmul(LU_wing_Rotation_matrix, LU_a_shift_in_wing)).squeeze().tolist()
    LU_a = np.array(np.sum(Test_bac. F_a, axis=0)).squeeze()

    LU_t_shift_in_wing = np.array([Test_bac. X_pos_t, Test_bac. Y_pos_t, 0])
    LU_t_shift = np.array(np.matmul(LU_wing_Rotation_matrix, LU_t_shift_in_wing)).squeeze().tolist()
    LU_drag = np.array(np.sum(Test_bac. F_t_drag, axis=0)).squeeze()

    LU_lift = np.array(np.sum(Test_bac. F_t_lift, axis=0)).squeeze()
    indicator_1 = indicator_1 + 1

end_time = time.time()

print('A single loop is complete!!')
print('Comsume time of ', end_time - start_time)