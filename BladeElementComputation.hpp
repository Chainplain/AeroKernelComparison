#ifndef _BladeElementComputation_H_
#define _BladeElementComputation_H_
#include"Eigen/Core"
#include"Eigen/Dense"
#include<vector>
#include<cmath>
#include<string>
#define data_cat double

// #include"BladeElementComputation.cpp"

using namespace std;

data_cat SMALL_CONSTANT_4_NORMALIZATION = 0.000000001;

/*This class is refractured from our previous python version,
https://github.com/Chainplain/Flapping_wing_Simu
the algorithm is completely the same. We just would like to check 
and compare the runtime of them.*/
class BladeAeroCalculator
{
    private:
        /*This value is simply determined by the input data.*/
        int Number_of_BladeElements = 50;
        Eigen::VectorX<data_cat>  x_data;
        Eigen::VectorX<data_cat>  y_data;

    public:
        Eigen::Matrix<data_cat, 3, 3> wing_rotation_relative_to_inertia;
        Eigen::Matrix<data_cat, 3, 3> virtural_Wing_Plane_relative_to_Wing;
        
        
        data_cat Filter_Constant = 0.2;

        //Please note that v_AoA is not the rate of AoA, but the wing plane.
        data_cat v_AoA = 0.0;
        data_cat BladeElementswidth = 0.0;

        Eigen::Vector<data_cat, 3> v_infi;
        Eigen::Vector<data_cat, 3> v_t;

        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> v_r;

        Eigen::MatrixX<data_cat> rotation_eFV;
        Eigen::MatrixX<data_cat> added_mass_FV;
        Eigen::MatrixX<data_cat> last_added_mass_FV;

        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> eFV;

        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> d_eFV;

        /*(Number_of_BladeElements X 1) - dim */
        Eigen::VectorX<data_cat> AoA;

        // data_cat v_AoA;

        Eigen::VectorX<data_cat> Zeros_vect_with_dim_of_Number_of_BladeElements;
        Eigen::VectorX<data_cat> Ones_vect_with_dim_of_Number_of_BladeElements;

        Eigen::MatrixX<data_cat> lasteFV;
        Eigen::MatrixX<data_cat> lastAoA;

        data_cat air_Density = 1.29;
        data_cat x_hat_0 = 0;
        Eigen::MatrixX<data_cat> posLeadingEdgeVec;
        //simulation time between two step, with the unit of second.
        data_cat length_period = 0.001;
        

        Eigen::Vector<data_cat, 3> CoM;
        data_cat avReynoldsNumber = 7000;

        //Rename from LERT, the new name is more explicitive.
        Eigen::Vector<data_cat, 3> leading_edge_relative_to_inertia;

        //Rename from chordD
        Eigen::Vector<data_cat, 3> chordwise_direction_relative_to_inertia;

        Eigen::Vector<data_cat, 3> leading_edge_relative_to_wing;

        Eigen::Vector<data_cat, 3> chordwise_direction_relative_to_wing;

        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> F_t_lift;
        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> F_t_drag;
        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> F_r;
        /*(Number_of_BladeElements X 3) - dim */
        Eigen::MatrixX<data_cat> F_a;

        data_cat X_pos_t =0.0;
        data_cat Y_pos_t =0.0;

        data_cat X_pos_r =0.0;
        data_cat Y_pos_r =0.0;

        data_cat X_pos_a =0.0;
        data_cat Y_pos_a =0.0;

        data_cat A_L = 1.966 - 3.94 *  pow(avReynoldsNumber, -0.429);
        data_cat A_D = 1.873 - 3.14 *  pow(avReynoldsNumber, -0.369);
        data_cat C_D_0 = 0.031 + 10.48 * pow(avReynoldsNumber, -0.764);
        

        BladeAeroCalculator(string x_filename, string y_filename);

        /*Vector  normilzer bt sqrt*/
        static Eigen::VectorX<data_cat>  vec_Normalize(Eigen::VectorX<data_cat> & input_vec);
        static Eigen::MatrixX<data_cat>  mat_Normalize(Eigen::MatrixX<data_cat> & input_mat);
         /*If the input is a matrix, then we render this as vertically stacked vector,
        and we return the */
        static Eigen::MatrixX<data_cat>  mat_SecNorm(Eigen::MatrixX<data_cat> & input_mat);
        static Eigen::MatrixX<data_cat>  mat_vec_cross(Eigen::MatrixX<data_cat> & input_mat, Eigen::VectorX<data_cat> & input_vec);
        static Eigen::MatrixX<data_cat>  vec_mat_cross(Eigen::VectorX<data_cat> & input_vec, Eigen::MatrixX<data_cat> & input_mat);

       
        int SetAverageReynoldsNumberAndDecideCoeffs(data_cat Input_Desire_avReynoldsNumber);
        int SetVirturalWingPlaneRelative2Wing(Eigen::Matrix<data_cat, 3, 3> \
            Input_wing_rotation_relative_to_inertia);
        int RequestWingOrientation2InertiaFrame(Eigen::Matrix<data_cat, 3, 3> \
            Input_virtural_Wing_Plane_Relative_to_Wing);
        int setlength_period( data_cat l_p);
        /*In the python version: 
        https://github.com/Chainplain/Flapping_wing_Simu
        this is function called RequestWingPlaneDirection,
        which is not named properly.*/
        int CalcWingPlaneDirection();
        int RequestVelocities(Eigen::Vector<data_cat, 3> vel_FreeFlow,\
                              Eigen::Vector<data_cat, 3> vel_body_translation,\
                              Eigen::VectorX<data_cat> vel_wing_rotation,\
                              data_cat vel_AoA);
        int CalcEffectiveVelocity();
        int CalcAoA();
        int CopmputeAerodynamicForce(); 
};

#endif