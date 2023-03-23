#include "BladeElementComputation.hpp"

/*Thanks to bro llohse: https://github.com/llohse/libnpy
we can load the npy documents */
#include "npy.hpp"
// #include <cassert>
using namespace std;
#define PI 3.1415926

BladeAeroCalculator::BladeAeroCalculator(string x_filename, string y_filename){

    vector<unsigned long> shape {};
    bool fortran_order;
    vector<double> tempx, tempy;
    npy::LoadArrayFromNumpy(x_filename, shape, fortran_order, tempx);
    npy::LoadArrayFromNumpy(y_filename, shape, fortran_order, tempy);

    vector<data_cat> tempx_new, tempy_new;
    for (int i = 0; i < tempx.size(); i++){
        tempx_new.push_back( data_cat(tempx.at(i)));
    }
    for (int i = 0; i < tempy.size(); i++){
        tempy_new.push_back( data_cat(tempy.at(i)));
    }
    this->x_data = Eigen::Map<Eigen::VectorX<data_cat>, \
        Eigen::Unaligned>(tempx_new.data(), tempx_new.size());
    this->y_data = Eigen::Map<Eigen::VectorX<data_cat>, \
        Eigen::Unaligned>(tempy_new.data(), tempy_new.size());
    this->Number_of_BladeElements = tempx.size();

    wing_rotation_relative_to_inertia = Eigen::Matrix<data_cat, 3, 3> ::Identity(3,3);
    virtural_Wing_Plane_relative_to_Wing  = Eigen::Matrix<data_cat, 3, 3> ::Identity(3,3);
    this->BladeElementswidth = (this->x_data.maxCoeff() - \
        this->x_data.minCoeff()) / (this->Number_of_BladeElements);
    
    this->v_infi << 0.0, 0.0, 0.0;
    this->v_t << 0.0, 0.0, 0.0;
    this->v_r = Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements, 3);
    this->eFV = Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements, 3);
    this->d_eFV = Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements, 3);

    this->added_mass_FV = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements, 3);
    this->AoA = Eigen::VectorX<data_cat>::Zero(this->Number_of_BladeElements);
    this->v_AoA = 0.0;

    this->posLeadingEdgeVec = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements, 3);
    this->Zeros_vect_with_dim_of_Number_of_BladeElements = \
        Eigen::VectorX<data_cat>::Zero(this->Number_of_BladeElements);
    this->Ones_vect_with_dim_of_Number_of_BladeElements = \
        Eigen::VectorX<data_cat>::Ones(this->Number_of_BladeElements);
    this->posLeadingEdgeVec<< this->x_data, this->Zeros_vect_with_dim_of_Number_of_BladeElements,\
        this->Zeros_vect_with_dim_of_Number_of_BladeElements;

    this->F_t_lift = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements,3);
    this->F_t_drag = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements,3);
    this->F_r = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements,3);
    this->F_a = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements,3);
    
    this->CoM<<0.0, 0.0, 0.0;
    Eigen::VectorX<data_cat> COM_x_with_dim_of_Number_of_BladeElements = this->CoM[0] *\
        this->Ones_vect_with_dim_of_Number_of_BladeElements;
    Eigen::VectorX<data_cat> COM_y_with_dim_of_Number_of_BladeElements = this->CoM[1] *\
        this->Ones_vect_with_dim_of_Number_of_BladeElements;
    Eigen::VectorX<data_cat> COM_z_with_dim_of_Number_of_BladeElements = this->CoM[2] *\
        this->Ones_vect_with_dim_of_Number_of_BladeElements;
    Eigen::MatrixX<data_cat> COM_with_dim_of_Number_of_BladeElements = \
        Eigen::MatrixX<data_cat>::Zero(this->Number_of_BladeElements, 3);
    COM_with_dim_of_Number_of_BladeElements<<\
        COM_x_with_dim_of_Number_of_BladeElements,\
        COM_y_with_dim_of_Number_of_BladeElements,\
        COM_z_with_dim_of_Number_of_BladeElements;
    
    this->posLeadingEdgeVec = this->posLeadingEdgeVec - COM_with_dim_of_Number_of_BladeElements;

    this->leading_edge_relative_to_wing << 1,0,0;
    this->chordwise_direction_relative_to_wing << 0,1,0; 

    cout<<"posLeadingEdgeVec:"<<posLeadingEdgeVec<<endl;
    cout<<"@STATE:Initial values settled!!!"<<endl;
}

Eigen::VectorX<data_cat>  BladeAeroCalculator:: vec_Normalize(Eigen::VectorX<data_cat> input){
    data_cat frac = SMALL_CONSTANT_4_NORMALIZATION + input.norm();
    return 1 / frac * input;
}

Eigen::MatrixX<data_cat> BladeAeroCalculator:: mat_Normalize(Eigen::MatrixX<data_cat> input_mat){
    Eigen::MatrixX<data_cat> output_mat = Eigen::MatrixX<data_cat>::Zero(input_mat.rows(), input_mat.cols());
    for (int i = 0; i < input_mat.rows(); i++){
        output_mat.row(i) = BladeAeroCalculator:: vec_Normalize((Eigen::VectorX<data_cat>)input_mat.row(i));
    }
    return output_mat;
}

Eigen::MatrixX<data_cat>  BladeAeroCalculator:: mat_SecNorm(Eigen::MatrixX<data_cat> input_mat){
    Eigen::MatrixX<data_cat> output(input_mat.rows(), input_mat.cols());
    output = input_mat.cwiseProduct(input_mat);
    Eigen::MatrixX<data_cat> outputsum(input_mat.rows(), 1);
    outputsum = output.rowwise().sum();
    outputsum = outputsum.cwiseSqrt();
    return outputsum;
}

Eigen::MatrixX<data_cat> BladeAeroCalculator:: mat_vec_cross(Eigen::MatrixX<data_cat> input_mat, \
    Eigen::VectorX<data_cat> input_vec){
  
    if ((input_vec.size() == 3) && (input_mat.cols() == 3)){
        Eigen::MatrixX<data_cat> output = Eigen::MatrixX<data_cat>::Zero(input_mat.rows(), input_mat.cols());
        Eigen::Vector<data_cat,3> temp;
        for (int i = 0; i < input_mat.rows(); i++){
            temp = ((Eigen::Vector<data_cat,3>) input_mat.row(i)).cross(\
                (Eigen::Vector<data_cat,3>)input_vec);
            output.row(i)=temp;
        }
        return  output;              
    }
    
     Eigen::MatrixX<data_cat> output(1,1);
     cout<<"@WARNNING:Please notice this message: something aboud cross product went wrong."<<endl;
     output<<0;
     return output;
}

Eigen::MatrixX<data_cat> BladeAeroCalculator:: vec_mat_cross(Eigen::VectorX<data_cat> input_vec, \
    Eigen::MatrixX<data_cat> input_mat){
        return - BladeAeroCalculator::mat_vec_cross(input_mat, input_vec);
}

int BladeAeroCalculator::SetAverageReynoldsNumberAndDecideCoeffs(data_cat Input_Desire_avReynoldsNumber){
    this->avReynoldsNumber = Input_Desire_avReynoldsNumber;
    this->A_L = 1.966 - 3.94 *  pow(avReynoldsNumber, -0.429);
    this->A_D = 1.873 - 3.14 *  pow(avReynoldsNumber, -0.369);
    this->C_D_0 = 0.031 + 10.48 * pow(avReynoldsNumber, -0.764);
    cout<<"@CHANGE: The average Reynolds Number(avRey) is set to" \
        << this->avReynoldsNumber << '.'<<endl;
    cout<<"@CHANGE: The aerodynamics coefficients A_L, A_D, C_D_0 are accordingly reset."<<endl;
    return 1;
}

int BladeAeroCalculator::SetVirturalWingPlaneRelative2Wing(Eigen::Matrix<data_cat, 3, 3> \
    Input_wing_rotation_relative_to_inertia){
    this->wing_rotation_relative_to_inertia = Input_wing_rotation_relative_to_inertia;
    return 1;
}

/*Frequently get the orientation，such that it is divided into request.*/
int BladeAeroCalculator::RequestWingOrientation2InertiaFrame(Eigen::Matrix<data_cat, 3, 3> \
    Input_virtural_Wing_Plane_Relative_to_Wing){
    this->virtural_Wing_Plane_relative_to_Wing = Input_virtural_Wing_Plane_Relative_to_Wing;
    return 1;
}

int BladeAeroCalculator::CalcWingPlaneDirection(){
    this->leading_edge_relative_to_inertia = this->wing_rotation_relative_to_inertia * \
        this->virtural_Wing_Plane_relative_to_Wing * this->leading_edge_relative_to_wing;
    this->chordwise_direction_relative_to_inertia = this->wing_rotation_relative_to_inertia * \
        this->virtural_Wing_Plane_relative_to_Wing * this->chordwise_direction_relative_to_wing;
    return 1;
}

int BladeAeroCalculator::RequestVelocities(Eigen::Vector<data_cat, 3> vel_FreeFlow,\
    Eigen::Vector<data_cat, 3> vel_body_translation,\
    Eigen::Vector<data_cat, 3> vel_wing_rotation,\
    data_cat vel_AoA){
    this->v_infi = vel_FreeFlow;
    this->v_t = vel_body_translation;

    Eigen::MatrixX<data_cat> posLeadingEdgeVec_in_Inertia_frame = this->wing_rotation_relative_to_inertia * \
        (this->posLeadingEdgeVec.transpose());
    this->v_r = this->vec_mat_cross(vel_wing_rotation, posLeadingEdgeVec_in_Inertia_frame.transpose());
    
    this->v_AoA = vel_AoA;
    // this->v_r = 
    // this->posLeadingEdgeVec
    // const int right_number_of_bladeelements = this->Number_of_BladeElements;
    return 1;
}

///Calculations


int BladeAeroCalculator::CalcEffectiveVelocity(){
    Eigen::MatrixX<data_cat> minus_v_r = - this->v_r;

    // (Number_of_BladeElements X 3) - dim 
    Eigen::MatrixX<data_cat> rawResultantFlow = minus_v_r.rowwise() + (this->v_infi - this->v_t).transpose(); 
    this->last_added_mass_FV = this->added_mass_FV;
    this->eFV = rawResultantFlow - (rawResultantFlow * this->leading_edge_relative_to_inertia) * \
        (this->leading_edge_relative_to_inertia).transpose();
    this->rotation_eFV = minus_v_r- (minus_v_r * this->leading_edge_relative_to_inertia) * \
        (this->leading_edge_relative_to_inertia).transpose();
    this->added_mass_FV = this->rotation_eFV;
    
    return 1;
}

int BladeAeroCalculator::setlength_period( data_cat l_p){
    this->length_period = l_p;
    return 1;
}

int BladeAeroCalculator::CalcAoA(){
    Eigen::VectorX<data_cat> the_quotient = Eigen::VectorX<data_cat>(this->mat_SecNorm(this->eFV)).array() + SMALL_CONSTANT_4_NORMALIZATION;
    Eigen::VectorX<data_cat> cos_between_effective_velocity_and_chordwise_direction_relative_to_inertia = \
        (this->eFV * this -> chordwise_direction_relative_to_inertia).cwiseQuotient( the_quotient );
    this->AoA = cos_between_effective_velocity_and_chordwise_direction_relative_to_inertia.array().acos();
    return 1;
}

int BladeAeroCalculator::CopmputeAerodynamicForce(){
    Eigen::VectorX<data_cat> C_Lt = this->A_L * (2 * this->AoA).array().sin();
    Eigen::VectorX<data_cat> C_Dt = this->A_D * (- 1 * (2 * this->AoA).array().cos() + 1) + this->C_D_0;

    //##################### 1st Force amplitude component #####################

    Eigen::VectorX<data_cat> First_Force_component_2nd_token = \
        C_Lt.cwiseProduct( Eigen::VectorX<data_cat> (this->mat_SecNorm(this->eFV).array().pow(2)));

    /* (Number_of_BladeElements X 1) - dim */    
    Eigen::VectorX<data_cat> F_t_lift_amp = 0.5 * this->air_Density * this->BladeElementswidth * \
        ((this->y_data).cwiseAbs()).cwiseProduct(First_Force_component_2nd_token);

    //##################### 2nd Force amplitude component #####################

    Eigen::VectorX<data_cat> Second_Force_component_2nd_token = \
        C_Dt.cwiseProduct( Eigen::VectorX<data_cat> (this->mat_SecNorm(this->eFV).array().pow(2)));
    Eigen::VectorX<data_cat> F_t_drag_amp = 0.5 * this->air_Density * this->BladeElementswidth * \
        ((this->y_data).cwiseAbs()).cwiseProduct(Second_Force_component_2nd_token);

     //##################### 3rd Force amplitude component #####################

    this->d_eFV = (1 - this-> Filter_Constant) * this->d_eFV + \
        this-> Filter_Constant / this->length_period * (this->added_mass_FV - this->last_added_mass_FV);
    data_cat C_r = PI * (0.75 - this->x_hat_0);
    Eigen::MatrixX<data_cat> Sec_Norm_eFV = this->mat_SecNorm(this->rotation_eFV);
    Eigen::MatrixX<data_cat> squared_ydata = (this->y_data).cwiseProduct(this->y_data); 

    Eigen::VectorX<data_cat> F_r_amp = 0.5 * this->air_Density * C_r * this->BladeElementswidth * \
        abs(this->v_AoA) * Sec_Norm_eFV.cwiseProduct(squared_ydata);
    
    //##################### 4th Force amplitude component #####################
    Eigen::MatrixX<data_cat> EF_tokken_1_1 = ((this->eFV.cwiseProduct(this->d_eFV))).cwiseQuotient(\
        (Eigen::MatrixX<data_cat>)(((this->mat_SecNorm(this->eFV))).array()+SMALL_CONSTANT_4_NORMALIZATION).replicate(1,3));
    
    Eigen::VectorX<data_cat> EF_tokken_1   = ((Eigen::MatrixX<data_cat>)(EF_tokken_1_1.rowwise().sum())).cwiseProduct(\
        (Eigen::MatrixX<data_cat>)((this->AoA).array().sin()));

    Eigen::VectorX<data_cat> EF_tokken_2_1 = (this->mat_SecNorm(this->d_eFV)).cwiseProduct(\
        (Eigen::MatrixX<data_cat>)((- this->AoA).array()+PI));
    Eigen::VectorX<data_cat> EF_tokken_2   =  EF_tokken_2_1.cwiseProduct(\
        (Eigen::MatrixX<data_cat>)((this->AoA).array().cos()));
    Eigen::VectorX<data_cat> EF_tokken = EF_tokken_1 + EF_tokken_2;

    Eigen::VectorX<data_cat> F_a_amp = 0.25 * this->air_Density * this->BladeElementswidth * \
        PI * squared_ydata.cwiseProduct(EF_tokken);
    
    
    // #################### Directions computation #####################

    Eigen::MatrixX<data_cat> _1st_comp_direc = this->mat_Normalize( this->vec_mat_cross( \
        this->leading_edge_relative_to_inertia, this->eFV));

    for (int i = 0; i < _1st_comp_direc.rows(); i++){
        if ((_1st_comp_direc.row(i) * (this->chordwise_direction_relative_to_inertia))\
            (0,0)> 0.0){
               _1st_comp_direc.row(i) = -_1st_comp_direc.row(i); 
        }
    }

    Eigen::MatrixX<data_cat> _2nd_comp_direc = this->mat_Normalize(this->eFV);

    Eigen::Vector3<data_cat> _3rd_and_4th_comp_direc_single = (this->leading_edge_relative_to_inertia).cross(\
        this->chordwise_direction_relative_to_inertia);

    Eigen::MatrixX<data_cat> _3rd_and_4th_comp_direc =  (_3rd_and_4th_comp_direc_single.transpose())\
        .replicate(this->Number_of_BladeElements,1);
    for (int i = 0; i < _3rd_and_4th_comp_direc.rows(); i++){
        // cout<<"_3rd_and_4th_comp_direc.row(i)"<<_3rd_and_4th_comp_direc.row(i)<<endl;
        // cout<<"this->eFV.transpose()"<<this->eFV.transpose()<<endl;
        if ((_3rd_and_4th_comp_direc.row(i) * (this->eFV.transpose()))\
            (0,0)< 0.0){
               _3rd_and_4th_comp_direc.row(i) = -_3rd_and_4th_comp_direc.row(i); 
        }
    }

    data_cat sum_abs_of_F_t_lift_amp = F_t_lift_amp.cwiseAbs().sum();
    data_cat sum_abs_of_F_r_amp = F_r_amp.cwiseAbs().sum();
    data_cat sum_abs_of_F_a_amp = F_a_amp.cwiseAbs().sum();

    data_cat sign_sum_of_F_t_lift_amp = (F_t_lift_amp.sum() > 0) ? 1 : -1;
    data_cat sign_sum_of_F_r_amp = (F_r_amp.sum() > 0) ? 1 : -1;
    data_cat sign_sum_of_F_a_amp = (F_a_amp.sum() > 0) ? 1 : -1;

    Eigen::VectorX<data_cat> Weighted_F_t_lift_amp = 1 / this->Number_of_BladeElements * Eigen::VectorX<data_cat>::Ones(this->Number_of_BladeElements);
    Eigen::VectorX<data_cat> Weighted_F_r_amp = 1 / this->Number_of_BladeElements * Eigen::VectorX<data_cat>::Ones(this->Number_of_BladeElements);
    Eigen::VectorX<data_cat> Weighted_F_a_amp = 1 / this->Number_of_BladeElements * Eigen::VectorX<data_cat>::Ones(this->Number_of_BladeElements);

    if (sum_abs_of_F_t_lift_amp > SMALL_CONSTANT_4_NORMALIZATION){
        Weighted_F_t_lift_amp = sign_sum_of_F_t_lift_amp / (sum_abs_of_F_t_lift_amp + SMALL_CONSTANT_4_NORMALIZATION) *\
            F_t_lift_amp;
    }
    if (sum_abs_of_F_r_amp > SMALL_CONSTANT_4_NORMALIZATION){
        Weighted_F_r_amp = sign_sum_of_F_r_amp / (sum_abs_of_F_r_amp + SMALL_CONSTANT_4_NORMALIZATION) *\
            F_r_amp;
    }
    if (sum_abs_of_F_a_amp > SMALL_CONSTANT_4_NORMALIZATION){
        Weighted_F_a_amp = sign_sum_of_F_a_amp / (sum_abs_of_F_a_amp + SMALL_CONSTANT_4_NORMALIZATION) *\
            F_a_amp;
    }

    this->F_t_lift = (F_t_lift_amp.replicate(1,_1st_comp_direc.cols())).cwiseProduct(_1st_comp_direc);

    this->F_t_drag = F_t_drag_amp.replicate(1,_2nd_comp_direc.cols()).cwiseProduct(_2nd_comp_direc);
    this->F_r = F_r_amp.replicate(1,_1st_comp_direc.cols()).cwiseProduct(_1st_comp_direc);
    this->F_a = F_a_amp.replicate(1,_3rd_and_4th_comp_direc.cols()).cwiseProduct(_3rd_and_4th_comp_direc);
    
    data_cat raw_X_pos_t = (Weighted_F_t_lift_amp.cwiseProduct(this->x_data)).sum();
    data_cat raw_X_pos_r = (Weighted_F_r_amp.cwiseProduct(this->x_data)).sum();
    data_cat raw_X_pos_a = (Weighted_F_a_amp.cwiseProduct(this->x_data)).sum();

    data_cat raw_Y_pos_t = (Weighted_F_t_lift_amp.cwiseProduct(this->y_data.cwiseAbs())).cwiseProduct(\
        1 / PI * Eigen::VectorX<data_cat> (- this->AoA.array() + PI).cwiseAbs()).sum();
    data_cat raw_Y_pos_r = - 0.5 * (Weighted_F_r_amp.cwiseProduct(this->y_data.cwiseAbs())).sum();
    data_cat raw_Y_pos_a = - 0.5625 * (Weighted_F_a_amp.cwiseProduct(this->y_data.cwiseAbs())).sum();

    Eigen::Vector3<data_cat> raw_pos_t;
    raw_pos_t << raw_X_pos_t, raw_Y_pos_t, 0;
    Eigen::Vector3<data_cat> raw_pos_r;
    raw_pos_r << raw_X_pos_r, raw_Y_pos_r, 0;
    Eigen::Vector3<data_cat> raw_pos_a;
    raw_pos_a << raw_X_pos_a, raw_Y_pos_a, 0;

    Eigen::Vector3<data_cat> raw_pos_t_in_wing_frame = (this->virtural_Wing_Plane_relative_to_Wing.transpose())*\
        raw_pos_t;
    Eigen::Vector3<data_cat> raw_pos_r_in_wing_frame = (this->virtural_Wing_Plane_relative_to_Wing.transpose())*\
        raw_pos_r;
    Eigen::Vector3<data_cat> raw_pos_a_in_wing_frame = (this->virtural_Wing_Plane_relative_to_Wing.transpose())*\
        raw_pos_a;

    this->X_pos_t = raw_pos_t_in_wing_frame[0];
    this->Y_pos_t = raw_pos_t_in_wing_frame[1];

    this->X_pos_r = raw_pos_r_in_wing_frame[0];
    this->Y_pos_r = raw_pos_r_in_wing_frame[1];

    this->X_pos_a = raw_pos_a_in_wing_frame[0];
    this->Y_pos_a = raw_pos_a_in_wing_frame[1];

    return 1;


}