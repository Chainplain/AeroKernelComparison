#include<iostream>
#include"Eigen/Core"
#include"Eigen/Dense"
#include"Eigen/Geometry"
#include<random>
#include<ctime>
#include<cmath>


#include"npy.hpp"
#include"BladeElementComputation.hpp"
#include"BladeElementComputation.cpp"

using namespace std;
Eigen::VectorX<data_cat> newnorm(Eigen::VectorX<data_cat> a){
    data_cat fac = 0.1 + a.norm();
    return fac * a;
 }

double generate_random(double dummy){
  static default_random_engine e(time(0));
  static normal_distribution<double> n(0,30);
  return n(e);
}


int main()
{
    int indicator = 0;
    int Compare_loop_length = 4000;

    BladeAeroCalculator Test_bac("SimpleFlapper0Wing40BLE_X_pos.npy",\
      "SimpleFlapper0Wing40BLE_y_pos.npy");

    Eigen::Matrix<data_cat, 3, 3> svwprw;
    svwprw<<1,0,0,0,1,0,0,0,1;
    Test_bac.SetVirturalWingPlaneRelative2Wing(svwprw);


    clock_t startTime,endTime;

    startTime = clock();
    while (indicator < Compare_loop_length){
        Eigen::Matrix<data_cat, 3, 3> LU_wing_Rotation_matrix;
        LU_wing_Rotation_matrix<<1,0,0,0,1,0,0,0,1;

        Eigen::Vector<data_cat, 3> vel_FreeFlow;
        vel_FreeFlow<<0,0,0;
        
        Eigen::Vector<data_cat, 6> LU_velVec;
        LU_velVec = Eigen::Vector<data_cat, 6>::Zero(6).unaryExpr(ptr_fun(generate_random));

        data_cat d_Real_motor_LU_wing_joint_sensor_value = generate_random(0.0);

        Test_bac.RequestWingOrientation2InertiaFrame(LU_wing_Rotation_matrix);
        Test_bac.CalcWingPlaneDirection();
        Test_bac.RequestVelocities(vel_FreeFlow, LU_velVec(Eigen::seq(0, 2)), LU_velVec(Eigen::seq(3, 5)), \
            d_Real_motor_LU_wing_joint_sensor_value);

        Test_bac.CalcEffectiveVelocity();
        Test_bac.CalcAoA();
        Test_bac.CopmputeAerodynamicForce();
        indicator++;
    }
    endTime = clock();
    cout<<"A single loop is complete!!"<<endl;
    cout<<"Comsume time of "<< (double(endTime-startTime))/CLOCKS_PER_SEC<<endl;
    // cout<< LU_velVec(Eigen::seq(0, 2)) <<endl;
    // cout<< LU_velVec(Eigen::seq(3, 5)) <<endl;
    // vector<unsigned long> shape {};
    // bool fortran_order;
    // vector<double> data_d;
  
    // const string path {"SimpleFlapper0Wing40BLE_X_pos.npy"};
    // npy::LoadArrayFromNumpy(path, shape, fortran_order, data_d);
    // vector<data_cat> data;
    // for (int i = 0; i < data_d.size(); i++){
    //     data.push_back( data_cat(data_d.at(i)));
    // }

    // Eigen::VectorX<data_cat> dara_eigen;
    
    // for (int i = 0; i < data.size(); i++)
    // {
    //     cout<< data[i]<<" ";
    // }
    // cout<< endl;

    // cout<< "Data size:" << data.size()<<endl;

    // dara_eigen = Eigen::Map<Eigen::VectorX<data_cat>, Eigen::Unaligned>(data.data(), data.size());

    // cout<< "dara_eigen:"<<dara_eigen<<endl;

    // Eigen::Index maxidex;
    // cout<< "max(dara_eigen):"<< dara_eigen.maxCoeff(&maxidex)<<endl;
    // cout<< "max(dara_eigen) Index:"<< maxidex<<endl;

    // Eigen::VectorX<data_cat> zeros = Eigen::VectorX<data_cat>::Zero(3);
    // Eigen::Matrix3<data_cat> mazeros;
    // mazeros<<zeros,zeros+ 12 * Eigen::VectorX<data_cat>::Ones(3), 1.1* Eigen::VectorX<data_cat>::Ones(3);
    
    // Eigen::Vector3<data_cat> test_cross, test_cross2;
    // test_cross<<1,-2,3; 
    // test_cross2<<1,2,3; 

    // cout<< "mazeros" << mazeros <<endl;
    // cout<< "mazeros.*mazeros" << mazeros.cwiseProduct(mazeros)<<endl;
    // cout<< "mazeros.*mazeros.sum" << mazeros.cwiseProduct(mazeros).rowwise().sum()<<endl;
    // cout<< "mazeros.*mazeros" << mazeros.cwiseProduct(mazeros).rowwise().sum().cwiseSqrt()<<endl;
    // cout<< "mazeros.*mazeros.sum.sqrt" << mazeros.cwiseProduct(mazeros).rowwise().sum().cwiseSqrt()<<endl;
    // cout<< "mazeros.*mazeros.sum.sqrt_" <<BladeAeroCalculator::mat_SecNorm(mazeros);
    
    // Eigen::MatrixX<data_cat> ss =  mazeros;

    // cout<<"ss"<<ss<<endl;
    // // int n = 1;
    // // assert(n == 2);


    // Eigen::MatrixX<data_cat> aa =  Eigen::MatrixX<data_cat>::Ones(3,3);
    // cout<<"aa:"<<aa<<endl;
    // cout<<"aa.normalize:"<<(aa.normalized())<<endl;
    // Eigen::MatrixX<data_cat> bb(6,3);
    // bb<<aa,aa;
    // cout<<"bb-test_cross:"<<bb.rowwise()-test_cross.transpose()<<endl;
    // Eigen::MatrixX<data_cat> hit = BladeAeroCalculator:: mat_vec_cross(bb, test_cross);
    // cout<<"cross:"<<hit<<endl;
    // Eigen::MatrixX<data_cat> git = BladeAeroCalculator:: vec_mat_cross(test_cross, bb);
    // cout<<"cross:"<<git<<endl;

    // cout<<"dot:"<<test_cross * (test_cross.transpose())<<endl;
    // cout<<"Quo:"<<test_cross.cwiseQuotient(test_cross) <<endl;
    // cout<<"Product:"<<test_cross.cwiseProduct(test_cross) <<endl;
    // cout<<"Scalar sum:"<<(test_cross.colwise().sum())<<endl;
    // cout<<"Scalar pow:"<<-(test_cross.array().pow(2))<<endl;
    // Eigen::VectorX<data_cat> testnormalize(4);
    // testnormalize<<1,2,3,4;
    // cout<<"vector normalize:"<<(testnormalize.normalized())<<endl;

    // Eigen::MatrixX<data_cat> ccc(3,3);
    // ccc<<1,1,1,2,2,6,7,8,9;
    // cout<<"ccc:"<< ccc <<endl;
    // Eigen::MatrixX<data_cat> jit =BladeAeroCalculator::mat_Normalize(ccc);
    // cout<<"matrix normalize:"<<jit<<endl;
    // // int index_ccc[]={2,3};
    // cout<<"Multi test_cross"<<test_cross.replicate(2,1)<<endl;
    // cout<<"Multi test_cross"<<test_cross.transpose().replicate(2,1)<<endl;
    // cout<<"Sum test_cross"<<test_cross.cwiseAbs().sum()<<endl;
    // cout<<"Pro test_cross"<<test_cross.replicate(1,3).cwiseProduct(ccc)<<endl;
    // cout<<"test_cross Pro test_cross"<<test_cross.cwiseProduct(test_cross)<<endl;
    // cout<<"test_cross[0][1][2]"<<test_cross[0]<<test_cross[1]<<test_cross[2]<<endl;
    // cout<<"sign"<<
    // cout<<"new"<<endl;
    // cout<<

    // Eigen::MatrixX<data_cat> bb =  Eigen::MatrixX<data_cat>::Identity(4,4);
    // cout<<"cross:"<<(Eigen::Vector4d (bb.row(1))).cross((Eigen::Vector4d) bb.row(2)) <<endl;
    // sout<<"test_cross cross test_cross"<<test_cross.cross(test_cross2)<<endl;
    // Eigen::Matrix<float, 3, 3>matrix_33;
    // Eigen::Matrix<data_cat, 3, 3> wing_rotation_relative_to_inertia;
    // wing_rotation_relative_to_inertia << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    // cout << wing_rotation_relative_to_inertia << endl;
    // float a = 0.000000001;
    // Eigen::Vector<data_cat, 3> d_eFV_1_mult_3;  
    // d_eFV_1_mult_3 << 1,2,3;  
    // cout<<wing_rotation_relative_to_inertia * d_eFV_1_mult_3<<endl;
    // Eigen::MatrixX<data_cat>  kk;
    // kk = * new Eigen::MatrixX<data_cat>(2,2);
    // kk <<1,2,3,4;
    // cout<<kk * kk.transpose()<<endl;
    // Eigen::Matrix<data_cat, 4, 2> AoA;
    // AoA << kk, kk;
    // cout<< 0.023 * AoA<<endl;
    // // kk(1,1) = 2;
    // // cout << a <<endl;
    // // system("pause");
    // Eigen::VectorX<data_cat> vv(4) ;
    // vv<<0,1,3,4;
    // cout<<"vv's norm:"<<newnorm(vv)<<endl;
    // Eigen::Matrix<data_cat, 2, 2> Svv;
    // Svv<<vv, vv;
    // cout<< "Svv"<<Svv<<endl;
    // int ll;
    // cin>>ll;
    return 0;
}

