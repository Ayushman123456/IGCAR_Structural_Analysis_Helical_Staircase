#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
using namespace std;
void input(string &type, double &R1, double &R2, double &b, double &h, double &beta_val, double &phi_val, double &ri, double &height, double &w, double &finish_thickness, double &material, double &material_wt, double &t, double &ro)
{

    cout << "Enter the live load of the building in KN/m" << endl;
    cin >> w;
    cout << "Enter the height to be covered by the staircase(in m)" << endl;
    cin >> height;
    cout << "Enter the width of the staircase(in mm)" << endl;
    cin >> b;
    if (b < 1500)
    {

        cout << "Entered value of width is less than minimum height. It has been increased to 1500mm for calculations" << endl;
        b = 1500;
    }
    if (type == "Assembly Buildings")
    {
        if (b < 2000)
        {

            cout << "Entered value of width is less than minimum height. It has been increased to 2000mm for calculations" << endl;
            b = 2000;
        }
    }
    cout << "Enter the rise of the staircase (in mm)" << endl;
    cin >> h;
    if (h > 150)
    {
        cout << "Entered value of riser is more than maximum standard height. It has been reduced to 150mm for calculations" << endl;
        h = 150;
    }
    cout << "Enter the inner radius of the staircase(in mm)" << endl;
    cin >> ri;
    cout << "Enter the value of beta in degrees(the total angle inscribed by the staircase)" << endl;
    cin >> beta_val;
    cout << "Enter the value of the angle made with the horizontal in degrees(phi)" << endl;
    cin >> phi_val;
    if (phi_val > 45)
    {
        cout << "Entered value of phi is more than maximum standard pitch. It has been reduced to 45 degree for calculations" << endl;
        phi_val = 45;
    }
    cout << "Enter the thickness of the finishes(in mm)" << endl;
    cin >> finish_thickness;
    cout << "Enter the compressive strength of the concrete(as M of concrete)" << endl;
    cin >> material;
    cout << "Enter the thickness of the waist(in mm)" << endl;
    cin >> t;
    cout << "Enter the material weight for the finishing material(in kN/m2)" << endl;
    cin >> material_wt;

    ro = ri + b;
    R1 = 2 * ((pow(ro, 3) - pow(ri, 3)) / (3 * (pow(ro, 2) - pow(ri, 2)))) / 1000;
    R2 = 0.5 * (ri + ro) / 1000;
}

int main()
{
    ofstream myfile;
    myfile.open("test.csv");

    int z = 0;
    int a = 0;
    string type;
    double R1, R2, b, h, beta_val, phi_val, ri, height, w, finish_thickness, material, material_wt, t, ro;
    input(type, R1, R2, b, h, beta_val, phi_val, ri, height, w, finish_thickness, material, material_wt, t, ro);

    double beta[] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360};
    double phi_k2[] = {40, 35, 30, 25, 20};
    double phik1k3[] = {40, 20};

    // writing the data
    // 1st graph
    double k2a[] = {0.24, 0, 0.27, 0, 0.28, 0.46, 0.52, 0.58, 0.64, 0.69, 0.55, 0.65, 0.75, 0.87, 1, 0.61, 0.72, 0.85, 1.02, 1.22, 0.65, 0.77, 0.92, 1.11, 1.38, 0.69, 0.83, 0.99, 1.21, 1.51, 0.73, 0.88, 1.07, 1.31, 1.64, 0.80, 0.95, 1.14, 1.41, 1.78, 0.85, 1.02, 1.23, 1.51, 1.91, 0.91, 1.09, 1.32, 1.63, 2.07, 0.97, 1.16, 1.41, 1.74, 2.22, 1.02, 1.23, 1.51, 1.87, 2.36};
    double k1a[] = {-0.01, 0, -0.01, 0, -0.01, 0, -0.01, 0, -0.01, 0, -0.02, -0.01, -0.03, -0.03, -0.07, -0.06, 0.12, -0.11, -0.2, -0.18, -0.31, -0.29, -0.45, -0.41};
    double k3a[] = {-0.02, -0.02, -0.02, -0.05, -0.02, -0.07, -0.03, -0.08, -0.05, -0.1, -0.09, -0.14, -0.15, -0.19, -0.26, -0.3, -0.45, -0.48, -0.7, -0.72, -1.06, -1.09, -1.57, -1.59};

    // 2nd graph
    double k2b[] = {0.07, 0, 0.09, 0, 0.09, 0.25, 0, 0.3, 0, 0.31, 0.41, 0.46, 0.5, 0, 0.54, 0.51, 0.58, 0.66, 0.74, 0.81, 0.59, 0.69, 0.79, 0.92, 1.05, 0.65, 0.76, 0.9, 1.07, 1.26, 0.71, 0.83, 1, 1.19, 1.46, 0.77, 0.91, 1.09, 1.32, 1.63, 0.83, 0.98, 1.18, 1.45, 1.81, 0.89, 1.06, 1.29, 1.57, 1.97, 0.96, 1.15, 1.39, 1.7, 2.15, 1.02, 1.22, 1.47, 1.82, 2.31};
    double k1b[] = {0, 0, 0.02, 0.03, 0.03, 0.06, 0.02, 0.08, 0.03, 0.07, 0, 0.06, -0.01, 0.03, -0.05, -0.02, -0.11, -0.08, -0.19, -0.15, -0.3, -0.26, -0.44, -0.39};
    double k3b[] = {-0.01, -0.01, -0.05, -0.08, -0.08, -0.15, -0.09, -0.21, -0.11, -0.25, -0.14, -0.29, -0.2, -0.34, -0.31, -0.43, -0.47, -0.56, -0.73, -0.81, -1.11, -1.17, -1.58, -1.61};

    // 3rd graph
    double k2c[] = {0.36, 0.00, 0.37, 0.00, 0.41, 0.60, 0.69, 0.79, 1.12, 0.91, 0.71, 0.82, 0.96, 0.86, 1.28, 0.75, 0.90, 1.06, 1.27, 1.52, 0.79, 0.94, 1.12, 1.37, 1.69, 0.83, 0.99, 1.19, 1.45, 1.82, 0.88, 1.04, 1.26, 1.54, 1.95, 0.92, 1.10, 1.33, 1.64, 2.08, 0.99, 1.17, 1.41, 1.74, 2.22, 1.05, 1.25, 1.51, 1.86, 2.37, 1.11, 1.33, 1.60, 1.98, 2.52, 1.16, 1.39, 1.68, 2.09, 2.65};
    double k1c[] = {0.00, 0.00, 0.00, 0.01, -0.01, 0.00, -0.03, -0.03, -0.04, -0.04, -0.08, -0.08, -0.12, -0.12, -0.17, -0.17, -0.24, -0.24, -0.33, -0.32, -0.47, -0.44, -0.61, 0.58};
    double k3c[] = {-0.03, -0.03, -0.04, -0.08, -0.04, -0.10, -0.04, -0.11, -0.07, -0.13, -0.12, -0.17, -0.20, -0.25, -0.33, -0.37, -0.53, -0.57, -0.84, -0.87, -1.26, -1.30, -1.60, -1.60};

    // 4th graph
    double k2d[] = {0.16, 0, 0.17, 0, 0.14, 0.39, 0, 0.42, 0, 0.38, 0.53, 0.59, 0.65, 0, 0.71, 0.65, 0.74, 0.84, 0.93, 1.01, 0.72, 0.83, 0.97, 1.13, 1.28, 0.79, 0.91, 1.08, 1.28, 1.52, 0.84, 0.99, 1.17, 1.41, 1.71, 0.9, 1.07, 1.27, 1.55, 1.9, 0.96, 1.14, 1.38, 1.67, 2.08, 1.03, 1.22, 1.48, 1.8, 2.26, 1.09, 1.31, 1.58, 1.93, 2.43, 1.15, 1.37, 1.66, 2.05, 2.59};
    double k1d[] = {0.01, 0.01, 0.02, 0.02, 0.02, 0.03, -0.01, 0.03, -0.03, 0.02, -0.05, -0.01, -0.1, -0.06, -0.15, -0.12, -0.22, -0.19, -0.32, -0.29, -0.44, -0.41, -0.61, -0.56};
    double k3d[] = {0.02, -0.03, -0.07, -0.11, -0.1, -0.19, -0.12, -0.26, -0.14, -0.31, -0.18, -0.35, -0.25, -0.42, -0.37, -0.52, -0.57, -0.69, -0.88, -0.96, -1.32, -1.37, -1.61, -1.61};

    // for 1st graph
    double matrixk1a[12][2];
    double matrixk2a[12][5];
    double matrixk3a[12][2];

    // for 2nd graph
    double matrixk1b[12][2];
    double matrixk2b[12][5];
    double matrixk3b[12][2];

    // for 3rd graph
    double matrixk1c[12][2];
    double matrixk2c[12][5];
    double matrixk3c[12][2];

    // for 4th graph
    double matrixk1d[12][2];
    double matrixk2d[12][5];
    double matrixk3d[12][2];

    // fillind the matrices as representations of digitized graphs
    for (int i = 0; i < 12; i++)
    {
        for (int k = 0; k < 2; k++)
        {
            matrixk1a[i][k] = k1a[z];
            matrixk3a[i][k] = k3a[z];
            matrixk1b[i][k] = k1b[z];
            matrixk3b[i][k] = k3b[z];
            matrixk1c[i][k] = k1c[z];
            matrixk3c[i][k] = k3c[z];
            matrixk1d[i][k] = k1d[z];
            matrixk3d[i][k] = k3d[z];
            z++;
        }

        for (int j = 0; j < 5; j++)
        {
            matrixk2a[i][j] = k2a[a];
            matrixk2b[i][j] = k2b[a];
            matrixk2c[i][j] = k2c[a];
            matrixk2d[i][j] = k2d[a];
            a++;
        }
    }

    double swt_waist, swt_step, swt_finish, total_dw, ratio1, ratio2, Mo, K, Mvs, n, ul_wt;
    swt_waist = t * material / 1000;
    swt_step = h / 2 * phi_val / 1000;
    swt_finish = finish_thickness / 1000 * material_wt;
    total_dw = swt_waist + swt_step + swt_finish;
    ul_wt = (1.5 * total_dw) + (1.5 * w);
    n = b * ul_wt / 1000;
    ratio1 = R1 / R2;
    ratio1 = 1.061;
    ratio2 = b / h;

    // initializing secondary variables
    double min_diff, k1_40, k1_20, k3_40, k3_20, k11, k21, k31, k13, k23, k33, k1, k2, k3, k1_13, k2_13, k3_13, k1_24, k2_24, k3_24, k14, k24, k34, k12, k22, k32;
    int beta_pos, flag, beta_min_pos, phik2_min_pos, phi_pos, beta_final_pos, phi_final_pos;
    min_diff = INT_MAX;

    for (int i = 0; i < 12; i++)
    {
        if (beta_val == beta[i])
        {
            beta_pos = i;
            flag = 1;
        }
    }
    if (flag == 0)
    {
        for (int i = 0; i < 11; i++)
        {
            if (beta[i] < beta_val && beta[i + 1] > beta_val && (abs(beta[i] - beta_val) < min_diff))
                beta_min_pos = i;
            min_diff = abs(beta[i] - beta_val);
        }
    }
    flag = 0;
    min_diff = INT_MAX;
    for (int i = 0; i < 5; i++)
    {
        if (phi_val == phi_k2[i])
        {
            phi_pos = i;
            flag = 1;
        }
    }
    if (flag == 0)
    {
        for (int i = 0; i < 4; i++)
        {
            if (phi_k2[i] < phi_val && phi_k2[i + 1] > phi_val && (abs(phi_k2[i] - phi_val) < min_diff))
                phik2_min_pos = i;
            min_diff = abs(phi_k2[i] - phi_val);
        }
    }

    if ((ratio1 >= 1.0 && ratio1 <= 1.1) && (5 <= ratio2 <= 10))
    {
        if (beta_pos != 0)
        {
            beta_final_pos = beta_pos;
        }
        if (phi_pos != 0)
        {
            phi_final_pos = phi_pos;
        }

        // get the values from the first graph for R1/R2=1 and b/h=10;

        k1_40 = matrixk1a[beta_final_pos][0];
        k1_20 = matrixk1a[beta_final_pos][1];
        k3_40 = matrixk3a[beta_final_pos][0];
        k3_20 = matrixk3a[beta_final_pos][1];

        // values of k1,k2,k3 from first graph

        k21 = matrixk2a[beta_final_pos][phi_final_pos];
        k11 = k1_20 - (phi_val - 20) * ((k1_40 - k1_20) / 20);
        k31 = k3_20 + (phi_val - 20) * ((k3_40 - k3_20) / 20);

        // get the values from the third graph R1/R2=1.1 and b/h=10;
        k1_40 = matrixk1c[beta_final_pos][0];
        k1_20 = matrixk1c[beta_final_pos][1];
        k3_40 = matrixk3c[beta_final_pos][0];
        k3_20 = matrixk3c[beta_final_pos][1];

        // values of k1,k2,k3 from first graph
        k23 = matrixk2c[beta_final_pos][phi_final_pos];
        k13 = k1_20 - (phi_val - 20) * ((k1_40 - k1_20) / 20);
        k33 = k3_20 + (phi_val - 20) * ((k3_40 - k3_20) / 20);

        // finding values of k1,k2,k3 from first and 3rd graph(b/h=10)
        k1_13 = k13 - (1.1 - ratio1) * ((k13 - k11) / 0.1);
        k2_13 = k23 - (1.1 - ratio1) * ((k23 - k21) / 0.1);
        k3_13 = k33 + (1.1 - ratio1) * ((k33 - k31) / 0.1);

        // get values from the second graph
        k1_40 = matrixk1b[beta_final_pos][0];
        k1_20 = matrixk1b[beta_final_pos][1];
        k3_40 = matrixk3b[beta_final_pos][0];
        k3_20 = matrixk3b[beta_final_pos][1];

        // values of k1,k2,k3 from 2nd graph
        k22 = matrixk2b[beta_final_pos][phi_final_pos];
        k12 = k1_20 - (phi_val - 20) * ((k1_40 - k1_20) / 20);
        k32 = k3_20 + (phi_val - 20) * ((k3_40 - k3_20) / 20);

        // get values from the 4th graph
        k1_40 = matrixk1d[beta_final_pos][0];
        k1_20 = matrixk1d[beta_final_pos][1];
        k3_40 = matrixk3d[beta_final_pos][0];
        k3_20 = matrixk3d[beta_final_pos][1];

        // values of k1,k2,k3 from the 4th graph
        k24 = matrixk2b[beta_final_pos][phi_final_pos];
        k14 = k1_20 - (phi_val - 20) * ((k1_40 - k1_20) / 20);
        k34 = k3_20 + (phi_val - 20) * ((k3_40 - k3_20) / 20);

        // k1,k2,k3 values from the 2nd and 4th graph(b/h=5)
        k1_24 = k13 - (1.1 - ratio1) * ((k14 - k12) / 0.1);
        k2_24 = k23 - (1.1 - ratio1) * ((k24 - k22) / 0.1);
        k3_24 = k3_20 + (1.1 - ratio1) * ((k34 - k32) / 0.1);

        // now finding the k1,k2,k3 values for any values of b/h and R1/R2
        k1 = k1_13 - (10 - ratio2) * ((k1_24 - k1_13) / 5);
        k2 = k2_13 + (10 - ratio2) * ((k1_24 - k1_13) / 5);
        k3 = k3_13 - (10 - ratio2) * ((k1_24 - k1_13) / 5);
    }

    // calcualting the secondary variables
    // double Mo, H, Mvs;
    double H;
    Mo = k1 * n * pow(R2, 2);
    H = k2 * n * R2;
    Mvs = k3 * n * pow(R2, 2);

    double vertical_moment[36];
    double lateral_moment[36];
    double torsion[36];
    double axial_force[36];
    double vertical_shear[36];
    double Radial_horizontal_shear[36];

    double theta = 0;
    double phi;

    a = 1;
    for (int i = 10; i <= 360; i += 10)
    {
        theta = 3.14159 / 180 * i;
        phi = 3.14159 / 180 * phi_val;
        vertical_moment[a] = (Mo * cos(theta)) + (H * R2 * theta * tan(phi) * sin(theta)) - (n * pow(R1, 2) * (1 - cos(theta)));
        lateral_moment[a] = (Mo * sin(theta) * sin(phi)) - (H * R2 * theta * tan(phi) * cos(theta) * sin(phi)) - (H * R2 * sin(theta) * cos(phi)) + (n * R1 * sin(phi) * ((R1 * sin(theta)) - (R2 * theta)));
        torsion[a] = ((Mo * sin(theta)) - (H * R2 * theta * cos(theta) * tan(phi)) + (n * pow(R1, 2) * sin(theta)) - (n * R1 * R2 * theta)) * cos(phi) + (H * R2 * sin(theta) * sin(phi));
        axial_force[a] = (-H * sin(theta) * cos(phi)) - (n * R1 * theta * sin(phi));
        vertical_shear[a] = (n * R1 * theta * cos(phi)) - (H * sin(theta) * sin(phi));
        Radial_horizontal_shear[a] = H * cos(theta);
        a++;
    }

    cout << " " << endl;
    myfile << "Angle theta"
           << ","
           << ","
           << "Vertical Moment"
           << ","
           << ","
           << "Lateral Moment"
           << ","
           << ","
           << "Torsion"
           << ","
           << ","
           << "Axial Force"
           << ","
           << ","
           << "Vertical Shear"
           << ","
           << ","
           << "Radial Horizontal shear" << endl;
    myfile << 0 << ","
           << "," << (Mo * cos(0)) + (H * R2 * 0 * tan(phi) * sin(0)) - (n * pow(R1, 2) * (1 - cos(0))) << ","
           << "," << 0 << ","
           << "," << 0 << ","
           << "," << 0 << ","
           << "," << 0 << ","
           << "," << H * cos(0) << endl;
    int angle = 10;
    for (int i = 1; i < a; i++)
    {
        myfile << angle << ","
               << "," << vertical_moment[i] << ","
               << "," << lateral_moment[i] << ","
               << "," << torsion[i] << ","
               << "," << axial_force[i] << ","
               << "," << vertical_shear[i] << ","
               << "," << Radial_horizontal_shear[i] << endl;
        angle = angle + 10;
    }
    myfile << " " << endl;
    if (height * 1000 / h > 12)
    {
        myfile << "Intermediate landing is required at height" << (12 * h) << endl;
    }
    myfile << " " << endl;
    if (b <= 1500)
    {
        myfile << "Handrail has to be provided on one side of the staricase" << endl;
    }
    if (b > 1500)
    {
        myfile << "Handrail has to be provided on both sides of the staricase" << endl;
    }

    return 0;
}