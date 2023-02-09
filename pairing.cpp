#include <sstream>
#include <iostream>
#include <NTL/ZZ.h>
#include <gmp.h>
extern "C" {
    #include <relic/relic.h>
}
using namespace std;
using namespace NTL;

void ZZ2bn(bn_t out, ZZ in)
{
    if(in==0)
    {
        bn_read_str(out,"0",1,10);
    }else
    {
        std::stringstream buffer;
        buffer << in;
        const char* str=strdup(buffer.str().c_str());
        bn_read_str(out,str,floor(log(in)/log(10))+1,10);
    }
}

void bntoZZ(ZZ out, bn_t in)
{
    int size = bn_size_str(in, 10);
    char *in_char = new char[size];
    bn_write_str(in_char, size, in, 10);
    out = conv<ZZ>(in_char);
}

int main(int, char **)
{
    core_init();
    ep_t g1_gen;
    ep2_t g2_gen;
    fp12_t gT_gen;

    ep_new(g1_gen);
    ep2_new(g2_gen);
    fp12_new(gT_gen);

    ep_curve_init();
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(RLC_EP_MTYPE);
    ep_curve_get_gen(g1_gen);
    ep2_curve_get_gen(g2_gen);
    pp_map_oatep_k12(gT_gen, g1_gen, g2_gen);

    ZZ g1_order_ZZ, g2_order_ZZ;
    bn_t g1_order, g2_order;
    bn_new(g1_order);
    bn_new(g2_order);
    ep_curve_get_ord(g1_order);
    ep2_curve_get_ord(g2_order);
    bntoZZ(g1_order_ZZ, g1_order);
    bntoZZ(g2_order_ZZ, g2_order);

    ZZ A_ZZ;
    bn_t A;
    ep_t g1_A;
    fp12_t gT_A;

    bn_new(A);
    ep_new(g1_A);
    fp12_new(gT_A);
    // RandomBnd(A, g2_order_ZZ);
    A_ZZ = 5;
    ZZ2bn(A, A_ZZ);
    ep_mul_gen(g1_A, A);               // g_1^A
    pp_map_oatep_k12(gT_A, g1_A, g2_gen); // g_T^A=e(g_1^A,g_2)


    ZZ y_0_ZZ, y_1_ZZ, Y_0_ZZ, Y_1_ZZ, y_ZZ;
    bn_t y,Y_0,Y_1;
    bn_new(y);
    bn_new(Y_0);
    bn_new(Y_1);

    y_0_ZZ = 10000; 
    y_1_ZZ = 80000;
    Y_0_ZZ = 350000;
    Y_1_ZZ = 700000;
    y_ZZ = y_1_ZZ - y_0_ZZ;
    ZZ2bn(y, y_ZZ);
    ZZ2bn(Y_0, Y_0_ZZ);
    ZZ2bn(Y_1, Y_1_ZZ);

    ep_t g1_Y0,g1_Y1;
    ep_new(g1_Y0);
    ep_new(g1_Y1);
    ep_mul_gen(g1_Y0, Y_0);
    ep_mul_gen(g1_Y1, Y_1);   

    fp12_t gT_Y0, right_side, left_side;
    fp12_new(gT_Y0);
    fp12_new(right_side);
    fp12_new(left_side);

    fp12_exp(right_side, gT_A, y); // g_T^Ay
    pp_map_oatep_k12(gT_Y0, g1_Y0, g2_gen); // g_T^Y_0 
    fp12_mul(right_side, right_side, gT_Y0); // g_T^Ay+Y0

    pp_map_oatep_k12(left_side, g1_Y1, g2_gen); // g_T^Y1
    if (fp12_cmp(left_side, right_side) == RLC_EQ)
    {
        cout << "Verification passed!\nThe value of y= " << y_ZZ << "\n";
    }
    else
    {
        printf("******************** ERROR ********************\n");
        fp12_print(left_side);
        printf("\n\n");
        fp12_print(right_side);
    }

    return 0;
}
