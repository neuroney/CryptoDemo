#include <sstream>
#include <iostream>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

typedef struct
{
    int l;
    ZZ N2; // N^2
    int k;
    ZZ Bmsg;
    ZZ Bsk;

    ZZ N;
    Vec<ZZ> D;
    ZZ d;
    Vec<ZZ> d_Bsk;

    ZZ pk; // for PaillierEG
    ZZ g;  // for PaillierEG

} Para;

typedef struct
{
    ZZ kprf;
    Vec<ZZ> d_Bsk;
} EK;

void Paillier_Gen(Para &param, ZZ &d);
void Pailler_Enc(ZZ &ct, Para param, ZZ x);
void Pailler_Dec(ZZ &x, Para param, ZZ ct);
void Setup(Para &param, EK &ek0, EK &ek1);
void Input(Vec<ZZ> &I, Para param, ZZ x);
void ConvertInput(Vec<ZZ> &Mx, int sigma, Para param, EK ek, Vec<ZZ> Ix);
void Mul(Vec<ZZ> &Mz, Para param, Vec<ZZ> Ix, Vec<ZZ> My);
void DDLog(ZZ &z, Para param, ZZ g);
void AddMemory(Vec<ZZ> &Mz, Para param, Vec<ZZ> Mx, Vec<ZZ> My);
void AddInput(Vec<ZZ> &Iz, Para param, Vec<ZZ> Ix, Vec<ZZ> Iy);
void evaluate(ZZ &y_b_ZZ, int b, Para param, EK ekb, Vec<ZZ> Ix);

void Paillier_Gen(Para &param, ZZ &d)
{
    ZZ p, q, Phi_N, temp1, temp2;
    GenGermainPrime(p, 1036); // safe prime
    GenGermainPrime(q, 1036);
    mul(param.N, p, q);
    mul(param.N2, param.N, param.N);

    temp1 = p - 1;
    temp2 = q - 1;
    mul(Phi_N, temp1, temp2);
    InvMod(temp1, Phi_N, param.N);
    mul(temp1, temp1, Phi_N);
    mul(temp2, param.N, Phi_N);
    rem(d, temp1, temp2);
    param.d = d;

    param.k = 683;
    power(param.Bmsg, 2, param.k);
    power(param.Bsk, 2, param.k);
}

void Pailler_Enc(ZZ &ct, Para param, ZZ x)
{
    ZZ r, temp1, temp2;
    RandomBnd(r, param.N2);
    PowerMod(temp1, r, param.N, param.N2);
    temp2 = param.N + 1;
    PowerMod(temp2, temp2, x, param.N2);
    MulMod(ct, temp1, temp2, param.N2);
}

void Pailler_Dec(ZZ &x, Para param, ZZ ct)
{
    ZZ temp;
    PowerMod(temp, ct, param.d, param.N2);
    x = (temp - 1) / param.N;
}

void PaillierEG_Gen(Para &param, ZZ &d)
{
    ZZ p, q;
    GenGermainPrime(p, 1036); // safe prime
    GenGermainPrime(q, 1036);
    mul(param.N, p, q);
    mul(param.N2, param.N, param.N);
    RandomBnd(param.g, param.N2);
    PowerMod(param.g, param.g, 2 * param.N, param.N2);
    RandomBnd(d, param.N2);
    param.d = d;
    PowerMod(param.pk, param.g, d, param.N2);
}

void PaillierEG_Enc(Vec<ZZ> &ct, Para param, ZZ x)
{
    if (!ct.length())
    {
        ct.SetLength(2);
    }
    ZZ r;
    RandomBnd(r, param.N);
    PowerMod(ct[0], param.g, r, param.N2);
    PowerMod(r, param.pk, r, param.N2);
    PowerMod(ct[1], param.N + 1, x, param.N2);
    MulMod(ct[1], r, ct[1], param.N2);
}

void PaillerEG_Dec(ZZ &x, Para param, Vec<ZZ> ct)
{
    InvMod(x, ct[0], param.N2);
    PowerMod(x, x, param.d, param.N2);
    MulMod(x, ct[1], x, param.N2);
    divide(x, x - 1, param.N);
}

void Setup(Para &param, EK &ek0, EK &ek1)
{
    ZZ d, temp1, temp2;

    Paillier_Gen(param, d);
    temp1 = d;
    while (temp1 != 0)
    {
        // temp2 = temp1 % Bsk; temp1 = temp1 / Bsk
        DivRem(temp1, temp2, temp1, param.Bsk);
        param.d_Bsk.append(temp2);
    }
    param.l = param.d_Bsk.length();

    temp1 = param.Bmsg * param.Bsk;
    param.D.SetLength(param.l);
    for (int i = 0; i < param.l; ++i)
    {
        RandomBnd(temp2, temp1);
        ek1.d_Bsk.append(temp2);
        ek0.d_Bsk.append(temp2 - param.d_Bsk[i]);
        Pailler_Enc(param.D[i], param, param.d_Bsk[i]);
    }
}

void Input(Vec<ZZ> &I, Para param, ZZ x)
{
    ZZ X, r, temp1, temp2;
    Pailler_Enc(X, param, x);
    I.append(X);
    for (int i = 0; i < param.l; ++i)
    {
        RandomBnd(r, param.N2);
        temp1 = PowerMod(r, param.N, param.N2);
        temp2 = PowerMod(param.D[i], x, param.N2);
        I.append(MulMod(temp1, temp2, param.N2));
    }
}

void ConvertInput(Vec<ZZ> &Mx, int sigma, Para param, EK ek, Vec<ZZ> Ix)
{
    ZZ temp1;
    RandomBnd(temp1, param.Bmsg);
    Vec<ZZ> M1;

    if (sigma == 0)
    {
        M1.append(temp1);
    }
    else
    {
        AddMod(temp1, temp1, ZZ(1), param.N);
        M1.append(temp1);
    }
    M1.append(ek.d_Bsk);
    Mx.SetLength(param.l + 1);
    Mul(Mx, param, Ix, M1);
}

void Mul(Vec<ZZ> &Mz, Para param, Vec<ZZ> Ix, Vec<ZZ> My)
{
    if (!Mz.length())
    {
        Mz.SetLength(param.l + 1);
    }

    ZZ temp1, temp2;
    temp2 = 0;
    for (int i = 0; i < param.l; ++i)
    {
        temp1 = power(param.Bsk, i);
        temp1 *= My[i + 1];
        temp2 += temp1;
    }

    for (int i = 0; i < param.l + 1; ++i)
    {
        PowerMod(temp1, Ix[i], temp2, param.N2);
        DDLog(Mz[i], param, temp1);
    }
}

void DDLog(ZZ &z, Para param, ZZ g)
{
    ZZ h1, h, temp1;
    DivRem(h1, h, g, param.N); // h = g % N; h1 = g / N
    temp1 = InvMod(h, param.N);
    MulMod(z, h1, temp1, param.N);
}

void AddMemory(Vec<ZZ> &Mz, Para param, Vec<ZZ> Mx, Vec<ZZ> My)
{
    if (!Mz.length())
    {
        Mz.SetLength(param.l + 1);
    }
    for (int i = 0; i < param.l + 1; ++i)
    {
        Mz[i] = Mx[i] + My[i];
    }
}

void AddInput(Vec<ZZ> &Iz, Para param, Vec<ZZ> Ix, Vec<ZZ> Iy)
{
    if (!Iz.length())
    {
        Iz.SetLength(param.l + 1);
    }
    for (int i = 0; i < param.l + 1; ++i)
    {
        MulMod(Iz[i], Ix[i], Iy[i], param.N2);
    }
}

void evaluate(ZZ &y_b_ZZ, int b, Para param, EK ekb, Vec<ZZ> Ix)
{
    Vec<ZZ> x_b;
    ConvertInput(x_b, b, param, ekb, Ix); // E(x) -> x_0
    for (int i = 0; i < 1; ++i) {
        Mul(x_b, param, Ix, x_b); // E(x) * x_0 = x^2_0
    }
    
    y_b_ZZ = x_b[0];
}
