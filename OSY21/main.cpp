#include "tool.h"
int main(int, char **)
{
    Para param;
    EK ek0, ek1;
    Setup(param, ek0, ek1);

    ZZ x, y_0, y_1;
    x = 100;

    Vec<ZZ> Ix;
    Input(Ix, param, x);
    evaluate(y_0, 0, param, ek0, Ix);
    evaluate(y_1, 1, param, ek1, Ix);
    cout << y_1 - y_0 << endl;

    return 0;
}