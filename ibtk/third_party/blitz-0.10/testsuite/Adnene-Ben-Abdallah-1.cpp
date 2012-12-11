#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<double,1> psiElem(3);
    psiElem =
      -0.121990517304243,-0.099118834514638,-0.164750336116265;

    double psi_max = blitz::max(psiElem);
    cout << "psiElem = " << psiElem << endl;
    cout << "psi_max = " << psi_max << endl;
}

