// Polymorphic arrays example

#include <blitz/array.h>

using namespace blitz;

class Material {
public:
    virtual double density() const = 0;
};

class Air {
public:
    virtual double density() const { return 0.291; }
};

class Water {
public:
    virtual double density() const { return 0.335; }
};

int main()
{
    Array<Material*,1> A(4);
}

