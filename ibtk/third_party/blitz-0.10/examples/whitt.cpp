#include <blitz/array.h>

using namespace blitz;       

class Base {
       public:
         virtual void multiply(Array<int,1>& array) const=0;
       };

       class Derived1 : public Base {
       public:
         Derived1(){}
         void multiply(Array<int,1>& array) const {array*=1;}
       };

       class Derived2 : public Base {
       public:
         Derived2(){}
         void multiply(Array<int,1>& array) const {array*=2;}
       };

       Array<int,1>& operator*=(Array<int,1>& array, const Derived1& object) {
         object.multiply(array);
         return array;
       }

int main()
{
       Array<int,1> array(5);
       array=1;

       Derived1 derived1;

       array*=derived1;
}
