/* mystructs.hxx

Hier baue ich meine Structs zusammen für Dinge, die mir sinnvoll erscheinen wie zum Beispiel Vektoren.
Es werden auch nützliche Operatoren überladen.

-------------------------------------------------------------------------------------------------------*/

#include <math.h>

struct vektor {
  double x;
  double y;
  double z;
  
  void set(double xin, double yin, double zin)
  {
    x = xin;
    y = yin;
    z = zin;
  };
  
  inline double norm()
  {
    return sqrt(x*x + y*y + z*z);
  };
  
  vektor cross(vektor const &other){
    vektor erg;
    erg.x = y*other.z-z*other.y;
    erg.y = z*other.x-x*other.z;
    erg.z = x*other.y-y*other.x;
    return erg;
  };
  
  vektor operator+=(vektor &other){
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  };
  
  vektor operator-(vektor &other){
    vektor erg;
    erg.x = x - other.x;
    erg.y = y - other.y;
    erg.z = z - other.z;
    return erg;
  };
  
  vektor operator*(double &other){
    vektor erg;
    erg.x = x*other;
    erg.y = y*other;
    erg.z = z*other;
    return erg;
  };
  
  vektor operator/(double &other){
    vektor erg;
    erg.x = x/other;
    erg.y = y/other;
    erg.z = z/other;
    return erg;
  };
  
  vektor operator/=(double &other){
    x /= other;
    y /= other;
    z /= other;
    return *this;
  };
};

