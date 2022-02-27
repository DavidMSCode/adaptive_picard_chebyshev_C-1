#ifndef ORBIT_H
#define ORBIT_H

#include <vector>

class Orbit
{
    private:
        std::vector<std::vector<double> > Soln;
        struct Properties{
            double Area;
            double Reflectivity;
            double Mass;
        };
        
    public:
        //Constructors
        Orbit(std::vector<std::vector<double > > Solution);
        Orbit(double area, double reflectivity, double mass);
        //Setters
        void SetSolution(std::vector<std::vector<double > > Solution);
        void SetProperties(double area, double reflectivity, double mass);
        //Getters
        std::vector<double> getTimes(){return Soln[0];};
        std::vector<double> getPositionX(){return Soln[1];};
        std::vector<double> getPositionY(){return Soln[2];};
        std::vector<double> getPositionZ(){return Soln[3];};
        std::vector<double> getVelocityX(){return Soln[4];};
        std::vector<double> getVelocityY(){return Soln[5];};
        std::vector<double> getVelocityZ(){return Soln[6];};
        //Internal properties struct declaration
        struct Properties properties;
};

#endif
