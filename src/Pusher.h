/**
 * \file    Pusher.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March 2019
 *
 * \brief   Declares the Pusher class.
 *
 *  
 */

#ifndef PUSHER_H_INCLUDED
#define PUSHER_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <stdexcept>
#include <cassert>
#include <sstream>
#include <fstream>
#include <string>
#include <omp.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "Particle.h"
#include "Vector.h"
// #include "Orbit.h"
#include "Constants.h"


template <class T> // templated to be used with Orbit or Mirror class geometry
class Pusher{

    public:
        Pusher(T& geo);
        ~Pusher();

        /**
         * \brief Pushes a single particle, returns final position
         * \param part The particle to be pushed
         * \param dt time step
         * \param iter Number of iterations
         * \bool write Whether the trajectory should be written to file
         */
        Vector pushSingle(Particle& part, double dt, int iter, bool write, std::ofstream& coord);

        /**
         * \brief Single particle pusher in cylindrical geometry
         * \note only part.move is changed for now. Maybe there's more?
         */
        Vector pushSingleCyl(Particle& part, double dt, int iter, bool write, std::ofstream& coord);

//-------------------------------------  Mirror   --------------------------------------------

        /**
         * \brief Simple particle push, with particles sourced at x=0
         */
        void midplaneBurst(double temperature, int spec, int nparts, int maxiter, bool write);

        /**
         * \brief A parallel particle pusher, collects lost and trapped particles
         */
        double losscone(double energy, bool spec, int nparts, int maxiter, bool write, std::string& suffix);

        /**
         * \brief Generates velocity vector from a Gaussian distribution.
         */
        Vector gaussian(double center, double vbar, std::default_random_engine& generator);

        /**
         * \brief Generates velocity vector from a flux distribution.
         */
        Vector flux(double center, double vbar, std::default_random_engine& generator);

//---------------------------------------  NBI   ---------------------------------------------


        /**
         * \brief burst out particles from spherical surface
         * \param radius Radius of curvature for surface source
         */
        void gridBurst(double radius, double ylim, int nsources, bool write);

        /**
         * \brief Burst out cones of particles from spherical surface. Basically a better 
         *        version of gridBurst
         * \param radius Radius of curvature for spherical source
         * \param ylim   half height of source plate
         * \param dtheta Divergence of gaussian cones
         * \param nsources Number of particle sources
         * \param partPerS Number of particles per source
         * \param write  Whether to write things to file
         *
         */
        void conicBurst(double radius, double ylim, double dtheta, int nsources, int partPerS, bool write);

        /**
         * \brief Generate a random pos vector on sphere surface
         * \note  Assumes that the canter of the sphere is at (radius, 0, 0)
         */
        Vector sphere(double radius, double ylim, std::default_random_engine& generator);


        /**
         * \brief Generate a normal vector on sphere surface, pointing towards center
         * \note  Assumes that the canter of the sphere is at (radius, 0, 0)
         */
        Vector sphereNormal(double radius, Vector pos);

        /**
         * \brief  Find a diverged vector on the conic surface. 
         *
         * \note   Calculates normal vector from radius and pos, generate a random vector on the tangent plane,
         *         then add to normal vector to find the diverged vector.
         *
         */
        Vector diverge(double radius, Vector& pos, double dtheta, std::default_random_engine& generator);

    private:
        T * geo_; // a pointer to the geometry class object

};

template <class T>
Pusher<T>::Pusher(T& geo)
        :geo_(nullptr)
{
    geo_ = &geo;
    return;
}

template <class T>
Pusher<T>::~Pusher()
{
    return;
}

template <class T>
Vector Pusher<T>::pushSingle(Particle& part, double dt, int iter, bool write, std::ofstream& coord)
{
    //coord.open("coordinates.out");
    //coord << std::setprecision(10);
    int lastCrossed = -1;// start at the last sightline (first to cross)
    for (int i = 0; i < iter; i++){
    //std::cerr << "pushing iteration: " << i << std::endl;    
        Vector posNow = part.pos();
        Vector BNow = (*geo_).getB(posNow);
        Vector ENow = (*geo_).getE(posNow);
        if (write && (i % 2 == 0)){
            coord << part.pos() << std::endl;
            double Ex = pow(part.vel().x(), 2) * 0.5 * part.mass()/QE;
	    double Eperp = part.energy() - Ex;    
	std::cout << Ex << "," << Eperp << "," << part.energy() << std::endl;
        }
            part.move(ENow, BNow, dt);
        //int lastCrossed = 26;// start at the last sightline (first to cross)
//		std::cerr << "E, "<< ENow << ", pos, " << posNow << ", v, " << part.vel() << std::endl;
            if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
            // std::cerr << "particle lost to limiter after" << i \
            // << "iterations." << std::endl;
                break;
            } else {
            // Next line for Mirror only
                    Vector position = part.pos();
                    (*geo_).addToBin(position);
            }
            if (write){ // always collect for density calculations
            // Next line for NBI only
            // lastCrossed = (*geo_).sightline(part, lastCrossed);
            }    
        }    
    return part.pos();
}

template <class T>
Vector Pusher<T>::pushSingleCyl(Particle& part, double dt, int iter, bool write, std::ofstream& coord)
{
    for (int i = 0; i < iter; i++){
        Vector posNow = part.pos();
            Vector BNow = (*geo_).getB(posNow);
            Vector ENow = (*geo_).getE(posNow);
        if (write && (i % 100 == 0)){
                    coord << part.pos() << std::endl;
            }
            part.moveCyl(ENow, BNow, dt);

            if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
            //std::cerr << "particle lost to limiter after" << i \
            << "iterations." << std::endl;
                break;
            }    
            if (write){
                // (*geo_).sightline(part);
            }
        }
    return part.pos();
}


//--------------------------------------------------------------------------------------------
//--------------------------- End of Basic Class Implementations -----------------------------
//-------------------------- Start geography specific calculations ---------------------------
//--------------------------------------------------------------------------------------------

//-------------------------------------  Mirror   --------------------------------------------
template <class T>
void Pusher<T>::midplaneBurst(double temperature, int spec, int nparts, int maxiter, bool write)
{
    std::ofstream coord;
    coord.open("midplaneBurst.out", ios::app);
    coord << std::setprecision(10);

    Doub mass = MI * (1 - spec) + ME * spec;
    Doub vbar = sqrt(temperature  *  EVTOJOULE / mass); // thermal velocity, <v^2> in distribution, sigma.

    // decoupling dt from magnetic field
    Doub Btypical = BMAGAXIS; // magnetic field on the order of unity Tesla.
    Doub fLamor = ( 1520 * (1 - spec) + 2.8E6 * spec ) * Btypical; // another logical, constants from NRL p28 ( Btypical is in Gauss. The only place this is the case.
    Doub TLamor = 1.0/fLamor;
    Doub dt = TLamor / NPERORBIT;
    std::cerr << "dt: " << dt;
    std::cerr <<  " iter: " << maxiter << std::endl;
    std::default_random_engine generator(int(time(NULL)));
    
    for (int ipart = 0; ipart < nparts; ipart++){

        // Vector veli = gaussian(0, vbar, generator);
//        Vector veli = flux(0, vbar, generator);
        Vector veli(69205.0, 169518.0, 0);
        Vector posi(0, 0, 0);
        Vector BNow = (*geo_).getB(posi);
        Vector vPara = veli.parallel(BNow);
        Vector vPerp = veli.perp(BNow);  
        Vector vGC = Vector(vPara.mod(), vPerp.mod(), 0);
        Particle part(posi, veli, spec);
        pushSingle(part, dt, maxiter, write, coord);
    }
    return;
}

template <class T>
double Pusher<T>::losscone(double energy, bool spec, int nparts, int maxiter, bool write, std::string& suffix)
{
    std::list<Vector> initVel;
    std::list<Vector> finlVel;
    std::list<Doub> paraVel;

    Doub mass = MI * (1 - spec) + ME * spec; // logical statement, choosing between ion and electron mass.
    Doub vbar = sqrt(energy  *  EVTOJOULE / mass); // thermal velocity, <v^2> in distribution, sigma.
    
    // decoupling dt from magnetic field
    Doub Btypical = BMAGAXIS; // magnetic field on the order of BMAGAXIS.
    Doub fLamor = ( 1520 * (1 - spec) + 2.8E6 * spec ) * Btypical; // another logical, constants from NRL p28
    Doub TLamor = 1/fLamor;
    Doub dt = TLamor / NPERORBIT;

    std::default_random_engine generator(int(time(NULL)));
//    omp_set_num_threads(1);
    #pragma omp parallel
    {
        generator.seed( int(time(NULL)) ^ omp_get_thread_num() ); // seed the distribution generator

        Vector veli, posNow, BNow, ENow, vPara, vPerp, vGC;
        Vector posi(0, 0, 0);
        Particle part(posi, veli, spec);
        std::list<Vector> initVel_private; // A list of (vpara, vperp)
        std::list<Vector> finlVel_private;
        std::list<Doub>   paraVel_private; // A list of parallel velocity at exit

        #pragma omp for private(part, veli, \
            posi, posNow, BNow, ENow, vPara, vPerp, vGC) 
            for (int i=0; i < nparts; ++i ){
    //            veli = gaussian(0, vbar, generator);
          //      veli = Vector(1000, 0, 0);
                veli = flux(0, vbar, generator);
                part.setPos(posi);
                part.setVel(veli);
                part.setSpec(spec);
                BNow = (*geo_).getB(posi);

                vPara = veli.parallel(BNow);
                vPerp = veli.perp(BNow);
                vGC = Vector(vPara.mod(), vPerp.mod(), 0); // Guiding Center velocity in (vpara, vperp)

                initVel_private.push_back(vGC);
                for (int step = 0; step < maxiter; ++step){ 
                    posNow = part.pos();
                    BNow = (*geo_).getB(posNow);
                    ENow = (*geo_).getE(posNow);

                    part.move(ENow, BNow, dt);
                    if ((*geo_).isLimiter(posNow)){
        //              std::cerr << "particle lost to limiter after" << step \
                        << "iterations." << std::endl;
                        part.lost();
                        finlVel_private.push_back(vGC); // a list of initial velocities that are lost
                        // collect the final paralell velocity here
                        vPara = part.vel().parallel(BNow);
                        paraVel_private.push_back(vPara.mod());
                        break;
                    }
                }
            }
        #pragma omp critical
            // collect everything from all threads back into the main structure
            initVel.insert(initVel.end(), initVel_private.begin(), initVel_private.end());
            finlVel.insert(finlVel.end(), finlVel_private.begin(), finlVel_private.end());
            paraVel.insert(paraVel.end(), paraVel_private.begin(), paraVel_private.end());
    }
    std::cerr << "initial velocities: " << initVel.size() << std::endl;    
    std::cerr << std::endl << "ones that were lost: " << finlVel.size() << std::endl;
    
    //TODO: FIX THIS!!! MOVE TO MAIN 

    // std::string suffix = "test.out";
    if (true){
        ofstream initial;
        initial.open("output_new/initial" + suffix);
        ofstream final;
        final.open("output_new/final" + suffix);
        while (!initVel.empty()){
            initial << initVel.front() << std::endl;
            initVel.pop_front();
        }    
        while (!finlVel.empty()){
            final << finlVel.front() << std::endl;
            finlVel.pop_front();
        }
    } 
    Doub gammaOut = 0;
    while(!paraVel.empty()){
        gammaOut += paraVel.front();
        paraVel.pop_front();
    }    
    return gammaOut;
}

template <class T>
Vector Pusher<T>::gaussian(double center, double vbar, std::default_random_engine& generator)
{
    std::normal_distribution<double> distribution(center, vbar); // generate a Gaussian distributed velocity
    double vx = distribution(generator); // generate 3 normal distributed velocities.
    double vy =  distribution(generator);
    double vz = distribution(generator);

    Vector vel(vx, vy, vz);
    return vel;
}

template <class T>
Vector Pusher<T>::flux(double center, double vbar, std::default_random_engine& generator)
{
    std::normal_distribution<double> gauss(center, vbar); // generate a Gaussian distributed velocity
    double vy =  gauss(generator);
    double vz = gauss(generator);

    std::bernoulli_distribution coinflip(0.5);
    bool coin = coinflip(generator); // choose which direction x is going
    int sign;
    if (coin){
        sign = 1;
    } else { sign = - 1; }
    
    std::uniform_real_distribution<double> uniform(0, 1); 
    double y = uniform(generator);
    double vx = sign * sqrt(- 2 * pow(vbar, 2) * log( 1 - y)); // from CDF transform.
	//std::cout << vx << std::endl;
    Vector vel(vx, vy, vz);
    return vel;
}

#endif // PUSHER_H_INCLUDED
