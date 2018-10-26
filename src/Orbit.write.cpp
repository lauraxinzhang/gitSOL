/**
 * \file    Orbit.write.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Implementations for outputting utilities of the Orbit class.
 *
 * \NOTE    See header file for detailed explainations for each 
 *          function \n
 *
 *          Uses Numerical Recipe nr3.h data containers for Vector and Matrix
 *          functionalities; uses Vector.h for 3D vectors with arithmetic operations.
 *      
 *          Uses nr3 routines for interpolations.     
 *
 *          Reads data from g-eqdsk files.
 *      
 *          Uses double instead of float to avoid floating number errors.
 */


#include "Orbit.h"

void Orbit::fieldLines(std::ofstream &rBList, std::ofstream &zBList, std::ofstream &phiBList, const Vector& init, Doub dl, int iter)
{
    if (!rBList.is_open()) {
    	std::cerr << "Unable to open file" << std::endl; 
    }

    Vector now = init;

	for (int i = 0; i < iter; ++i){
		Doub rNow = now.x();
		Doub phiNow = now.y();
		Doub zNow = now.z();

		if (rNow >= rleft_ + rdim_ || zNow >= zmid_ + (zdim_/2) || zNow <= zmid_ - (zdim_/2)){
			break;
		}
		rBList   << std::setprecision(10) << rNow   << std::endl;
    	phiBList << std::setprecision(10) << phiNow << std::endl;
    	zBList   << std::setprecision(10) << zNow   << std::endl;

    	Vector bNow = getB(now);
		Vector bHat = bNow.normalize();

		Vector dB = bHat * dl;
		now = now + dB;
	}

	return;
}

void Orbit::eField(Doub Ti, Doub Te, Vector init, Doub dl, int iter, std::ofstream &output)
{
	if (!output.is_open()) {
    	std::cerr << "Unable to open output file" << std::endl; 
    }
	if (Rratio_ == nullptr){
		configMirror();
	}

	// Interpolate mirror ratio onto entire RZ plane
	INTERP2D mirrorRZ((*rGrid_), (*zGrid_), (*Rratio_));

    Vector  now  = init;
    Doub    lTot = 0;

    VecDoub lList(iter);
    VecDoub potList(iter);

    Doub maxIter(0);

	for (int i = 0; i < iter; ++i){
		Doub rr = now.x();
		Doub zz = now.z();

		std::cout << rr << "," << zz << std::endl;

		if (isLimiter(rr, zz)){
			std::cerr << "limiter reached after " << i << "iterations " << std::endl;
			break;
		} else {

			Doub ratioNow  = mirrorRZ.interp(rr, zz);

			if (std::isnan(ratioNow)){
				std::cerr << "mirror ratio undefined before limiter is reached. Break after " << i << "iterations " << std::endl;
				break;
			}
			Doub potential = pastukhov(Ti, Te, ratioNow);

			Vector bNow = getB(now);
			Vector bHat = bNow.normalize(); // find the direction of the magnetic field

			Vector dBl = bHat * dl;
			now = now + dBl;  // Each step's arc length is exactly dl

			lTot +=dl; // Add to arclength

			lList[i]   = lTot;
			potList[i] = potential;
			// std::cerr << lList[i] << " , " << potList[i]<< " , " << ratioNow << std::endl;
		}
		maxIter++;
	}

	// Done with calculating potential

	if (maxIter < iter){
	// move the elements over to a smaller container to not mess up the interpolation

		VecDoub lListNew(maxIter);
		VecDoub potListNew(maxIter);

		for (int i = 0; i < maxIter; ++i){
			lListNew[i] = lList[i];
			potListNew[i] = potList[i];
		}
		lList = lListNew;
		potList = potListNew;
	}
	// I don't know what this is doing anymore. Why are there trailing zeros in front?	
	int zeroHead(0);
	while(potList[zeroHead] == 0){
		++zeroHead;
	}
	potList.erase(potList.begin(), potList.begin() + zeroHead);
	

    VecDoub eList(potList.size());


	// Now calculate electric field via numerical differentiation.

	eFieldHelp phiOfL(lList, potList);
	Doub err;

	for (int i = 0; i < potList.size(); ++i){

		Doub lc = lList[i]; // get current arc location
		if (potList[i] == 0){
			eList[i] = NAN;
		} else {
	    	eList[i] = dfridr(phiOfL, lc, 0.2, err);
	    	// << eList[i] << ", " << err << std::endl;
	    	// std::cerr << lc << ", " << potList[i] << ", "<< phiOfL(lc) << ", " << eList[i] << ", " << err << std::endl;
	    	// And then write results to file.
	    }
    	output << std::setprecision(10) << lc << "," << potList[i] << "," << eList[i] << std::endl;
    }
	return;
}

void Orbit::writeModB(std::ofstream &output)
{
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*Bmod_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writeFlux(std::ofstream &output)
{
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*psiRZ_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writeMirrorRatio(std::ofstream &output)
{
	if (Rratio_ == nullptr){
		configMirror();
	}
	// MatDoub ratio(nw_, nh_);
	// mirrorRatio(ratio);
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*Rratio_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writePastukhov( Doub Ti, Doub Te, std::ofstream &output)
{
	std::cerr << "calling setPastukhov" << std::endl;
	setPastukhov(Ti, Te);
	std::cerr << "writing to file." << std::endl;

	if (!output.is_open()) {
    	std::cerr << "Unable to open output file" << std::endl; 
    	exit(0);
    }

	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
	    	output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*Phi_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::midPlanePhi( Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output)
{
	INTERP2D pasInterp((*rGrid_), (*zGrid_), pas);

	for(Doub r = rllmtr_; r < rrlmtr_; r += 0.002){
		Doub x = pasInterp.interp(r, 0);
		if (!std::isnan(x)) {
			output << r << ',' << x << std::endl;
		}
	}
	return;
}

void Orbit::temperature(Doub Ti_start, Doub dT, int iter, Doub R)
{
	ofstream temp;
	temp.open("./output/phivsT.out");
	temp << std::setprecision(10);

	int i(0);
	Doub TiNow = Ti_start;
	while (i <= iter){
		Doub phi   = pastukhov(TiNow, 1, R);
		temp << TiNow << ',' << phi << std::endl;
		i++;
		TiNow += dT;
	}
	temp.close();

	return;
}

void Orbit::printData()
{
	/*
	ofstream modB;
	modB.open("./output/modB.out");

	orbit.writeModB(modB);

	modB.close();

	ofstream psiout;
	psiout.open("./output/flux.out");
	orbit.writeFlux(psiout);
	psiout.close();

	ofstream RRout;
	RRout.open("./output/mirrorRatio.out");
	orbit.writeMirrorRatio(RRout);
	RRout.close();

	/*
	MatDoub pass(orbit.nw_, orbit.nh_);

	*/
	ofstream passout1; // LOLOLOLOL
	passout1.open("./output/passTest.out");
	// input Ti, Te here. keep Te 1.
	writePastukhov(0.2, 1, passout1); // pass modified
	passout1.close();

	/*
	ofstream midplane1;
	midplane1.open("./output/midplane02.out");
	orbit.midPlanePhi(0.2, 1, pass, midplane1);
	midplane1.close();

	ofstream eFieldLine2;
	eFieldLine2.open("./output/eField2.out");

	Doub rStart = orbit.rllmtr_ + 0.2;
	Doub zStart = abs(orbit.getzLimiter(rStart)) - 0.01;
	// Vector start( 0.575, 0, 0);
	Vector start(rStart, 0, zStart);
	orbit.eField(2, 1, start, 0.005, 15000, eFieldLine2);
	*/
}