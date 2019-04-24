/** 
 * \file    Beam.cc
 * \author  Xin Zhang, PPPL
 * \date    March 2019
 *
 * \brief   Deprecated part of the Mirror class that was used to produce NBI simulation
 * \note    Move back into their original class to reproduce NBI simulations; DO NOT compile
 *          by itself.
 */

class Mirror
{
	public:

		/**
		 * \brief Check whether the particle is crossing a sightline, return the index of that sightline
		 * \return Index of sightline that Vector pos is on; 0 if not on any sightlines
		 */
		int sightline(Particle& part, int lastCrossed);

		/**
		 * \brief Parse input sightline definitions
		 *
		 * \param inputSL a string that includes the to the input file
		 */
		void setSightlines(const std::string& inputSL, int rows);

	private:
		std::vector<double> sightlines_; //default constructed
};



int Mirror::sightline(Particle& part, int lastCrossed)
{
	std::string prefix  = "./SLoutput/sl";
	std::string suffix  = ".out";
	ofstream output;

	if (sightlines_.size() == 0){
		std::string path = "./input/viewlineBL.csv";
		setSightlines(path, 26);
	}

	if ( abs( part.pos().z() ) < 0.008 ){ // pos is in the plane of the sightlines
		// std::cerr << "in the right plane" << std::endl;
		for (int i = lastCrossed + 1; i < 26; i++){
		// i is index for rows, goes from 0 to 25
			double xinter = sightlines_[3 * i + 1];
			double slope  = sightlines_[3 * i + 2];

			double yexpect = slope * ( part.pos().x() - xinter );
			//std::cerr << "i = " << i << ", xinter" <<xinter <<  ", slope" << slope << ", x: " << part.pos().x()<< std::endl;
			double deltaY = abs( yexpect - part.pos().y() );
			//std::cerr << ", deltaY = " << deltaY << ", yexpect: " << yexpect << ", pos.y " << part.pos().y() << std::endl;
			if ( deltaY < 0.008 ){
				// pos is on sightline numbered i+1
				//std::cerr << "crossing sightline #" << i+1 << std::endl;

				Vector sl(1, slope, 0); // define vector direction for current sightline
				Vector vpara = part.vel().parallel(sl); // find parallel velocity

				double parallel = vpara.dot(sl.normalize());
				output.open(prefix + std::to_string(i+1) + suffix, ios::app);
				output << parallel << std::endl;
				//output.close();
				//output.clear();
				return i;
			}
		}
	}
	return lastCrossed;
}

void Mirror::setSightlines(const std::string& inputSL, int rows)
{
	ifstream input;
	input.open(inputSL);
	if (!input.is_open()) {
    	std::cerr << "Unable to open sightline file" << std::endl; 
    }

    double val;
    for (int row = 1; row <= 3 * rows; row++){
    	input >> val;
    //	std::cerr << "reading val = "<< val << std::endl;
	sightlines_.push_back(val);
    }
    //for (int i = 0; i < 26; i++){
//	std::cerr << sightlines_[3*i + 1] << "," << sightlines_[3 * i + 2] << std::endl;
  //  }
    return;
}

//---------------------------------------  NBI   ---------------------------------------------


template <class T>
void Pusher<T>::gridBurst(double radius, double ylim, int nsources, bool write)
{
	std::ofstream coord;
	coord.open("coordBurst.out");
	coord << std::setprecision(10);

	std::default_random_engine generator(int(time(NULL)));
	// std::uniform_real_distribution<double> distribution(-1 * ylim, ylim);

	for (int isource = 0; isource < nsources; isource++){
		// Vector posi(xCalc, yRand, zRand);
		// Vector veli((-1*xCalc + radius), -1* yRand, -1*zRand);
		Vector posi = sphere(radius, ylim, generator);
		Vector veli = sphereNormal(radius, posi).normalize();

		//std::cerr << posi << std::endl;
		//std::cerr << veli << std::endl;
		Particle part(posi, veli, 1, 0); // one particle per source for now
		pushSingle(part, 0.001, 4000, write, coord);

	}

	return;
}

template <class T>
void Pusher<T>::conicBurst(double radius, double ylim, double dtheta, int nsources, int partPerS, bool write)
{
	std::ofstream conic;
	conic.open("conicBurst.out");
	conic << std::setprecision(10);

	std::default_random_engine generator(int(time(NULL))); // initialize outside loop to avoid overseeding
	
	// std::uniform_real_distribution<double> location(-1 * ylim, ylim); // for particle sources
	// std::normal_distribution<double> pitchAngle(0, dtheta);

	for (int iS = 0; iS < nsources; iS++){
		Vector posi = sphere(radius, ylim, generator);
	//	std::cerr << "source #" << iS << std::endl;
		for (int n = 0; n < partPerS; n++){
			Vector veli = diverge(radius, posi, dtheta, generator);
			//std::cerr << posi << std::endl;
			//std::cerr << veli << std::endl;
			
			Particle part(posi, veli, 1, 0);
			pushSingle(part, 0.001, 3000, write, conic);
		}
	}
	return;
}

template <class T>
Vector Pusher<T>::sphere(double radius, double ylim, std::default_random_engine& generator)
{
	std::uniform_real_distribution<double> distribution(-1 * ylim, ylim);

	double yRand = distribution(generator);	
	double zRand = distribution(generator);
	while (yRand * yRand + zRand * zRand >= ylim * ylim){
		// if it's not in the circle, try again.
		yRand = distribution(generator);
        zRand = distribution(generator);
    }
	// double xCalc = -1 * sqrt(radius * radius - yRand * yRand - zRand * zRand) + radius;
    double xCalc = sqrt(radius * radius - yRand * yRand - zRand * zRand) - radius; // flip beam source to the right

	Vector posi(xCalc, yRand, zRand);
	return posi;
}

template <class T>
Vector Pusher<T>::sphereNormal(double radius, Vector pos)
{
	// double x = radius - pos.x();
	double x = 0 - radius - pos.x(); // flip beam source to the right
	double y = -1 * pos.y();
	double z = -1 * pos.z();
	Vector result(x, y, z);
	return result;
}

template <class T>
Vector Pusher<T>::diverge(double radius, Vector& pos, double dtheta, std::default_random_engine& generator)
{
	std::normal_distribution<double> pitchAngle(0, dtheta);
	std::uniform_real_distribution<double> uni(-1, 1);

	Vector posNorm = pos.normalize();
	double x0 = posNorm.x();
	double y0 = posNorm.y();
	double z0 = posNorm.z();

	double b = uni(generator);
	double c = uni(generator);
	// double a = (radius * x0 - b * y0 - c * z0) / (x0 - radius);
	double a = (0 - radius * x0 - b * y0 - c * z0) / (x0 + radius); // flipping beam source


	Vector tangent(a - x0, b - y0, c - z0); // a randomly generated vector tangent to sphere at pos
	Vector velNorm = sphereNormal(radius, pos).normalize(); // normal vector of sphere at pos

	double theta = pitchAngle(generator);
	Vector vperp = tangent.normalize() * tan(theta);

	Vector result = velNorm + vperp;

	return result.normalize();
}



