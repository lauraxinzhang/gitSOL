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


