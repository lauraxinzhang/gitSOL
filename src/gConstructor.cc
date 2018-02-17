Orbit::Orbit(const std::string& path, int flag)
{
	ifstream myfile;
    myfile.open(path); // open eqdsk file

    if (!myfile.is_open()) {
    	std::cerr << "Unable to open file" << std::endl; 
    } else {
    	std::cerr << "Eqdsk successfully loaded" << std::endl;
    }

    std::string a;

    while (myfile >> a){
    	// truncates all the crap in front until a number is reached. Assuming
    	// 
    	if (std::isdigit(a[0])){
    		int dum = std::stoi(a);
    		if (dum == 0){ break; }
    	}	    
	}

	// std::cout << 1;
	// std::cout << std::endl;

	// while (myfile >> a){
	// 	if (std::isdigit(a[0]) || a[0]=='-'){
	// 	    	double d = std::stod (a);
	// 	    	std::cout << std::setprecision(10) << d << ',' << std::endl;
	// 	} 
	// }

	// line 1:
	// int nw_, nh_;
	myfile >> std::setprecision(10) >> nw_ >> nh_;

	// line 2:
	double rdim, zdim, rcentr, rleft, zmid;
	myfile >> std::setprecision(10) >> rdim >> zdim >> rcentr >> rleft >> zmid;
	// std::cout << rdim << zdim << rcentr << rleft << zmid << std::endl;
	rdim_ = rdim;
	zdim_ = zdim;
	rleft_ = rleft;
	zmid_ = zmid;

	// std::cout << std::setprecision(10) << nx << ' ' << nz << std::endl;
	// std::cout << std::setprecision(10) << xdim <<' ' << zdim << ' ' << rcentr << ' ' << rleft << ' ' << zmid << std::endl;

	// std::cout << 2;
	// std::cout << std::endl;

	// line 3:
	double rmaxis, zmaxis, simag, sibry, bcentr;
	myfile >> std::setprecision(10) >> rmaxis >> zmaxis >> simag >> sibry >> bcentr;
	simag_ = simag;
	sibry_ = sibry; 

	// line 4:
	double current, simag1, xdum, rmaxis1, xdum1;
	myfile >> std::setprecision(10) >> current >> simag1 >> xdum >> rmaxis1 >> xdum1;

	// line 5:
	double zmaxis1, xdum2, sibry1, xdum3, xdum4;
	myfile >> std::setprecision(10) >> zmaxis1 >> xdum2 >> sibry1 >> xdum3 >> xdum4;

	// std::cout << (rmaxis == rmaxis1 && zmaxis == zmaxis1 && sibry == sibry1) << std::endl;

    // while (myfile >> a)
    // {
    //     printf("%f ", a);
    // }

 //    std::cout << 3;
	// std::cout << std::endl;

	VecDoub fpol(nw_);
	VecDoub * siGrid = new VecDoub(nw_);
	dsi_ = (sibry_ - simag_)/nw_;

    for (int i = 0; i < nw_; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	fpol[i] = d;
    	(*siGrid)[i] = simag_ + i * dsi_;
    }
    siGrid_ = siGrid;

    VecDoub pres(nw_);
    for (int i = 0; i < nw_; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	pres[i] = d;
    }

    VecDoub ffprim(nw_);
    for (int i = 0; i < nw_; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	ffprim[i] = d;
    }


    VecDoub pprime(nw_);
    for (int i = 0; i < nw_; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	pprime[i] = d;
    }

 //    std::cout << 4;
	// std::cout << std::endl;

    MatDoub psirz(nw_, nh_);

    VecDoub * rGrid = new VecDoub(nw_);
    VecDoub * zGrid = new VecDoub(nh_);

    for (int i = 0; i < nw_; ++i){     // 'column' index, in r
    	for (int j = 0; j < nh_; ++j){ // 'row' index, in h
    		double d;
    		myfile >> std::setprecision(10) >> d;
    		psirz[i][j] = d;
    	}
    }

    for (int i = 0; i < nw_; ++i){     // 'column' index
    	(*rGrid)[i] = rleft + i * (rdim_/nw_);
    }

    for (int j = 0; j < nh_; ++j){ // 'row' index
		(*zGrid)[j] = (zmid - zdim_/2.0) + j * (zdim_/nh_);
	}

    rGrid_ = rGrid;
    zGrid_ = zGrid;

 //    std::cout << 5;
	// std::cout << std::endl;

    VecDoub qpsi(nw_);
    for (int i = 0; i < nw_; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	qpsi[i] = d;
    }    

    int nbbbs, limitr; // number of boundary points and limiter points
    myfile >> std::setprecision(10) >> nbbbs >> limitr;

    VecDoub rbbbs(nbbbs), zbbbs(nbbbs);
    for (int i = 0; i < nbbbs; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	rbbbs[i] = d;
    } 

    for (int i = 0; i < nbbbs; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	zbbbs[i] = d;
    } 

    /**
    // VecDoub rlim(limitr), zlim(limitr);
    VecDoub * rlim = new VecDoub(limitr);
    VecDoub * zlim = new VecDoub(limitr);

    for (int i = 0; i < limitr; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	(*rlim)[i] = d;
    }

    for (int i = 0; i < limitr; ++i){
    	double d;
    	myfile >> std::setprecision(10) >> d;
    	(*zlim)[i] = d;
    } 

    rlim_ = rlim;
    zlim_ = zlim; 
	*/

    // fpol_ = fpol;
    // psirz_ = psirz;
    // rlim_ = rlim;
    // zlim_ = zlim;

    setupB(fpol, psirz); // we don't need to store fpol and psirz, just calculate magnetic field
 //    std::cout << 6;
	// std::cout << std::endl;

  return;
}



void Orbit::setupB( VecDoub& fpol, MatDoub& psirz)
{
    int nr = rGrid_ -> size();
    int nz = zGrid_ -> size();

    MatDoub * Br   = new MatDoub(nr, nz);
    MatDoub * Bz   = new MatDoub(nr, nz);
    MatDoub * Btor = new MatDoub(nr, nz);

    // Setup interpolations for Btor:
    Linear_interp fpolPsi( (*siGrid_), fpol );
    Bilin_interp psiRZ( (*rGrid_), (*zGrid_), psirz) ;


    for (int ir = 0; ir < nr; ir++ ){
        // iterate over each r
        for (int iz = 0; iz < nz; iz++){
            // iterate over each z

            Doub rr = (*rGrid_)[ir];
            Doub zz = (*zGrid_)[iz];

            // Setup Bpol
                // construct partial functions
            Psir psir( (*rGrid_), (*zGrid_), psirz, zz ); 
            Psiz psiz( (*rGrid_), (*zGrid_), psirz, rr );

            Doub err;
            const Doub h( rdim_/(nr*10.0) ); // differentiation stepsize = x_c/10;

            (*Br)[ir][iz] = dfridr(psiz, zz, h, err) / (-2.0 * PI * rr); // differentiate psi(z) at z = zz
            (*Bz)[ir][iz] = dfridr(psir, rr, h, err) / ( 2.0 * PI * rr); // differentiate psi(r) at r = rr

            // Setup Btor
            Doub psi = psiRZ.interp(rr, zz);
            Doub frz = fpolPsi.interp(psi); // found f(r,z)

            (*Btor)[ir][iz] = frz/rr;

        }
    }
    Br_   = Br;
    Bz_   = Bz;
    Btor_ = Btor;

    return;
}