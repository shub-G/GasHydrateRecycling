/*
 * driver.hh
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

#ifndef OPERATORS_DRIVER_HH_
#define OPERATORS_DRIVER_HH_

template<typename GV,typename PTree>
void driver( const GV& gv, 				// GridView
		 	 const PTree& ptree, 		// Input Parameters-Tree
			 std::string output_path,
			 Dune::MPIHelper& helper){

	// CHOOSE DOMAIN AND RANGE FIELD TYPE
	using Coord = typename GV::Grid::ctype;
	const int dim = GV::dimension;

	double time = 0.0;
	double dt = 0.0;

	// MATERIAL PROPERTIES, NUMERICAL AND TEST PARAMETERS, CONSTITUTIVE RELATIONSHIPS
	using Parameters = PropertiesAndParameters<GV,PTree>;
	Parameters parameter(gv,ptree,&time,&dt);

	std::string pathName = output_path;
	auto pathExt = ptree.get("output.path_name",(std::string)"test0");
	pathName += pathExt;
	pathName += "/";
	auto fileName = ptree.get("output.file_name",(std::string)"test");

	Operations operation;
	auto Skin = parameter.matrix.S_kin();
	auto Se   = parameter.matrix.S_eqb();
	auto Sep  = operation.S_eqb_primary();
	auto Ses  = operation.S_eqb_secondary();
	auto Seqb_secondary_inverse  = operation.S_eqb_secondary_inverse();
	auto Seqb_primary_star = operation.S_eqb_primary_star();
	auto component_matrix  = operation.ComponentMatrix();

	std::string name = "Sk";
	operation.printMatrix(name,Skin);
	name = "Se";
	operation.printMatrix(name,Se);
	name = "Sep";
	operation.printMatrix(name,Sep);
	name = "Ses";
	operation.printMatrix(name,Ses);
	name = "Ses_inv";
	operation.printMatrix(name,Seqb_secondary_inverse);
	name = "Sep_star";
	operation.printMatrix(name,Seqb_primary_star);
	name = "U";
	operation.printMatrix(name,component_matrix);

	operation.checkComponentMatrix();

	//print matrices
	std::string file_name = pathName+fileName ;
	std::string fn = file_name;
	fn += "_matrices.txt";
	std::fstream outFile;
	outFile.open(fn, std::fstream::out | std::fstream::trunc);
	if (!outFile.is_open()){
		std::cout << "Could not open file:" << fn << std::endl;
		exit(0);
	}
	outFile	<< "S_kin:" << '\n';
	for( int i=0;i<Skin.size();i++ ){
		for( int j=0;j<Skin.at(0).size();j++ ){
			outFile	<< Skin[i][j]  << " ";
		}
		outFile	<< '\n';
	}
	outFile	<< '\n';

	if(parameter.tag.Ne>0){
		outFile	<< "Seqb_secondary_inverse:" << '\n';
		for( int i=0;i<Seqb_secondary_inverse.size();i++ ){
			for( int j=0;j<Seqb_secondary_inverse.at(0).size();j++ ){
				outFile	<< Seqb_secondary_inverse[i][j]  << " ";
			}
			outFile	<< '\n';
		}
		outFile	<< '\n';

		outFile	<< "Seqb_primary_star:" << '\n';
		for( int i=0;i<Seqb_primary_star.size();i++ ){
			for( int j=0;j<Seqb_primary_star.at(0).size();j++ ){
				outFile	<< Seqb_primary_star[i][j]  << " ";
			}
			outFile	<< '\n';
		}
		outFile	<< '\n';

		outFile	<< "component_matrix:" << '\n';
		for( int i=0;i<component_matrix.size();i++ ){
			for( int j=0;j<component_matrix.at(0).size();j++ ){
				outFile	<< component_matrix[i][j]  << " ";
			}
			outFile	<< '\n';
		}
	}
	outFile.close();

	parameter.tag.check();

	/************************************************************************************************/
	dt  = ptree.get("time.dt_initial",(double)1.0); //annum
	// total simulation-time
	double t_END_0  = ptree.get("time.time_end_0",(double)1.e4); //annum
	double t_END_1  = ptree.get("time.time_end_1",(double)1.e4); //annum
	double t_END  	= t_END_1; //annum
	// output time interval
	double t_OP_0 = ptree.get("output.time_interval_0",(double)100.0); //annum
	double t_OP_1 = ptree.get("output.time_interval_1",(double)100.0); //annum
	double t_OP = t_OP_0;
	//adaptive time control
	bool is_adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
	double dt_min = ptree.get("adaptive_time_control.dt_min",(double)1.e-6); //annum
	double dt_max_0 = ptree.get("adaptive_time_control.dt_max_0",(double)2.); //annum
	double dt_max_1 = ptree.get("adaptive_time_control.dt_max_1",(double)2.); //annum
	double dt_max = dt_max_0; //annum
	int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)6);
	int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)4);

	double dtstart = dt;
	double time_op = time;
	double clock_time_elapsed = 0.;

	// prepare vtk writer
	bool plot_paraview = ptree.get("output.plot_paraview",(bool)false);
	int subsampling = 1;
	using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
	VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
	VTKSEQUENCEWRITER vtkSequenceWriter( gv,fileName,pathName,"",Dune::VTK::nonconforming);

	/************************************************************************************************/
	// MAIN
	/************************************************************************************************/
	// FEM FOR PRIMARY VARIABLES
	typedef typename GV::Grid::ctype Coord;
	using CON0 = Dune::PDELab::P0ParallelConstraints;
	using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;
	using FEMP0 = Dune::PDELab::QkDGLocalFiniteElementMap<Coord,double,0,dim,Dune::PDELab::QkDGBasisPolynomial::lagrange>;
	FEMP0 femp0;

	//1. GFS
	using GFS0 = Dune::PDELab::GridFunctionSpace<GV, FEMP0, CON0, VBE0>;
	GFS0 gfs0(gv, femp0);
	using GFS = Dune::PDELab::PowerGridFunctionSpace< GFS0,
													  Tags::nEqns,
													  Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>,
													  Dune::PDELab::EntityBlockedOrderingTag >;
	GFS gfs(gfs0);

	using CC = typename GFS::template ConstraintsContainer<double>::Type;
	CC cc;
	cc.clear();

	//2. VECTOR CONTAINER FOR THE SOLUTION
	using U = Dune::PDELab::Backend::Vector<GFS,double>;
	U u_old(gfs,0.0);
	U u_new(gfs,0.0);

	//3. INITIAL CONDITIONS
	InitialConditions<GV,PTree,Parameters> ic(gv,ptree,parameter);
	InitialValues<GV,InitialConditions<GV,PTree,Parameters>,GFS,U> icvalues( gv,ic,gfs,&u_old );
	icvalues.evaluate();
	//initialize u_new
	u_new = u_old;

	//4. BOUNDARY CONDITIONS
	using BC = BoundaryConditions<GV,PTree,Parameters>;
	BC bc(gv,ptree,parameter) ;

	//5. LOCAL OPERATORS
	using LOP = LocalOperator<GV,Parameters,BC>;
	LOP lop( gv, parameter, bc,
			 &Skin, &Seqb_secondary_inverse, &Seqb_primary_star, &component_matrix,
			 &time, &dt );
	using TOP = TimeOperator<GV,Parameters>;
	TOP top( gv, parameter,
			 &Seqb_secondary_inverse, &Seqb_primary_star, &component_matrix,
			 &time, &dt );
//
	//6. GRID OPERATORS
	using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
	MBE mbe(35);
	using GOLOP = Dune::PDELab::GridOperator< GFS, GFS, LOP, MBE, double, double, double, CC, CC>;
	GOLOP golop(gfs, cc, gfs, cc, lop, mbe);
	using GOTOP = Dune::PDELab::GridOperator< GFS, GFS, TOP, MBE, double, double, double, CC, CC>;
	GOTOP gotop(gfs, cc, gfs, cc, top, mbe);

	using IGO = Dune::PDELab::OneStepGridOperator< GOLOP, GOTOP >;
	IGO igo( golop, gotop );

	//7. LINEAR SOLVER
//	using LS = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO>;
//	LS ls(gfs,100,1,true,true);
	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
	LS ls;

	//8. NONLINEAR SOLVER
//	using PDESOLVER = Dune::PDELab::NewtonMethod<IGO,LS>;
//	PDESOLVER newton(igo,ls);
//	newton.setParameters(ptree.sub("newton"));
	typedef Dune::PDELab::Newton< IGO, LS, U > PDESOLVER;
	PDESOLVER newton( igo, ls );
	newton.setLineSearchStrategy(PDESOLVER::Strategy::noLineSearch);
	newton.setReassembleThreshold(0.0);
	newton.setVerbosityLevel(ptree.get("newton.verbosity_0",(int)2));
	newton.setReduction(ptree.get("newton.reduction_0",(double)1e-6));
	newton.setMinLinearReduction(ptree.get("newton.min_linear_reduction_0",(double)1e-6));
	newton.setMaxIterations(ptree.get("newton.max_iterations_0",(int)20));
	newton.setForceIteration(ptree.get("newton.force_iteration_0",(bool)true));
	newton.setAbsoluteLimit(ptree.get("newton.absolute_error_0",(double)1.e-5));

	//9. ONE STEP METHOD
	// time stepping methods
    Dune::PDELab::ImplicitEulerParameter<double> method1;
    Dune::PDELab::OneStepThetaParameter<double> method2(0.0);//Crank-Nicholson -> 0.5, Implicit Euler -> 1.0, Explicit Euler -> 0.0
    Dune::PDELab::Alexander2Parameter<double> method3;
    Dune::PDELab::Alexander3Parameter<double> method4;
    Dune::PDELab::FractionalStepParameter<double> method5;
    Dune::PDELab::HeunParameter<double> method6;
    Dune::PDELab::Shu3Parameter<double> method7;
    Dune::PDELab::RK4Parameter<double> method8;

	Dune::PDELab::TimeSteppingParameterInterface<double> *method = &method1;
	Dune::PDELab::OneStepMethod< double, IGO, PDESOLVER, U, U > osm( *method, igo, newton );
	osm.setVerbosityLevel(2);

	//10.VTK-NAMES
	if(plot_paraview){
		for(int i=0; i<parameter.tag.nEqns; i++){
			gfs.child(i).name(parameter.tag.getPrimaryVariablesName(i));
		}
		Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs,u_new);
	}

	//11.POST PROCESSES
	std::string qoi_file = pathName;
	qoi_file +=fileName;
	PostProcess<GV,Parameters,GFS,U> postprocess( qoi_file, gv, parameter,
										  	  	  gfs, &u_new,
												  &Skin, &Seqb_secondary_inverse, &Seqb_primary_star, &component_matrix,
												  &time, &dt);
	for(int i=0; i<=parameter.tag.NPP;i++)
		postprocess.print_qoi(i);

	/***********************************************/
	//	BEGIN TIME LOOP
	/***********************************************/
	int opcount = 1;
	double t_OP_next = t_OP;
	double dtLast = dtstart;
	int dtFlag = 0;
	bool isPart0 = true;
	bool exceptionCaught = false;
	unsigned int newton_iterations = 0;

	while( time < t_END - 1e-8 ){

		if( exceptionCaught==false ){
			dt = std::max(dt,dt_min);
		}

		if(helper.rank()==0){
		std::cout<< "_____________________________________________________" <<std::endl;
		std::cout<< " current opcount = " << opcount - 1 << std::endl;
		}

		clock_t start = clock();
		try{
			std::cout<< "––––––––––––––––––––––––––" << std::endl;
			std::cout<< "SOLVING:" << std::endl;
			std::cout<< "––––––––––––––––––––––––––" << std::endl;
			osm.apply( time, dt, u_old, u_new );
			newton_iterations = osm.getPDESolver().result().iterations;

			exceptionCaught = false;

		}catch ( Dune::Exception &e ) {
			exceptionCaught = true;
			if( dt > 1e-15 ){

				if(helper.rank()==0){
				std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}

				u_new = u_old;
				newton_iterations = 0;

				dt *= 0.5;
				continue;
			}
			else
			{
				if(helper.rank()==0){
					std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
				}
				exit(0);
			}
		}
		clock_t end = clock();
		double clock_time_this_step = (double) (end-start) / CLOCKS_PER_SEC;
		clock_time_elapsed += clock_time_this_step;

		if(helper.rank()==0){
		std::cout<<"DONE"<<std::endl;
		std::cout<<"_____________________________________________________"<<std::endl;
		}

		/*********************************************************************************************
		 * OUTPUT
		 *********************************************************************************************/
		/* At each time step: **Statistics**
		 * t_new,
		 * dt,
		 * fixed point iterations,
		 * newton iterations per fixed point iteration,
		 * total newton terations
		 */
		if(helper.rank()==0){
		std::string statistics_file = pathName;
		statistics_file +=fileName;
		statistics_file +="_statistics";
		statistics_file += ".txt";
		parameter.ReportStatistics( statistics_file,
									time,
									dt,
									newton_iterations,
									clock_time_elapsed );
		}

		/* At t_OP
		 *
		 */
//		if( ( time+dt > t_OP_next - dt_min ) and ( time+dt < t_OP_next + dt_min )  )
//		{
//			// EVALUATE SECONDARY VARS
//			postprocess.print_qoi();
//
//			// WRITE OUTPUT
//			if(plot_paraview){
//				vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
//			}
//
//			if(helper.rank()==0){
//			std::cout<< " ******************************************************************* " << std::endl;
//			std::cout<< " OUTPUT WRITTEN " << opcount << std::endl;
//			std::cout<< " ******************************************************************* " << std::endl;
//			std::cout<< std::flush;
//			}
//
//			if( isPart0==true and time+dt>t_END_0){
//				t_OP = t_OP_1;
//				dt_max = dt_max_1;
//				newton.setVerbosityLevel(ptree.get("newton.verbosity_1",(int)2));
//				newton.setReduction(ptree.get("newton.reduction_1",(double)1e-6));
//				newton.setMinLinearReduction(ptree.get("newton.min_linear_reduction_1",(double)1e-6));
//				newton.setMaxIterations(ptree.get("newton.max_iterations_1",(int)20));
//				newton.setForceIteration(ptree.get("newton.force_iteration_1",(bool)true));
//				newton.setAbsoluteLimit(ptree.get("newton.absolute_error_1",(double)1.e-5));
//				isPart0 = false;
//			}
//			t_OP_next = time+dt+t_OP;
//			opcount = opcount+1;
//		}
		if( ( time+dt > t_OP_next - dt_min ) and ( time+dt < t_OP_next + dt_min )  )
		{
			// if( time+dt > t_END_0 ){
				// EVALUATE SECONDARY VARS
				for(int i=0; i<=parameter.tag.NPP;i++)
					postprocess.print_qoi(i);

				// WRITE OUTPUT
				if(plot_paraview){
					vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
				}

				if(helper.rank()==0){
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< " OUTPUT WRITTEN " << opcount << std::endl;
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< std::flush;
				}
			// }

			if( isPart0==true and time+dt>t_END_0){
				t_OP = t_OP_1;
				dt_max = dt_max_1;
				newton.setVerbosityLevel(ptree.get("newton.verbosity_1",(int)2));
				newton.setReduction(ptree.get("newton.reduction_1",(double)1e-6));
				newton.setMinLinearReduction(ptree.get("newton.min_linear_reduction_1",(double)1e-6));
				newton.setMaxIterations(ptree.get("newton.max_iterations_1",(int)20));
				newton.setForceIteration(ptree.get("newton.force_iteration_1",(bool)true));
				newton.setAbsoluteLimit(ptree.get("newton.absolute_error_1",(double)1.e-5));
				isPart0 = false;
			}
			t_OP_next = time+dt+t_OP;
			opcount = opcount+1;
		}

		// PREPARE FOR NEXT TIME INTEGRATION
		//1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		u_old = u_new;

		//2. ADVANCE TIME:
		time += dt;

		if(helper.rank()==0){
		std::cout<<" "<< std::endl;
		std::cout<< " time = " << time ;
		std::cout<< std::flush;
		}

		if( is_adaptive_time_control ){
			if(newton_iterations>maxAllowableIterations){
				dt=std::max(dt*0.9,dt_min);
			}
			else if(newton_iterations<=minAllowableIterations){
				dt=std::min(dt*1.1,dt_max);
			}
		}
		else{
			dt = dtstart;
		}

		if(helper.rank()==0){
		std::cout << " , time+dt = " << (time + dt)
				  << " , opTime = "  << t_OP_next ;
		std::cout<< std::flush;
		}

		if( time + dt  > t_OP_next){
			dtLast = dt;
			dt 	 = t_OP_next - time ;

			if(helper.rank()==0){
			std::cout<< " , because timeNext > opNext , dt set to : " << dt << std::endl;
			std::cout<< std::flush;
			}

			dtFlag = 0;
		}
		dtFlag += 1;

		if(helper.rank()==0){
		std::cout<< " , dt  : " << dt << std::endl;
		std::cout<<" "<< std::endl;
		std::cout << " READY FOR NEXT ITERATION. " << std::endl;
		std::cout<< std::flush;
		}

	}

}

#endif /* OPERATORS_DRIVER_HH_ */