/*
 * main.cc
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

/*********************************************************************/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<exception>
#include<chrono>
/*********************************************************************/
#define DIMENSION 2
/*********************************************************************/
#include"../duneincludes.hh"
#include"parameters.hh"
#include"initial.hh"
#include"boundary.hh"
#include"../operators/operations.hh"
#include"../operators/initialvaluefunction.hh"
#include"../operators/boundaryvaluefunction.hh"
#include"../operators/localoperator.hh"
#include"../operators/timeoperator.hh"
#include"postprocess.hh"
#include"../operators/driver.hh"
/*********************************************************************/

int main(int argc, char** argv)
{
	try{
		//initialize MPI
		Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
	    if(helper.rank()==0){
	    std::cout << "This is project GasHydrateRecycling." << std::endl;
	    }
	    if(Dune::MPIHelper::isFake){
	      std::cout<< "This is a sequential program." << std::endl;
	    }
	    else {
	    	if(helper.rank()==0){
	    		std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()<<" processes!"<<std::endl;
	    	}
	    }
		/**************************************************************************************************/
		// INPUTS
	    if (argc!=3)
	    {
	    	if(helper.rank()==0){
	    		std::cout << "usage: ./problem <user-name> <input-file> " << std::endl;
	    	}
	        return 1;
	    }
        /**************************************************************************************************/
		// USER-NAME
	    char USER_NAME[100];
	    sscanf(argv[1],"%99s", USER_NAME);
	    // DUNE MODEL PATH
	    std::string MODEL_PATH = "/home/";
		MODEL_PATH += USER_NAME;
		MODEL_PATH += "/dune_2_8/GasHydrateRecycling/src/";
	    // PROBLEM NAME
	    std::string PROBLEM_NAME = "problem/";
	    // INPUT PATH NAME
	    std::string INPUT_PATH  = MODEL_PATH + PROBLEM_NAME + "inputs/";
	    // OUTPUT PATH NAME
	    std::string OUTPUT_PATH = MODEL_PATH + "outputs/" + PROBLEM_NAME;
	    // INI-FILE FOR USER-DEFINED INPUTS
	    char INI_FILE[100];
	    sscanf(argv[2],"%99s", INI_FILE);
	    std::string input_file = INPUT_PATH;
	    input_file += INI_FILE;
	    input_file += ".ini";
        if(helper.rank()==0){
	    std::cout<< "input file: " << input_file << std::endl ;
        }
        /**************************************************************************************************/
        // PARAMETER TREE
	    Dune::ParameterTree ptree;
	    Dune::ParameterTreeParser ptreeparser;
	    ptreeparser.readINITree(input_file,ptree);
	    ptreeparser.readOptions(argc,argv,ptree);
		/**************************************************************************************************/
		// MESH
#if DIMENSION==2
		const int dim = 2;
#elif DIMENSION==3
		const int dim = 3;
#else
		std::cout<< "Invalid DIMENSION:" << DIMENSION
				 << " in LINE:" << __LINE__
				 << " of FILE:" << __FILE__
				 << std::endl;
		exit(0);
#endif
		//WITH UG MESH
#ifdef USE_UG
		using GridType = Dune::UGGrid<dim>;
		GridType grid_type;
#if DIMENSION==2
		const std::string grid_name = ptree.get("domain.ug.name",(std::string)"default_1x1");
#elif DIMENSION==3
		const std::string grid_name = ptree.get("domain.ug.name",(std::string)"default_1x1x1");
#else
		std::cout<< "Invalid DIMENSION:" << DIMENSION
				 << " in LINE:" << __LINE__
				 << " of FILE:" << __FILE__
				 << std::endl;
		exit(0);
#endif
		auto grid_file = MODEL_PATH + PROBLEM_NAME + "gmsh/" + grid_name + ".msh";
		Dune::GmshReader<GridType> gmshreader;
		std::shared_ptr<GridType> grid(gmshreader.read(grid_file,true,false));

		using GV = GridType::LeafGridView;
		GV gv = grid->leafGridView();
        grid->loadBalance();
#else
        //WITH YASP (DEFAULT)
		Dune::FieldVector<double,dim> L(1.0);	// L represents the right top node of the rectangular/cuboidal domain
		std::array<int,dim> N;

		double X = ptree.get("domain.yasp.Xmax",(double)1.0);//in m
		X *= 100.; //convert to cm
		const int X_cells = ptree.get("domain.yasp.nX",(int)1) ;
		L[0] = X;
		N[0] = X_cells ;
#if DIMENSION==2
		double Z = ptree.get("domain.yasp.Zmax",(double)400.0);//in m
		Z *= 100.;// convert to cm
		const int Z_cells = ptree.get("domain.yasp.nZ",(int)1000) ;
		L[1] = Z;
		N[1] = Z_cells ;
#elif DIMENSION==3
		double Y = ptree.get("domain.yasp.Ymax",(double)400.0);//in m
		Y *= 100.;// convert to cm
		const int Y_cells = ptree.get("domain.yasp.nY",(int)2000) ;
		double Z = ptree.get("domain.yasp.Zmax",(double)400.0);//in m
		Z *= 100.;// convert to cm
		const int Z_cells = ptree.get("domain.yasp.nZ",(int)2000) ;
		L[1] = Y;
		L[2] = Z;
		N[1] = Y_cells ;
		N[2] = Z_cells ;
#endif
		std::bitset<dim> periodic(false);
		int overlap=1;
		Dune::YaspGrid<dim> grid(L,N,periodic,overlap,helper.getCommunicator());
		typedef Dune::YaspGrid<dim>::LeafGridView GV;
		const GV& gv=grid.leafGridView();
        grid.loadBalance();

		/**************************************************************************************************/
		// DRIVER
		driver( gv,
				ptree,
        		OUTPUT_PATH,
				helper);
		/**************************************************************************************************/
#endif

	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
