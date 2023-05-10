/*
 * boundaryvaluefunction.hh
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

#ifndef OPERATORS_BOUNDARYVALUEFUNCTION_HH_
#define OPERATORS_BOUNDARYVALUEFUNCTION_HH_

template<typename GV, typename PGMap, class ProblemBCT >
class BoundaryTypes
{
private :
	const GV& gv ;
	const PGMap& pgmap ;
	const ProblemBCT& bct;

public :

	// ! construct from gridview
	BoundaryTypes ( const GV& gv_, const PGMap& pgmap_, const ProblemBCT& bct_ )
	: gv ( gv_ ), pgmap ( pgmap_ ), bct( bct_ ) {}

	// ! return bctype for Pw at point on intersection
	template<typename I>
	inline void evaluate ( I& i ,
						   const Dune::FieldVector<double,GV::dimension-1>& xlocal ,
						   double t,
						   double dt,
						   int tag,
						   int& y ) const {

	    y=bct.evaluate(i,xlocal,t,dt,tag);
		return ;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }

} ;



template<typename GV, typename PGMap, class ProblemBCV >
class BoundaryValues
: public Dune::PDELab::BoundaryGridFunctionBase< Dune::PDELab::BoundaryGridFunctionTraits<GV, double , 1 , Dune::FieldVector<double,1> >,
  	  	  	  	  	  	  	  	  	  	  	  	 BoundaryValues<GV,PGMap,ProblemBCV> >
{
public :

	// ! construct from gridview
	BoundaryValues ( const GV& gv_, const PGMap& pgmap_, const ProblemBCV& bcv_ )
	: gv ( gv_ ) ,
	  pgmap ( pgmap_ ),
	  bcv( bcv_)
	{}

	// ! return bcvalue of species 'tag' at point on intersection
	template<typename I>
	inline void evaluate (	I& intersection,
							const Dune::FieldVector<double,GV::dimension-1>& xlocal ,
							double t,
							double dt,
							int tag,
							double& y ) const {

	    y = bcv.evaluate(intersection,xlocal,t,dt,tag) ;
	    return ;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }

private :
	const GV& gv ;
	const PGMap& pgmap ;
	const ProblemBCV& bcv;
} ;

#endif /* OPERATORS_BOUNDARYVALUEFUNCTION_HH_ */
