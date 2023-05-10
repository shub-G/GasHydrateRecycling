/*
 * initialvaluefunction.hh
 *
 *  Created on: Aug 30, 2021
 *      Author: sgupta
 */

#ifndef OPERATORS_INITIALVALUEFUNCTION_HH_
#define OPERATORS_INITIALVALUEFUNCTION_HH_

template< class GV,
		  class ProblemICV,
		  class GFS,
		  typename U >
class InitialValues{
private:
	const GV&  gv;
	const ProblemICV& icv;
	GFS gfs;
	U *u;

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	InitialValues(	const 	GV& gv_,
					const ProblemICV icv_,
					GFS	gfs_,
					U	*u_ )
	: gv(gv_),
	  icv(icv_),
	  gfs(gfs_),
	  u(u_)
	{}

	virtual ~InitialValues()
	{}

	void evaluate(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS > LFS;
		LFS lfs(gfs);

		typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
		LFSCache lfs_cache(lfs);
		typedef typename U::template LocalView<LFSCache> VectorView;
		VectorView u_view( (*u) );

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);

	        lfs.bind(*self);
			lfs_cache.update();
			u_view.bind(lfs_cache);
	        std::vector<double> ul(lfs.size());

	        for(int i = 0. ; i < lfs.size() ; i++){
	        	ul[lfs.child(i).localIndex(0)] = icv.evaluate(cell,cell_center_local,i) ;
	        }

			u_view.write( ul );
			u_view.commit();
			u_view.unbind();

		}//END:iterate over each volume

	}
};

#endif /* OPERATORS_INITIALVALUEFUNCTION_HH_ */
