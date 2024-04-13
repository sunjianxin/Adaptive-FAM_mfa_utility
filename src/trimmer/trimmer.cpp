//--------------------------------------------------------------
// example of evaluating and optionally differentiating a single point of an MFA
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include "opts.h"

#include "block.hpp"

using namespace std;

int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file

    // default command line arguments
    int MAX_DIM     = 10;                       // maximum domain dimensionality (temporary, only for input parameter)
    vector<real_t> param(MAX_DIM);              // point to evaluate, initialized to 0s by default
    int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    bool help;                                  // show help

    // get command line arguments
    opts::Options ops;
    ops >> opts::Option('p', "param",   param,   " parameters of point to evaluate");
    ops >> opts::Option('d', "deriv",   deriv,   " which derivative to take (1 = 1st, 2 = 2nd, ...)");
    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // echo args
    fprintf(stderr, "\n--------- Input arguments ----------\n");
    cerr << "param (first two dims.) = [ " << param[0] << " " << param[1] << " ]" << endl;
    cerr << "deriv   = "    << deriv << endl;
#ifdef MFA_TBB
    cerr << "threading: TBB" << endl;
#endif
#ifdef MFA_KOKKOS
    cerr << "threading: Kokkos" << endl;
#endif
#ifdef MFA_SYCL
    cerr << "threading: SYCL" << endl;
#endif
#ifdef MFA_SERIAL
    cerr << "threading: serial" << endl;
#endif
    fprintf(stderr, "-------------------------------------\n\n");

    // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX");               // used for blocks to be moved out of core
    diy::Master      master(world,
            -1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy,
            &storage,
            &Block<real_t>::save,
            &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), -1);   // number of blocks set by read_blocks()

    // read MFA model
    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
    int nblocks = master.size();
    std::cout << nblocks << " blocks read from file "<< infile << "\n";
    int dom_dim = master.block<Block<real_t>>(0)->dom_dim;  // dimensionality of input domain
    int pt_dim  = master.block<Block<real_t>>(0)->pt_dim;   // dimensionality of output point

	std::cout << dom_dim << std::endl;
	std::cout << pt_dim << std::endl;

    void* bblock = master.block(0);
	std::cout << "0" << std::endl;
	Block<real_t>* b = (Block<real_t>*)(bblock);
	std::cout << "1" << std::endl;
	
	// mfa::MFA_Data<real_t>* mfa_data_2 = b->vars[0].mfa_data;
	std::cout << "size: " << b->vars.size() << std::endl;
	// std::cout << b->vars.at(0).mfa_data->p(0) << std::endl;
	// std::cout << b->vars[0].mfa_data->p(0) << std::endl;
	// std::cout << b->vars[0].mfa_data->tmesh.all_knots[0].size() << std::endl;
	std::cout << "1.5" << std::endl;
#if 1
	mfa::MFA_Data<real_t>& mfa_data = *(b->vars[0].mfa_data);
	std::cout << "2" << std::endl;
	// TensorProduct<real_t>&  tt = mfa_data_2->tmesh.tensor_prods[0];
	TensorProduct<real_t>&  tt = mfa_data.tmesh.tensor_prods[0];
	std::cout << "3" << std::endl;
	// mfa::Decoder<real_t> decoder(*mfa_data_2, 0);
	mfa::Decoder<real_t> decoder(mfa_data, 0);
	std::cout << "4" << std::endl;
  	mfa::FastDecodeInfo<real_t> di(decoder);
	VectorX<real_t> extents = b->bounds_maxs - b->bounds_mins;

	std::cout << "extents: " << extents << std::endl;

	di.ResizeDers(1);
	VectorX<real_t> param_fast(3);
	VectorX<real_t> out_pt_fast(1);
    param_fast[0] = 0.3;
    param_fast[1] = 0.3;
    param_fast[2] = 0.3;
	decoder.FastVolPt(param_fast, out_pt_fast, di, tt);
	std::cout << "fast value queried: " << out_pt_fast[0] << std::endl;

	VectorX<real_t> grad(3);
	decoder.FastGrad(param_fast, di, tt, grad);
	std::cout << "fast grad queried: " << grad << std::endl;
	

	unsigned int knots_size = mfa_data.tmesh.all_knots[0].size();
	unsigned int degree = mfa_data.p(0);
	unsigned int ctrl_pts_size = tt.ctrl_pts.rows(); // = (knots_size - 2- 1)^3
	unsigned int ctrl_pts_size_per_dim = std::cbrt(ctrl_pts_size);

	
	std::cout << "-------Info-------" << std::endl;
	std::cout << "degree: " << mfa_data.p(0) << ", " << mfa_data.p(1) << ", "<< mfa_data.p(2) << std::endl;
	std::cout << "knots size: " << mfa_data.tmesh.all_knots[0].size() << "+" << mfa_data.tmesh.all_knots[1].size() << "+" << mfa_data.tmesh.all_knots[2].size() << std::endl;
	std::cout << "ctrl pts: " << ctrl_pts_size_per_dim << "x" << ctrl_pts_size_per_dim << "x" << ctrl_pts_size_per_dim << std::endl;
	// std::cout << "Info: ctrl pts: " << tt.ctrl_pts.rows() << "x" << tt.ctrl_pts.cols() << std::endl;

	std::cout << "-------Verification-------" << std::endl;
	std::cout << "all_knots[0][34]: " << mfa_data.tmesh.all_knots[1][34] << std::endl;
	std::cout << "tt.ctrl_pts(0): " << tt.ctrl_pts(0) << std::endl;
	std::cout << "tt.ctrl_pts(100000): " << tt.ctrl_pts(100000) << std::endl;
	std::cout << "tt.ctrl_pts(ctrl_pts_size - 1): " << tt.ctrl_pts(ctrl_pts_size - 1) << std::endl;

#endif


#if 1

	std::ofstream wf("test.mfab", std::ios::out | std::ios::binary);
   	if(!wf) {
    	cout << "Cannot open file!" << endl;
    	return 1;
  	}

	// Write degree
    wf.write((char *)&degree, sizeof(unsigned int));
    // Write number of knots on each dimension
    wf.write((char *)&knots_size, sizeof(unsigned int));
	// Write all the knots
   	for(unsigned int i = 0; i < knots_size; i++)
    	wf.write((char *)&mfa_data.tmesh.all_knots[0][i], sizeof(float));
	for(unsigned int i = 0; i < knots_size; i++)
    	wf.write((char *)&mfa_data.tmesh.all_knots[1][i], sizeof(float));
   	for(unsigned int i = 0; i < knots_size; i++)
    	wf.write((char *)&mfa_data.tmesh.all_knots[2][i], sizeof(float));
	// Write all ctrl points
	for(unsigned int i = 0; i < ctrl_pts_size; i++)
    	wf.write((char *)&tt.ctrl_pts(i), sizeof(float));
   	wf.close();

	// Only write control points
	std::ofstream wff("test.cpts", std::ios::out | std::ios::binary);
   	if(!wff) {
    	cout << "Cannot open file!" << endl;
    	return 1;
  	}
	// Write all ctrl points
	for(unsigned int i = 0; i < ctrl_pts_size; i++)
    	wff.write((char *)&tt.ctrl_pts(i), sizeof(float));
   	wff.close();

	std::cout << "mfab file saved, size: " << (1 + 1 + mfa_data.tmesh.all_knots[0].size()+mfa_data.tmesh.all_knots[1].size()+mfa_data.tmesh.all_knots[2].size() + ctrl_pts_size)*4 << " bytes" << std::endl;


#if 0 // Verification
	std::ifstream fd;
	std::string mfab_file = "test.mfab";
	fd.open(mfab_file, std::ios::binary);
	assert(fd);

	std::cout << "-------Readback:-------" << std::endl;
	unsigned int degree_readback;
	fd.read(reinterpret_cast<char *>(&degree_readback), sizeof(unsigned int));
	std::cout << "degree_readback: " << degree_readback << std::endl;


	unsigned int knots_size_readback;
	fd.read(reinterpret_cast<char *>(&knots_size_readback), sizeof(unsigned int));
	std::cout << "knots_size_readback: " << knots_size_readback << std::endl;

	float all_knots_readback[79*3];
	fd.read(reinterpret_cast<char *>(all_knots_readback), sizeof(float)*79*3);
	std::cout << "all_knots_readback[0]: " << all_knots_readback[0] << std::endl;
	std::cout << "all_knots_readback[34]: " << all_knots_readback[34] << std::endl;
	std::cout << "all_knots_readback[79*3 - 1]: " << all_knots_readback[79*3 - 1] << std::endl;

	float ctrl_pts_readback[ctrl_pts_size];
	fd.read(reinterpret_cast<char *>(ctrl_pts_readback), sizeof(float)*ctrl_pts_size);
	std::cout << "ctrl_pts_readback[0]: " << ctrl_pts_readback[0] << std::endl;
	std::cout << "ctrl_pts_readback[100000]: " << ctrl_pts_readback[100000] << std::endl;
	std::cout << "ctrl_pts_readback[ctrl_pts_size - 1]: " << ctrl_pts_readback[ctrl_pts_size - 1] << std::endl;
#endif
#endif

	// Data formate: knot_size(unsigned int, 4 bytes); knot_value(float, 4 bytes)

}
