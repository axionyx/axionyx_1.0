#include "Nyx.H"
#include "Nyx_F.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;


void Nyx::swfft_test(MultiFab& rhs, MultiFab& soln, Geometry& geom, int verbose)
{
    const BoxArray& ba = soln.boxArray();
    amrex::Print() << "BA " << ba << std::endl;
    const DistributionMapping& dm = soln.DistributionMap();

    if (rhs.nGrow() != 0 || soln.nGrow() != 0)
       amrex::Error("Current implementation requires that both rhs and soln have no ghost cells");

    // Define pi and (two pi) here
    Real  pi = 4 * std::atan(1.0);
    Real tpi = 2 * pi;

    // We assume that all grids have the same size hence
    // we have the same nx,ny,nz on all ranks
    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
    int nz = ba[0].size()[2];

    Box domain(geom.Domain());

    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
    int nbz = domain.length(2) / nz;
    int nboxes = nbx * nby * nbz;
    if (nboxes != ba.size())
       amrex::Error("NBOXES NOT COMPUTED CORRECTLY");
    amrex::Print() << "Number of boxes:\t" << nboxes << std::endl;

    Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    DistributionMapping dmap = rhs.DistributionMap();

    for (int ib = 0; ib < nboxes; ++ib)
    {
        int i = ba[ib].smallEnd(0) / nx;
        int j = ba[ib].smallEnd(1) / ny;
        int k = ba[ib].smallEnd(2) / nz;

        // This would be the "correct" local index if the data wasn't being transformed
        // int local_index = k*nbx*nby + j*nbx + i;

        // This is what we pass to dfft to compensate for the Fortran ordering
        //      of amrex data in MultiFabs.
        int local_index = i*nby*nbz + j*nbz + k;

        rank_mapping[local_index] = dmap[ib];
        if (verbose)
          amrex::Print() << "LOADING RANK NUMBER " << dmap[ib] << " FOR GRID NUMBER " << ib
                         << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }


    Real h = geom.CellSize(0);
    Real hsq = h*h;

    Real start_time = amrex::second();

    // Assume for now that nx = ny = nz
    int Ndims[3] = { nbz, nby, nbx };
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d);

    //some constants which need to be passed from the init file
    double hbar = 1.;
    double k0 = 1.;
    double m = 1.;
    double dt = 0.1;

    //get the potential
    MultiFab& phi = get_new_data(PhiGrav_Type);

    //loop for the drift
    for (MFIter mfi(rhs,false); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        size_t local_size  = dfft.local_size();

        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

        a.resize(nx*ny*nz);
        b.resize(nx*ny*nz);

        dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

        // *******************************************
        // Copy real data from Rhs into real part of a -- no ghost cells and
        // put into C++ ordering (not Fortran)
        // *******************************************
        complex_t zero(0.0, 0.0);
        size_t local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
            for(size_t j=0; j<(size_t)ny; j++) {
                for(size_t i=0; i<(size_t)nx; i++) {

                    complex_t temp(rhs[mfi].dataPtr()[local_indx],0.); // the imaginary part should be filled too
                    a[local_indx] = temp;
                    local_indx++;

                }
            }
        }
        // *******************************************
        // Same thing as above, for the potentail -only need the real part
        // put into C++ ordering (not Fortran)
        // *******************************************
        local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
            for(size_t j=0; j<(size_t)ny; j++) {
                for(size_t i=0; i<(size_t)nx; i++) {

                    complex_t temp(phi[mfi].dataPtr()[local_indx],0.);
                    b[local_indx] = temp;
                    local_indx++;

                }
            }
        }



        //  *******************************************
        //  Compute the forward transform
        //  *******************************************
        dfft.forward(&a[0]);

        //  *******************************************
        //  drift - one full time step
        //  *******************************************
        const int *self = dfft.self_kspace();
        const int *local_ng = dfft.local_ng_kspace();
        const int *global_ng = dfft.global_ng();
        const std::complex<double> imagi(0.0,1.0);

        local_indx = 0;
        for(size_t k=0; k<(size_t)local_ng[0]; k++) {
            size_t global_k = local_ng[2]*self[0] + k; //maybe need another condition for if (k < local_ng[2]/2)?

            for(size_t j=0; j<(size_t)local_ng[1]; j++) {
                size_t global_j = local_ng[1]*self[1] + j;

                for(size_t i=0; i<(size_t)local_ng[2]; i++) {
                    size_t global_i = local_ng[0]*self[2] + i;

                    if (global_i == 0 && global_j == 0 && global_k == 0) {
                        a[local_indx] = 0;
                    }
                    else {
                        double k2 = k0 * k0 * ( double(global_i) * double(global_i)
                                 + double(global_j) * double(global_j)
                                 + double(global_k) * double(global_k) );

                        std::complex<double> cpsi(std::real(a[local_indx]),std::imag(a[local_indx]));
                        cpsi *= std::exp( -imagi * hbar * k2/2.0/m  * dt );

                        a[local_indx]= { cpsi.real() /local_ng[0]/local_ng[1]/local_ng[2], cpsi.imag() /local_ng[0]/local_ng[1]/local_ng[2] };

                    }
                    local_indx++;
                }
            }
        }

        // *******************************************
        // Compute the backward transformdfft.global_size
        // *******************************************
        dfft.backward(&a[0]);


        //  *******************************************
        //  kick - one full time step
        //  *******************************************
        local_indx = 0;
    	for(size_t k=0; k<(size_t)local_ng[2]; k++) {
    		for(size_t j=0; j<(size_t)local_ng[1]; j++) {
    			for(size_t i=0; i<(size_t)local_ng[0]; i++) {

                    std::complex<double> cpsi(std::real(a[local_indx]),std::imag(a[local_indx]));

    				cpsi *= std::exp( -imagi * b[local_indx] / hbar  * dt );

                    a[local_indx] = { cpsi.real(), cpsi.imag() };
                    local_indx++;
    			}
    		}
    	}


        size_t global_size  = dfft.global_size();
        std::cout << "GOBAL SIZE " << global_size << std::endl;
        double fac = hsq / global_size;

        local_indx = 0;
        for(size_t k=0; k<(size_t)nz; k++) {
           for(size_t j=0; j<(size_t)ny; j++) {
               for(size_t i=0; i<(size_t)nx; i++) {

                   // here, we need to keep the imag part of a as well
                   soln[mfi].dataPtr()[local_indx] = std::real(a[local_indx])/global_size;
                   local_indx++;

                }
            }
        }

    }//MFIter loop
}//swfft_test function
