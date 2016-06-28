#ifndef __STRUCTURE_FACTOR_H__
#define __STRUCTURE_FACTOR_H__

/*! \file StructureFactor.h
    \brief Declares the StructureFactor class
 */
#include <hoomd/hoomd.h>

#include </home/morsedc/tghasi/Sq-plugin/util/AutoCorrelation.tpp>
#include </home/morsedc/tghasi/Sq-plugin/util/DArray.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <complex>

// need to declare these classes with __host__ __device__ qualifiers when building in nvcc
// HOSTDEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define HOSTDEVICE __host__ __device__
#else
#define HOSTDEVICE
#endif

typedef std::complex<double> Complex;
using namespace Util;

//! Collective variable for studying phase transitions in block copolymer systems
class StructureFactor : public Analyzer
    {
    public:
        /*! Constructs the structure factor plugin
            \param sysdef The system definition
            \param mode The per-type coefficients of the Fourier mode
            \param lattice_vectors The Miller indices of the mode vector
            \param filename_sq Filename for the log file
            \param filename_vh Filename for the Van-Hove results
            \param overwrite Whether the log file should be overwritten
         */
        StructureFactor(boost::shared_ptr<SystemDefinition> sysdef,
                               const std::vector<Scalar>& mode,
                               const std::vector<int3>& lattice_vectors,
                               const unsigned int& vanhove_endtime,
                               const std::string& filename_sq,
                               const std::string& filename_vh,
                               bool overwrite=false
                               );
        virtual ~StructureFactor() {}

        void analyze(unsigned int timestep);

    protected:
        std::string m_filename_sq;            //!< The name of the log file
        std::string m_filename_vh;            //!< The name of the Van-Hove results
        std::ofstream m_file_sq;              //!< The log file handle
        std::ofstream m_file_vh;              //!< The VanHove file handle
        std::string m_delimiter;              //!< Record delimiter
        bool m_appending;                     //!< Whether we are appending to the log file
        std::vector<int3> m_lattice_vectors;  //!< Stores the list of miller indices
        unsigned int m_vanhove_endtime; //!< Stores the last step of VanHove
        std::vector<Scalar> m_mode;           //!< Stores the per-type mode coefficients

        GPUArray<Scalar3> m_wave_vectors;     //!< GPUArray of wave vectors
        GPUArray<Scalar2> m_fourier_modes;    //!< Fourier modes
        GPUArray<Scalar> m_phases;            //!< Phase shifts
        Util::DArray< Util::AutoCorrelation<Complex, Complex> > accumulator_;

        //! Helper function to update the wave vectors
        void calculateWaveVectors();

        //! Calculates the current value of the collective variable
        virtual void computeFourierModes(unsigned int timestep);

    private:
        void openOutputFiles();

    };

// ------------ Vector math functions --------------------------
//! Comparison operator needed for export of std::vector<int3>
HOSTDEVICE inline bool operator== (const int3 &a, const int3 &b)
    {
    return (a.x == b.x &&
            a.y == b.y &&
            a.z == b.z);
    }

//! Export StructureFactor to python
void export_StructureFactor();

#endif // __STRUCTURE_FACTOR_H__
