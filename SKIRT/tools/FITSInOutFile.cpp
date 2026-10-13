/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FITSInOutFile.hpp"
#include "FITSInOut.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "H5Lib.hpp"
#include "Log.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

void FITSInOutFile::read(const SimulationItem* item, const string& filename, Array& data, int& nx, int& ny, int& nz)
{
    Log* log = item->find<Log>();

    // determine the candidate locations for the input file: a regular FITS file and a bundle in an HDF5 file
    auto [filePath, bundlePath] = item->find<FilePaths>()->input(filename);

    // determine whether the regular file exists; its name may be followed by an extension name between brackets
    bool fileExists = false;
    if (!filePath.empty())
    {
        auto bracket = filePath.find('[', filePath.find_last_of("/\\") + 1);
        fileExists = System::isFile(filePath.substr(0, bracket));
    }

    // read from the regular file if it exists; if not, read from the bundle if there is a bundle candidate
    if (fileExists || bundlePath.empty())
    {
        FITSInOut::read(filePath, data, nx, ny, nz);
        log->info(item->typeAndName() + " read FITS file " + filePath);
    }
    else
    {
        H5DatasetR dataset = H5Lib::openBundle(bundlePath).getDataset("data");

        // the dimensions are stored with the index that varies least rapidly listed first
        vector<size_t> shape = dataset.shape();
        size_t rank = shape.size();
        if (rank != 2 && rank != 3)
            throw FATALERROR("The dataset 'data' in HDF5 bundle " + bundlePath
                             + " must have two or three dimensions rather than " + std::to_string(rank));
        for (size_t size : shape)
            if (size > static_cast<size_t>(INT_MAX))
                throw FATALERROR("The dataset 'data' in HDF5 bundle " + bundlePath + " is too large");
        nx = static_cast<int>(shape[rank - 1]);
        ny = static_cast<int>(shape[rank - 2]);
        nz = rank == 3 ? static_cast<int>(shape[0]) : 1;

        data = dataset.readDouble();
        log->info(item->typeAndName() + " read FITS bundle " + bundlePath);
    }

    // log the data dimensions
    if (nz > 1) log->info("  Data cube has " + std::to_string(nz) + " frames");
    log->info("  Frame dimensions are " + std::to_string(nx) + " x " + std::to_string(ny));
}

////////////////////////////////////////////////////////////////////
