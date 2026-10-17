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
#include "ProcessManager.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the grid points of a regularly spaced axis with the specified number of points, increment, and center
    // value, computed as they follow from the CRPIX, CRVAL, and CDELT cards in a FITS file
    Array regularGrid(int n, double increment, double center)
    {
        Array grid(n);
        double refpix = (n + 1) / 2.;
        for (int i = 0; i != n; ++i) grid[i] = center + increment * (i + 1 - refpix);
        return grid;
    }

    // creates a dataset with the specified name, values, unit, and shape in the specified bundle
    void writeDataset(const H5BundleW& bundle, const string& name, const Array& values, const string& unit,
                      const vector<size_t>& shape = vector<size_t>())
    {
        bundle.createDoubleDataset(name, values, shape).setStringAttribute("unit", unit);
    }

    // writes the description and the observer information, if any, as attributes of the specified bundle
    void setBundleAttributes(const H5BundleW& bundle, const string& description, const FITSInOut::ObserverInfo* obsInfo)
    {
        bundle.setStringAttribute("description", description);
        if (obsInfo)
        {
            bundle.setDoubleAttribute("inclination", obsInfo->inclination);
            bundle.setDoubleAttribute("azimuth", obsInfo->azimuth);
            bundle.setDoubleAttribute("roll", obsInfo->roll);
            bundle.setDoubleAttribute("redshift", obsInfo->redshift);
            bundle.setDoubleAttribute("luminosity_distance", obsInfo->luminosityDistance);
            bundle.setDoubleAttribute("angular_diameter_distance", obsInfo->angularDiameterDistance);
            bundle.setStringAttribute("distance_unit", obsInfo->distanceUnits);
        }
    }
}

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
        // if there is no bundle either, report both locations
        item->find<FilePaths>()->verifyBundleInput(bundlePath);
        if (!filePath.empty() && !H5Lib::hasBundle(bundlePath))
            throw FATALERROR("Could not find the FITS input file " + filename + ": there is no regular file " + filePath
                             + " and no bundle " + bundlePath);

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

void FITSInOutFile::write(const SimulationItem* item, const string& description, const string& filename,
                          const Array& data, const string& dataUnits, int nx, int ny, double incx, double incy,
                          double xc, double yc, const string& xyUnits, const Array& z, const string& zUnits,
                          const FITSInOut::ObserverInfo* obsInfo)
{
    // only the root process writes the output file
    if (!ProcessManager::isRoot()) return;

    // determine the candidate locations for the output file: a regular FITS file and a bundle in an HDF5 file
    auto [filePath, bundlePath] = item->find<FilePaths>()->output(filename + ".fits");

    // write a regular FITS file if there is no bundle candidate
    if (bundlePath.empty())
    {
        FITSInOut::write(filePath, data, dataUnits, nx, ny, incx, incy, xc, yc, xyUnits, z, zUnits, obsInfo);
        item->find<Log>()->info(item->typeAndName() + " wrote " + description + " to FITS file " + filePath);
    }

    // otherwise write a bundle
    else
    {
        size_t nz = z.size();
        if (data.size() != static_cast<size_t>(nx) * static_cast<size_t>(ny) * (nz ? nz : 1))
            throw FATALERROR("Inconsistent data size when creating FITS file bundle " + bundlePath);

        item->find<FilePaths>()->verifyBundleOutput(bundlePath);
        H5BundleW bundle = H5Lib::createBundle(bundlePath);
        setBundleAttributes(bundle, description, obsInfo);
        writeDataset(bundle, "x", regularGrid(nx, incx, xc), xyUnits);
        writeDataset(bundle, "y", regularGrid(ny, incy, yc), xyUnits);
        if (nz) writeDataset(bundle, "z", z, zUnits);

        // a data cube with a single frame is stored as a frame, next to a z-axis with a single grid point
        vector<size_t> shape;
        if (nz > 1) shape.push_back(nz);
        shape.push_back(ny);
        shape.push_back(nx);
        writeDataset(bundle, "data", data, dataUnits, shape);

        item->find<Log>()->info(item->typeAndName() + " wrote " + description + " to bundle " + bundlePath);
    }
}

////////////////////////////////////////////////////////////////////

void FITSInOutFile::writeMap(const SimulationItem* item, const string& description, const string& filename,
                             const Array& data, const string& dataUnits, const Array& x, const Array& y,
                             const string& xUnits, const string& yUnits, const FITSInOut::ObserverInfo* obsInfo)
{
    // only the root process writes the output file
    if (!ProcessManager::isRoot()) return;

    // determine the candidate locations for the output file: a regular FITS file and a bundle in an HDF5 file
    auto [filePath, bundlePath] = item->find<FilePaths>()->output(filename + ".fits");

    // write a regular FITS file if there is no bundle candidate
    if (bundlePath.empty())
    {
        FITSInOut::writeMap(filePath, data, dataUnits, x, y, xUnits, yUnits, obsInfo);
        item->find<Log>()->info(item->typeAndName() + " wrote " + description + " to FITS file " + filePath);
    }

    // otherwise write a bundle
    else
    {
        if (data.size() != x.size() * y.size())
            throw FATALERROR("Inconsistent data size when creating FITS file bundle " + bundlePath);

        item->find<FilePaths>()->verifyBundleOutput(bundlePath);
        H5BundleW bundle = H5Lib::createBundle(bundlePath);
        setBundleAttributes(bundle, description, obsInfo);
        writeDataset(bundle, "x", x, xUnits);
        writeDataset(bundle, "y", y, yUnits);
        writeDataset(bundle, "data", data, dataUnits, {y.size(), x.size()});

        item->find<Log>()->info(item->typeAndName() + " wrote " + description + " to bundle " + bundlePath);
    }
}

////////////////////////////////////////////////////////////////////
