#include "mutation++.h"
#include <string>
#ifndef INPUT_H
#define INPUT_H

using namespace std;
using namespace Mutation::Utilities;

/// Input class providing access to global options.
class Input
{
public:
    Input();

    Input(const std::string &filename);

    Input(const char *filename);

    Input(IO::XmlElement &element);

    ~Input();

    void setDefaultOptions();

    /// Gets the data directory.
    const std::string &dataDirectory();

    /// Sets the data directory.
    void dataDirectory(const std::string &dir);

    /// Gets the data directory.
    const std::string &workingDirectory();

    /// Sets the data directory.
    void workingDirectory(const std::string &dir);

    /// Gets the mixture name.
    const std::string &mixturename();

    /// Sets the mixture name.
    void mixturename(const std::string &name);

    /// Gets the init pressure.
    const double &pressureinit();

    /// Sets the init pressure.
    void pressureinit(const double &pressure);

    /// Gets the init temperature.
    const double &temperatureinit();

    /// Sets the init temperature.
    void temperatureinit(const double &temperature);

    /// Gets the init distance.
    const double &distanceinit();

    /// Sets the init temperature.
    void distanceinit(const double &distance);

    /// Gets the init maxstep.
    const int &maxstepinit();

    /// Sets the init temperature.
    void maxstepinit(const int &maxstep);

    /// Gets the init increment species density.
    const double &pertminit();

    /// Sets the init increment species density.
    void pertminit(const double &pert_m);

    /// Gets the init increment temperature.
    const double &pertTinit();

    /// Sets the init increment temperature.
    void pertTinit(const double &pert_T);

        /// Gets the tolerance for Newton solver.
    const double &tolinit();

    /// Sets the tolerance for Newton solver.
    void tolinit(const double &tol);

    /// Switch to iteration history of Newton solver.
    const std::string& iNewtonhistory();

    /// Switch to iteration history of Newton solver.
    void iNewtonhistory(const std::string& inewtonhistory);

    /// Gets the file separator character.
    const char separator();

    /// Sets the file separator character.
    void separator(char sep);
    void loadFromFile();

    void loadFromXmlElement(IO::XmlElement &element);

private:
    /**
     * Returns the value for the environment variable key.  If the key is not
     * found, then an empty string is returned.
     */
    std::string getEnvironmentVariable(const std::string &key)
    {
        char *value = std::getenv(key.c_str());
        return std::string(value == NULL ? "" : value);
    }

private:
    /// Data directory
    std::string m_input_file;
    /// Data directory
    std::string m_data_directory;

    /// Working directory
    std::string m_working_directory;

    ///Mixture name
    std::string m_mixture_name;

    ///Pressure init
    double m_pressure_init;

    ///temperature init
    double m_temperature_init;

    ///distance init
    double m_distance;

    ///Max step for Newton iteration
    int m_maxstep;

    ///increment for species density performs in Newton solver
    double m_pert_m;

    ///increment for temperature performs in Newton solver
    double m_pert_T;

    ///Tolerance in Newton solver
    double m_tol;

    ///Whether write history in Newton solver
    std::string m_iNewtonhistory;

    /// File separator character
    char m_separator;

}; // class Input

#endif // Input_H