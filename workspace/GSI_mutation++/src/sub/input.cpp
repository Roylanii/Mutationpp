#include "input.h"
using namespace std;
using namespace Mutation::Utilities;
Input::Input()
{
    setDefaultOptions();
}
void Input::setDefaultOptions()
{
    m_input_file = "input.xml";
    m_data_directory = getEnvironmentVariable("MPP_DATA_DIRECTORY");
    m_working_directory = "data";
    m_mixture_name = "";
    m_pressure_init = 101325;
    m_temperature_init = 3000;
    m_distance = 1e-3;
    m_maxstep = 100;
#ifdef _WIN32
    m_separator = '\\';
#else
    m_separator = '/';
#endif
}

Input::Input(const std::string &filename)
{
    setDefaultOptions();
    m_input_file = filename;
    loadFromFile();
}

Input::Input(const char *filename)
{
    setDefaultOptions();
    m_input_file = filename;
    loadFromFile();
}

Input::Input(IO::XmlElement &element)
{
    // Initalize to the default options
    setDefaultOptions();
    loadFromXmlElement(element);
}

Input::~Input() {}
const std::string &Input::dataDirectory()
{
    return m_data_directory;
}
void Input::dataDirectory(const std::string &dir)
{
    m_data_directory = dir;
}
const std::string &Input::workingDirectory()
{
    return m_working_directory;
}
void Input::workingDirectory(const std::string &dir)
{
    m_working_directory = dir;
}
const std::string &Input::mixturename()
{
    return m_mixture_name;
}

/// Sets the mixture name.
void Input::mixturename(const std::string &name)
{
    m_mixture_name = name;
}

/// Gets the init pressure.
const double &Input::pressureinit()
{
    return m_pressure_init;
}

/// Sets the init pressure.
void Input::pressureinit(const double &pressure)
{
    m_pressure_init = pressure;
}

/// Gets the init temperature.
const double &Input::temperatureinit()
{
    return m_temperature_init;
}

/// Sets the init temperature.
void Input::temperatureinit(const double &temperature)
{
    m_temperature_init = temperature;
}

/// Gets the init distance.
const double &Input::distanceinit()
{
    return m_distance;
}

/// Sets the init distance.
void Input::distanceinit(const double &distance)
{
    m_distance = distance;
}

/// Gets the init maxstep.
const int &Input::maxstepinit()
{
    return m_maxstep;
}

/// Sets the init maxstep.
void Input::maxstepinit(const int &maxstep)
{
    m_maxstep = maxstep;
}

/// Gets the file separator character.
const char Input::separator()
{
    return m_separator;
}

/// Sets the file separator character.
void Input::separator(char sep)
{
    m_separator = sep;
}

void Input::loadFromFile()
{

    // Now load the XML file
    IO::XmlDocument input_doc(m_input_file);
    IO::XmlElement root = input_doc.root();

    // Now we can load from the mixture XmlElement
    loadFromXmlElement(root);
}

void Input::loadFromXmlElement(IO::XmlElement &element)
{
    if (element.tag() != "input")
        element.parseError(
            "XmlElement is not of 'input' type!");
    
    IO::XmlElement::const_iterator iter;
    for (iter = element.begin(); iter != element.end(); ++iter)
    {

    }

    element.getAttribute<std::string>("data_directory", m_data_directory, m_data_directory);

    element.getAttribute<std::string>("working_directory", m_working_directory, m_working_directory);

    element.getAttribute<std::string>("mixture_name", m_mixture_name, "A mixture name must be given");

    // Get the thermal conductivity algorithm
    element.getAttribute<double>("pressure_init", m_pressure_init, "pressure_init must be given");

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<double>("temperature_init", m_temperature_init, "temperature_init must be given");

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<double>("distance", m_distance, "distance must be given");

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<int>("maxstep", m_maxstep, "maxstep must be given");

    // IO::XmlElement::const_iterator iter;
    // for (iter = element.begin(); iter != element.end(); ++iter)
    // {
    //     // Load the species list
    //     if (iter->tag() == "species")
    //         m_species_descriptor = String::trim(iter->text());
    //     else if (iter->tag() == "element_compositions")
    //         loadElementCompositions(*iter);
    // }
}
