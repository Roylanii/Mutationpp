#include "input.h"
#include <getopt.h>
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
    m_maxiteration = 50;
    m_subiteration = 5;
    m_pert_m = 1e-2;
    m_pert_T = 1.0;
    m_eps = 1.0e-12;
    m_iNewtonhistory = "false";
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
const int &Input::maxiterationsinit()
{
    return m_maxiteration;
}

const int &Input::subiterationsinit()
{
    return m_subiteration;
}

/// Sets the init maxstep.
void Input::maxiterationsinit(const int &maxstep)
{
    m_maxiteration = maxstep;
}

void Input::subiterationsinit(const int &substep)
{
    m_subiteration = substep;
}

/// Gets the init increment species density.
const double &Input::pertminit()
{
    return m_pert_m;
}

/// Sets the init increment species density.
void Input::pertminit(const double& pert_m)
{
    m_pert_m = pert_m;
}

/// Gets the init increment temperature.
const double &Input::pertTinit()
{
    return m_pert_T;
}

/// Sets the init increment temperature.
void Input::pertTinit(const double& pert_T)
{
    m_pert_T = pert_T;
}

const double &Input::tolinit()
{
    return m_eps;
}

/// Sets the init increment temperature.
void Input::tolinit(const double& tol)
{
    m_eps = tol;
}

const std::string &Input::iNewtonhistory()
{
    return m_iNewtonhistory;
}

/// Sets the init increment temperature.
void Input::iNewtonhistory(const std::string& inewtonhistory)
{
    m_iNewtonhistory = inewtonhistory;
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
    element.getAttribute<int>("maxiterations", m_maxiteration, m_maxiteration);

    element.getAttribute<int>("subiterations", m_subiteration, m_subiteration);

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<double>("pert_m", m_pert_m, m_pert_m);

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<double>("pert_T", m_pert_T, m_pert_T);

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<double>("resnorm", m_eps, m_eps);

    // Get the type of Gas-Surface Interaction for the wall
    element.getAttribute<std::string>("iNewtonhistory", m_iNewtonhistory, m_iNewtonhistory);

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

void getOptLong(int argc, char *argv[], std::map<int, std::string> &input)
{
    int opt;              // getopt_long() 的返回值
    int digit_optind = 0; // 设置短参数类型及是否需要参数

    int option_index = 0;
    // 设置短参数类型及是否需要参数
    const char *optstring = "a";
    /*
    struct option {
             const char * name;  // 参数的名称
             int has_arg; // 是否带参数值，有三种：no_argument， required_argument，optional_argument
             int * flag; // 为空时，函数直接将 val 的数值从getopt_long的返回值返回出去，
                         当非空时，val的值会被赋到 flag 指向的整型数中，而函数返回值为0
             int val; // 用于指定函数找到该选项时的返回值，或者当flag非空时指定flag指向的数据的值
        };
     */
    static struct option long_options[] = {
        {"workdir", required_argument, NULL, 'a'},
        {"mixture", required_argument, NULL, 'b'},
        {"pressure", required_argument, NULL, 'c'},
        {"temperature", required_argument, NULL, 'd'},
        {0, 0, 0, 0} // 添加 {0, 0, 0, 0} 是为了防止输入空值
    };

    while ((opt = getopt_long(argc,
                              argv,
                              optstring,
                              long_options,
                              &option_index)) != -1)
    {

        // printf("opt = %c\n", opt);                           // 命令参数，亦即 -a -b -n -r
        // printf("optarg = %s\n", optarg);                     // 参数内容
        // printf("optind = %d\n", optind);                     // 下一个被处理的下标值
        // printf("argv[optind - 1] = %s\n", argv[optind - 1]); // 参数内容
        // printf("option_index = %d\n", option_index);         // 当前打印参数的下标值
        // printf("\n");
        if (opt == 'a')
            input.insert(std::pair<int, std::string>(0, optarg));
        else if (opt == 'b')
            input.insert(std::pair<int, std::string>(1, optarg));
        else if (opt == 'c')
            input.insert(std::pair<int, std::string>(2, optarg));
        else if (opt == 'd')
            input.insert(std::pair<int, std::string>(3, optarg));
    }
}
