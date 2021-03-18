#include <cstring>

#include "misc.h"
#include "h5readout.h"


ArgumentParser::ArgumentParser() {
    m_progname = "h5readout";
    m_progdesc = "Readout DDAS event data to HDF5 format.";
    m_chunk_dims = new hsize_t[2];
    init_options();
}

ArgumentParser::~ArgumentParser() {
    delete m_chunk_dims;
}

std::string ArgumentParser::get_progname() {
    return m_progname;
}

std::string ArgumentParser::get_progdesc() {
    return m_progdesc;
}

std::string ArgumentParser::get_default(std::string name) {
    return m_default_params[name];
}

void ArgumentParser::init_options() {
    m_parser.add_params({
            "-o", "--output",
            "-n", "--events",
            "-s", "--chunk-size",
            "-c", "--compress",
            "--compress-level",
            "-v", "--verbose",
            "-h", "--help",
            "--version"});
    m_default_params["output"] = "<input>.h5";
    m_default_params["events"] = std::to_string(INT_MAX);
    m_default_params["chunk-size"] = "0x0";
    m_default_params["compress"] = "gzip";
    m_default_params["compress-level"] = std::to_string(8);
}

void ArgumentParser::parse(int argc, char** argv) {
    m_parser.parse(argc, argv);

    m_parser(0) >> m_prog;

    // gzip compress level
    m_parser("compress-level", get_default("compress-level")) >> m_gzip_comp_level;

    // compress method
    m_parser({"c", "compress"}, get_default("compress")) >> m_comp_method;

    // chunk dims
    std::stringstream ss(m_parser({"s", "chunk-size"}, get_default("chunk-size")).str());
    int i, idx = 0;
    while (ss >> i) {
        m_chunk_dims[idx++] = i;
        if (ss.peek() == 'x') ss.ignore();
    }

    // max event number
    m_parser({"n", "events"}, get_default("events")) >> m_max_evt;

    // verbose?
    m_verbose = m_parser[{"v", "verbose"}];

    // input data source
    m_parser(1) >> m_ifname; // first positional arg

    // output file name
    m_parser({"o", "output"}, get_default("output")) >> m_ofname;
    if (m_ofname == "<input>.h5") {
        m_ofname = m_ifname + ".h5";
    }

}

bool ArgumentParser::validate_input_data_source() {
    char full_datasource[PATH_MAX + 1];
    if (realpath(m_ifname.c_str(), full_datasource) == nullptr) {
        return false;
    } else {
        m_ifname_full = std::string(full_datasource);
        return true;
    }
}

bool ArgumentParser::is_verbose() {
    return m_verbose;
}

std::string ArgumentParser::get_version() {
    return VERSION;
}

bool ArgumentParser::print_version_if_possible() {
    if (m_parser["version"]) {
        std::cout << VERSION << std::endl;
        return true;
    }
    return false;
}

bool ArgumentParser::print_help_if_possible() {
    if (m_parser[{"h", "help"}] || m_ifname.empty()) {
        print_help();
        std::cout << std::endl;
        print_examples();
        return true;
    }
    return false;
}

int ArgumentParser::get_gzip_compress_level() {
    return m_gzip_comp_level;
}

std::string ArgumentParser::get_compress_method() {
    return m_comp_method;
}

hsize_t* ArgumentParser::get_chunk_dims() {
    return m_chunk_dims;
}

std::string ArgumentParser::get_input_data_source() {
    return m_ifname;
}

std::string ArgumentParser::get_input_data_source_as_uri() {
    return "file://" + m_ifname_full;
}

uint64_t ArgumentParser::get_max_evt() {
    return m_max_evt;
}

std::string ArgumentParser::get_output_filepath() {
    char outfile[PATH_MAX + 1];
    realpath(m_ofname.c_str(), outfile);
    m_ofname_full.assign(outfile);
    if (is_dir(outfile)) {
        // valid directory?
        // + auto filename
        char str_tmp[m_ifname_full.length() + 1];
        strcpy(str_tmp, m_ifname_full.c_str());
        m_ofname_full.assign(
                std::string(outfile) + "/" + std::string(basename(str_tmp)) + ".h5");
    }
    return m_ofname_full;
}

void ArgumentParser::print_all_args() {
    hsize_t* pdim = get_chunk_dims();
    std::cout << "Input data source: " << get_input_data_source() << "\n"
              << "Input data source (URI): " << get_input_data_source_as_uri() << "\n"
              << "Output data path: " << get_output_filepath() << "\n"
              << "Max physics events: " << get_max_evt() << "\n"
              << "Chunk dims: " << pdim[0] << "x" << pdim[1] << "\n"
              << "Compress method: " << get_compress_method() << "\n"
              << "Gzip compress level: " << get_gzip_compress_level() << "\n"
              << "Verbose on? : "  << is_verbose() << "\n"
              << "Show version? : " << print_version_if_possible() << "\n"
              << "Show help? : " << print_help_if_possible() << "\n";
}

/**
 * Print out examples of using h5readout.
 *
 */
void ArgumentParser::print_examples() {
    std::cout << "Examples:" << "\n"
              << "  # Default output h5 filepath:\n"
              << "  " << m_prog << " /home/devuser/data.evt # output /home/devuser/data.h5\n"
              << "  # Default output h5 filepath, to a directory:\n"
              << "  " << m_prog << " /home/devuser/data.evt -o /home/devuser/h5data\n"
              << "  # Output h5 file w/o compression:\n"
              << "  " << m_prog << " /home/devuser/data.evt -c none\n"
              << "  # Output h5 file w/ szip compression:\n"
              << "  " << m_prog << " /home/devuser/data.evt -c szip"
              << std::endl;
}

void ArgumentParser::print_help() {
    printf("%s\n\n", m_progdesc.c_str());
    printf("Usage:\n");
    printf("  %s [OPTION...] INPUT\n\n", m_prog.c_str());
    printf("  %-28s   %s\n", "INPUT", "URI for input event data source");
    printf("  %-28s   %s\n", "-o, --output arg", "File path for HDF5 data, if assigned with");
    printf("  %-28s   %s\n", "", "a valid directory, apply default filename");
    printf("  %-28s   %s\n", "", "under that directory (default: <INPUT>.h5)");
    printf("  %-28s   %s\n", "-n, --events arg", "Maximum number of physics events to readout");
    printf("  %-28s   %s\n", "", std::string("(default: " + std::to_string(INT_MAX) + ")").c_str());
    printf("  %-28s   %s\n", "-s, --chunk-size arg", "Chunk size MxN for HDF5 data (default: 0x0)");
    printf("  %-28s   %s\n", "-c, --compress arg", "Compression method, 'szip' or 'gzip' (default: gzip)");
    printf("  %-28s   %s\n", "    --compress-level arg", "GZip Compression level(0-9) (default: 8)");
    printf("  %-28s   %s\n", "-v, --verbose", "Show verbose message");
    printf("  %-28s   %s\n", "    --version", "Show version info");
    printf("  %-28s   %s\n", "-h, --help", "Print this message");
}
