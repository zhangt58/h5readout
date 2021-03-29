#include <cstring>
#include <ctime>

#include <algorithm>

#include "misc.h"
#include "h5readout.h"
#include "processor.h"

ArgumentParser::ArgumentParser()
{
    m_progname = "h5readout";
    m_progdesc = "Readout event data to HDF5 format.";
    m_chunk_dims = new hsize_t[2];
    init_options();
}

ArgumentParser::~ArgumentParser()
{
    delete m_chunk_dims;
    delete m_exclude_types;
    delete m_module_list;
}

std::string ArgumentParser::get_progname()
{
    return m_progname;
}

std::string ArgumentParser::get_progdesc()
{
    return m_progdesc;
}

std::string ArgumentParser::get_default(std::string name)
{
    return m_default_params[name];
}

void ArgumentParser::init_options()
{
    m_parser.add_params({"-o", "--output",
                         "-n", "--events",
                         "-s", "--chunk-size",
                         "-c", "--compress",
                         "--compress-level",
                         "-t", "--source-type",
                         "-m", "--modules",
                         "-v", "--verbose",
                         "-h", "--help",
                         "--exclude",
                         "--version"});
    m_default_params["output"] = "<input>.h5";
    m_default_params["events"] = std::to_string(INT_MAX);
    m_default_params["chunk-size"] = "0x0";
    m_default_params["compress"] = "szip";
    m_default_params["compress-level"] = std::to_string(8);
    m_default_params["controller-type"] = "ddas";
    m_default_params["verbosity"] = std::to_string(0);
    m_default_params["exclude"] = "";
    m_default_params["modules"] = "";
}

std::vector<std::string> *ArgumentParser::get_module_list()
{
    return m_module_list;
}

std::vector<uint16_t> *ArgumentParser::get_exclude_types()
{
    return m_exclude_types;
}

std::string ArgumentParser::get_ctrl_type()
{
    std::transform(m_ctrl_type.begin(), m_ctrl_type.end(), m_ctrl_type.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return m_ctrl_type;
}

void ArgumentParser::parse(int argc, char **argv)
{
    m_parser.parse(argc, argv);

    m_parser(0) >> m_prog;

    // gzip compress level
    m_parser("compress-level", get_default("compress-level")) >> m_gzip_comp_level;

    // compress method
    m_parser({"c", "compress"}, get_default("compress")) >> m_comp_method;

    // chunk dims
    std::stringstream ss(m_parser({"s", "chunk-size"}, get_default("chunk-size")).str());
    int i, idx = 0;
    while (ss >> i)
    {
        m_chunk_dims[idx++] = i;
        if (ss.peek() == 'x')
            ss.ignore();
    }

    // exclude type list
    std::string buf;
    std::stringstream ss1(m_parser("exclude", get_default("exclude")).str());
    while (getline(ss1, buf, ','))
    {
        auto it = ITEM_TYPE.find(buf);
        if (it != ITEM_TYPE.end())
        {
            m_exclude_types->push_back(it->second);
        }
    }

    // module list
    std::stringstream ss2(m_parser("modules", get_default("modules")).str());
    while (getline(ss2, buf, ','))
    {
        auto it = VME_MODULES.find(buf);
        if (it != VME_MODULES.end())
        {
            m_module_list->push_back(*it);
        }
    }

    // max event number
    m_parser({"n", "events"}, get_default("events")) >> m_max_evt;

    // verbose level
    m_parser({"v", "verbose"}, get_default("verbosity")) >> m_verbosity;

    // source type
    m_parser({"t", "controller-type"}, get_default("controller-type")) >> m_ctrl_type;

    // input data source
    m_parser(1) >> m_ifname; // first positional arg

    // output file name
    m_parser({"o", "output"}, get_default("output")) >> m_ofname;
    if (m_ofname == "<input>.h5")
    {
        m_ofname = m_ifname + ".h5";
    }
}

bool ArgumentParser::validate_input_data_source()
{
    char full_datasource[PATH_MAX + 1];
    if (realpath(m_ifname.c_str(), full_datasource) == nullptr)
    {
        return false;
    }
    else
    {
        m_ifname_full = std::string(full_datasource);
        return true;
    }
}

int ArgumentParser::get_verbosity()
{
    return m_verbosity;
}

std::string ArgumentParser::get_version()
{
    return APP_VERSION;
}

bool ArgumentParser::print_version_if_possible()
{
    if (m_parser["version"])
    {
        std::cout << APP_VERSION << std::endl;
        return true;
    }
    return false;
}

bool ArgumentParser::print_help_if_possible()
{
    if (m_parser[{"h", "help"}] || m_ifname.empty())
    {
        print_help();
        std::cout << std::endl;
        print_examples();
        return true;
    }
    return false;
}

int ArgumentParser::get_gzip_compress_level()
{
    return m_gzip_comp_level;
}

std::string ArgumentParser::get_compress_method()
{
    return m_comp_method;
}

hsize_t *ArgumentParser::get_chunk_dims()
{
    return m_chunk_dims;
}

std::string ArgumentParser::get_input_data_source()
{
    return m_ifname;
}

std::string ArgumentParser::get_input_data_source_as_uri()
{
    return "file://" + m_ifname_full;
}

uint64_t ArgumentParser::get_max_evt()
{
    return m_max_evt;
}

std::string ArgumentParser::get_output_filepath()
{
    char outfile[PATH_MAX + 1];
    realpath(m_ofname.c_str(), outfile);
    m_ofname_full.assign(outfile);
    if (is_dir(outfile))
    {
        // valid directory?
        // + auto filename
        char str_tmp[m_ifname_full.length() + 1];
        strcpy(str_tmp, m_ifname_full.c_str());
        m_ofname_full.assign(
            std::string(outfile) + "/" + std::string(basename(str_tmp)) + ".h5");
    }
    return m_ofname_full;
}

void ArgumentParser::print_all_args()
{
    hsize_t *pdim = get_chunk_dims();
    std::cout << "Input data source: " << get_input_data_source() << "\n"
              << "Input data source (URI): " << get_input_data_source_as_uri() << "\n"
              << "Device controller type: " << get_ctrl_type() << "\n"
              << "Output data path: " << get_output_filepath() << "\n"
              << "Max physics events: " << get_max_evt() << "\n"
              << "Chunk dims: " << pdim[0] << "x" << pdim[1] << "\n"
              << "Compress method: " << get_compress_method() << "\n"
              << "Gzip compress level: " << get_gzip_compress_level() << "\n"
              << "Verbose level : " << get_verbosity() << "\n"
              << "Show version? : " << print_version_if_possible() << "\n"
              << "Show help? : " << print_help_if_possible() << "\n";

    std::cout << "Exclude types: ";
    for (std::vector<uint16_t>::iterator it = m_exclude_types->begin(); it != m_exclude_types->end(); ++it)
    {
        std::cout << *it << ",";
    }
    std::cout << std::endl;

    std::cout << "Modules: ";
    for (auto it = m_module_list->begin(); it != m_module_list->end(); ++it)
    {
        std::cout << *it << ",";
    }
    std::cout << std::endl;
}

/**
 * Print out examples of using h5readout.
 *
 */
void ArgumentParser::print_examples()
{
    std::cout << "Examples:"
              << "\n"
              << "  1-Default output h5 filepath:\n"
              << "  $ " << m_prog << " /home/devuser/data.evt # output /home/devuser/data.h5\n"
              << "  2-Default output h5 filepath, to a directory:\n"
              << "  $ " << m_prog << " /home/devuser/data.evt -o /home/devuser/h5data\n"
              << "  3-Output h5 file w/o compression:\n"
              << "  $ " << m_prog << " /home/devuser/data.evt -c none\n"
              << "  4-Output h5 file w/ szip compression:\n"
              << "  $ " << m_prog << " /home/devuser/data.evt -c szip\n"
              << "  5-Skip parsing physics events:\n"
              << "  $ " << m_prog << " data.evt --exclude PHYSICS_EVENT"
              << std::endl;
}

void ArgumentParser::print_help()
{
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
    printf("  %-28s   %s\n", "-c, --compress arg", "Compression method, 'szip' or 'gzip'");
    printf("  %-28s   %s\n", "", "(default: szip)");
    printf("  %-28s   %s\n", "    --compress-level arg", "GZip Compression level (0-9) (default: 8)");
    printf("  %-28s   %s\n", "-t, --controller-type arg", "Type of device controller, 'ddas' (default), 'vme'");
    printf("  %-28s   %s\n", "-m, --modules arg", "A string of VME modules to read, separated by ','");
    printf("  %-28s   %s\n", "", "only required when controller type is not DDAS");
    printf("  %-28s   %s\n", "    --exclude arg", "A string of ringitem types to skip, separated by ','");
    printf("  %-28s   %s\n", "", "e.g. PHYSICS_EVENT (skip events), or PERIODIC_SCALERS,PHYSICS_EVENT");
    printf("  %-28s   %s\n", "", "(skip scalers and events), by default is empty");
    printf("  %-28s   %s\n", "-v, --verbose arg", "Level of verbosity (0-2) (default: 0)");
    printf("  %-28s   %s\n", "", "1: Show Scaler info");
    printf("  %-28s   %s\n", "", "2: Show Scaler and Event info");
    printf("  %-28s   %s\n", "", "3: Show all info");
    printf("  %-28s   %s\n", "    --version", "Show version info");
    printf("  %-28s   %s\n", "-h, --help", "Print this message");
}

bool write_metadata(RunMetaData &metadata, H5::H5File *group)
{
    // run number
    H5::Attribute *number = new H5::Attribute(group->createAttribute(META_DATA_RUN_NUMBER,
                                                                     H5::PredType::NATIVE_INT,
                                                                     H5::DataSpace(H5S_SCALAR)));
    number->write(H5::PredType::NATIVE_INT, &metadata.number);

    // title
    H5::StrType stype(H5::PredType::C_S1, 80);
    H5::Attribute *title = new H5::Attribute(group->createAttribute(META_DATA_TITLE,
                                                                    stype,
                                                                    H5::DataSpace(H5S_SCALAR)));
    title->write(stype, metadata.title);

    // ts
    H5::Attribute *ts0 = new H5::Attribute(group->createAttribute(META_DATA_TIMESTAMP_0,
                                                                  H5::PredType::NATIVE_INT,
                                                                  H5::DataSpace(H5S_SCALAR)));
    ts0->write(H5::PredType::NATIVE_INT, &metadata.ts0);

    H5::Attribute *ts1 = new H5::Attribute(group->createAttribute(META_DATA_TIMESTAMP_1,
                                                                  H5::PredType::NATIVE_INT,
                                                                  H5::DataSpace(H5S_SCALAR)));
    ts1->write(H5::PredType::NATIVE_INT, &metadata.ts1);

    // date
    H5::Attribute *date0 = new H5::Attribute(group->createAttribute(META_DATA_DATETIME_0,
                                                                    stype,
                                                                    H5::DataSpace(H5S_SCALAR)));
    date0->write(stype, metadata.date0);

    H5::Attribute *date1 = new H5::Attribute(group->createAttribute(META_DATA_DATETIME_1,
                                                                    stype,
                                                                    H5::DataSpace(H5S_SCALAR)));
    date1->write(stype, metadata.date1);

    // dt
    H5::Attribute *dt = new H5::Attribute(group->createAttribute(META_DATA_ELAPSEDTIME,
                                                                 H5::PredType::NATIVE_INT,
                                                                 H5::DataSpace(H5S_SCALAR)));
    dt->write(H5::PredType::NATIVE_INT, &metadata.dt);

    // ring format
    H5::StrType stype1(H5::PredType::C_S1, 6);
    H5::Attribute *fmt = new H5::Attribute(group->createAttribute(META_DATA_RING_FORMAT,
                                                                  stype1,
                                                                  H5::DataSpace(H5S_SCALAR)));
    fmt->write(stype1, metadata.fmt);

    // controller type
    H5::Attribute *ctrl_type = new H5::Attribute(group->createAttribute(META_DATA_CTRL_TYPE,
                                                                        stype1,
                                                                        H5::DataSpace(H5S_SCALAR)));
    ctrl_type->write(stype1, metadata.ctrl_type);

    // total events
    H5::Attribute *n_events = new H5::Attribute(group->createAttribute(META_DATA_TOTAL_EVENTS,
                                                                       H5::PredType::NATIVE_LONG,
                                                                       H5::DataSpace(H5S_SCALAR)));
    n_events->write(H5::PredType::NATIVE_LONG, &metadata.n_events);

    // total fragments
    H5::Attribute *n_frags = new H5::Attribute(group->createAttribute(META_DATA_TOTAL_FRAGMENTS,
                                                                      H5::PredType::NATIVE_LONG,
                                                                      H5::DataSpace(H5S_SCALAR)));
    n_frags->write(H5::PredType::NATIVE_LONG, &metadata.n_frags);

    delete number;
    delete title;
    delete ts0;
    delete ts1;
    delete date0;
    delete date1;
    delete dt;
    delete fmt;
    delete n_events;
    delete n_frags;

    return true;
}

bool write_fragdata(std::vector<FragmentData> *pfragdata, H5::Group *group)
{
    hsize_t dsize = pfragdata->size();
    if (dsize == 0)
    {
        return false;
    }

    // create mem data type for FragmentData
    const H5::CompType frag_dtype(sizeof(FragmentData));
    frag_dtype.insertMember(FRAGMENT_DATA_FRAGMENT_ID, HOFFSET(FragmentData, fragment_id), H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_EVENT_ID, HOFFSET(FragmentData, event_id), H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_TIMESTAMP, HOFFSET(FragmentData, timestamp), H5::PredType::NATIVE_DOUBLE);
    frag_dtype.insertMember(FRAGMENT_DATA_COARSE_TIME, HOFFSET(FragmentData, coarse_time), H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_ENERGY, HOFFSET(FragmentData, energy), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TRACE_LENGTH, HOFFSET(FragmentData, trace_length), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CRATE_ID, HOFFSET(FragmentData, crate_id), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_SLOT_ID, HOFFSET(FragmentData, slot_id), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_ID, HOFFSET(FragmentData, channel_id), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_FINISH_CODE, HOFFSET(FragmentData, finish_code), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_OVER_UNDER_FLOW, HOFFSET(FragmentData, adc_over_underflow), H5::PredType::NATIVE_SHORT);
    frag_dtype.insertMember(FRAGMENT_DATA_CFD_FAIL_BIT, HOFFSET(FragmentData, cfd_fail_bit), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_OVERFLOW_CODE, HOFFSET(FragmentData, overflow_code), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_MODMSPS, HOFFSET(FragmentData, adc_frequency), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_RESOLUTION, HOFFSET(FragmentData, adc_resolution), H5::PredType::NATIVE_INT);

    hsize_t dim[] = {dsize};
    H5::DataSpace *dspace = new H5::DataSpace(FRAGMENT_DATA_RANK, dim);

    // create a dataset under the new group
    H5::DataSet *dset = new H5::DataSet(group->createDataSet(FRAGMENTS_DSET_NAME, frag_dtype, *dspace));

    // write dataset: fragments
    dset->write(pfragdata->data(), frag_dtype);

    // clean up
    delete pfragdata;
    delete dset;
    delete dspace;
    //
    return true;
}

bool write_tracedata(std::vector<uint16_t> *ptracedata, H5::Group *group, uint64_t &frag_cnt,
                     hsize_t *chunk_dims, std::string &comp_meth, int &gfactor)
{
    hsize_t dsize = ptracedata->size();
    if (dsize == 0)
    {
        return false;
    }

    // create 2d array for trace data
    int dim0, dim1;
    int trace_length = dsize / frag_cnt; // frag_cnt (nrows), trace_length (ncolumns)
    dim0 = frag_cnt;
    dim1 = trace_length;
    int nrows_sub = T_NROWS_PER_WRITE; // how many rows write per time
    int current_row_id = 0;            // starting at the first fragment
    uint16_t tracedata_subarr[nrows_sub][dim1];

    // auto chunk if 0x0
    if (chunk_dims[0] == 0)
    {
        chunk_dims[1] = trace_length;
        chunk_dims[0] = 1.0 * 1024 * 1024 / sizeof(uint16_t) / trace_length; // 1MB
    }

    // create a new dataset under the defined group, TraceData
    H5::IntType trace_dtype(H5::PredType::NATIVE_SHORT);
    trace_dtype.setOrder(H5T_ORDER_LE);

    // extensible dataset for Traces
    hsize_t trace_dims[2]; // initial dset shape
    trace_dims[0] = nrows_sub;
    trace_dims[1] = dim1;
    // hsize_t max_trace_dims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t max_trace_dims[2] = {(hsize_t)frag_cnt, (hsize_t)trace_length};
    H5::DataSpace trace_dspace(TRACE_DATA_RANK, trace_dims, max_trace_dims);

    // modify dataset creation properties, e.g. enable chunking
    H5::DSetCreatPropList cprops;
    cprops.setChunk(TRACE_DATA_RANK, chunk_dims);

    if (comp_meth == "szip")
    {
        cprops.setSzip(H5_SZIP_NN_OPTION_MASK, 16);
    }
    else if (comp_meth == "gzip")
    {
        cprops.setDeflate(gfactor);
    }

    // create Traces dataset
    H5::DataSet trace_dset = group->createDataSet(TRACES_DSET_NAME, trace_dtype, trace_dspace, cprops);

    H5::DataSpace fspace;
    hsize_t size[2];
    size[0] = 0;
    size[1] = dim1;
    hsize_t offset[2];
    offset[1] = 0;
    hsize_t trace_dims1[2];
    trace_dims1[1] = dim1;

    while (current_row_id < dim0)
    {
        // prep data
        for (int i = 0; i < nrows_sub; i++)
        {
            for (int j = 0; j < dim1; j++)
            {
                tracedata_subarr[i][j] = (*ptracedata)[(i + current_row_id) * dim1 + j];
            }
        }

        // extend size along dim0
        size[0] += nrows_sub;

        // extend the dataset
        trace_dset.extend(size);

        // select a hyperslab
        fspace = trace_dset.getSpace();
        offset[0] = current_row_id;
        trace_dims1[0] = nrows_sub;
        fspace.selectHyperslab(H5S_SELECT_SET, trace_dims1, offset);

        trace_dset.write(tracedata_subarr, trace_dtype, trace_dspace, fspace);

        current_row_id += nrows_sub;
    }

    // clean up
    delete ptracedata;

    //
    return true;
}

bool write_scalerdata(std::vector<uint32_t> *pscalerdata, H5::Group *group,
                      std::vector<uint32_t> *pscalerlen, std::vector<time_t> *pscalerts,
                      std::string &stype)
{
    bool r;
    if (stype == "DDAS")
    {
        r = _write_ddas_scalerdata(pscalerdata, group, pscalerlen, pscalerts);
    }
    else if (stype == "VME")
    {
        r = _write_vme_scalerdata(pscalerdata, group, pscalerlen, pscalerts);
    }
    return r;
}

bool _write_ddas_scalerdata(std::vector<uint32_t> *pscalerdata, H5::Group *group,
                            std::vector<uint32_t> *pscalerlen, std::vector<time_t> *pscalerts)
{
    hsize_t dsize = pscalerts->size();
    if (dsize == 0)
    {
        return false;
    }

    // create mem data type for ScalerData
    const H5::CompType scaler_dtype(sizeof(DDASScalerData));
    scaler_dtype.insertMember(SCALER_TIMESTAMP, HOFFSET(DDASScalerData, timestamp), H5::PredType::NATIVE_LONG);
    scaler_dtype.insertMember(SCALER_DATETIME, HOFFSET(DDASScalerData, datetime), H5::StrType(H5::PredType::C_S1, 80));
    scaler_dtype.insertMember(SCALER_CRATE_ID, HOFFSET(DDASScalerData, crate_id), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_SLOT_ID, HOFFSET(DDASScalerData, slot_id), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH00, HOFFSET(DDASScalerData, inc_raw_ch00), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH00, HOFFSET(DDASScalerData, inc_val_ch00), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH01, HOFFSET(DDASScalerData, inc_raw_ch01), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH01, HOFFSET(DDASScalerData, inc_val_ch01), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH02, HOFFSET(DDASScalerData, inc_raw_ch02), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH02, HOFFSET(DDASScalerData, inc_val_ch02), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH03, HOFFSET(DDASScalerData, inc_raw_ch03), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH03, HOFFSET(DDASScalerData, inc_val_ch03), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH04, HOFFSET(DDASScalerData, inc_raw_ch04), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH04, HOFFSET(DDASScalerData, inc_val_ch04), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH05, HOFFSET(DDASScalerData, inc_raw_ch05), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH05, HOFFSET(DDASScalerData, inc_val_ch05), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH06, HOFFSET(DDASScalerData, inc_raw_ch06), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH06, HOFFSET(DDASScalerData, inc_val_ch06), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH07, HOFFSET(DDASScalerData, inc_raw_ch07), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH07, HOFFSET(DDASScalerData, inc_val_ch07), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH08, HOFFSET(DDASScalerData, inc_raw_ch08), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH08, HOFFSET(DDASScalerData, inc_val_ch08), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH09, HOFFSET(DDASScalerData, inc_raw_ch09), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH09, HOFFSET(DDASScalerData, inc_val_ch09), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH10, HOFFSET(DDASScalerData, inc_raw_ch10), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH10, HOFFSET(DDASScalerData, inc_val_ch10), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH11, HOFFSET(DDASScalerData, inc_raw_ch11), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH11, HOFFSET(DDASScalerData, inc_val_ch11), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH12, HOFFSET(DDASScalerData, inc_raw_ch12), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH12, HOFFSET(DDASScalerData, inc_val_ch12), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH13, HOFFSET(DDASScalerData, inc_raw_ch13), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH13, HOFFSET(DDASScalerData, inc_val_ch13), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH14, HOFFSET(DDASScalerData, inc_raw_ch14), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH14, HOFFSET(DDASScalerData, inc_val_ch14), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH15, HOFFSET(DDASScalerData, inc_raw_ch15), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH15, HOFFSET(DDASScalerData, inc_val_ch15), H5::PredType::NATIVE_INT);

    int seg_num = 33;

    // create 2d array for scaler data
    int s_dim0 = dsize, s_dim1 = (*pscalerlen)[0]; // assume all scaler of the same length
    uint32_t max_slot_num = 1 + s_dim1 / seg_num;  // max slot number (ddas only)
    uint32_t slot_num = 2;                         // first slot number (ddas only)
    time_t ts_tmp;
    char date_tmp[80];
    int ind0;

    // process scalerdata
    std::vector<DDASScalerData> *pscaler = new std::vector<DDASScalerData>();
    for (int i = 0; i < s_dim0; i++)
    {
        ts_tmp = (*pscalerts)[i];
        std::strftime(date_tmp, sizeof(date_tmp), "%c", std::localtime(&ts_tmp));
        slot_num = 2;
        while (slot_num <= max_slot_num)
        {
            // create a new DDASScalerData every 33, 33, ...
            // crate_id + 16 * 2 channel readings
            ind0 = (slot_num - 2) * seg_num + i * s_dim1;
            DDASScalerData i_scaler_data;
            i_scaler_data.slot_id = slot_num;
            i_scaler_data.crate_id = (*pscalerdata)[ind0];
            i_scaler_data.timestamp = ts_tmp;
            strcpy(i_scaler_data.datetime, date_tmp);
            i_scaler_data.inc_raw_ch00 = (*pscalerdata)[ind0 + 1];
            i_scaler_data.inc_val_ch00 = (*pscalerdata)[ind0 + 2];
            i_scaler_data.inc_raw_ch01 = (*pscalerdata)[ind0 + 3];
            i_scaler_data.inc_val_ch01 = (*pscalerdata)[ind0 + 4];
            i_scaler_data.inc_raw_ch02 = (*pscalerdata)[ind0 + 5];
            i_scaler_data.inc_val_ch02 = (*pscalerdata)[ind0 + 6];
            i_scaler_data.inc_raw_ch03 = (*pscalerdata)[ind0 + 7];
            i_scaler_data.inc_val_ch03 = (*pscalerdata)[ind0 + 8];
            i_scaler_data.inc_raw_ch04 = (*pscalerdata)[ind0 + 9];
            i_scaler_data.inc_val_ch04 = (*pscalerdata)[ind0 + 10];
            i_scaler_data.inc_raw_ch05 = (*pscalerdata)[ind0 + 11];
            i_scaler_data.inc_val_ch05 = (*pscalerdata)[ind0 + 12];
            i_scaler_data.inc_raw_ch06 = (*pscalerdata)[ind0 + 13];
            i_scaler_data.inc_val_ch06 = (*pscalerdata)[ind0 + 14];
            i_scaler_data.inc_raw_ch07 = (*pscalerdata)[ind0 + 15];
            i_scaler_data.inc_val_ch07 = (*pscalerdata)[ind0 + 16];
            i_scaler_data.inc_raw_ch08 = (*pscalerdata)[ind0 + 17];
            i_scaler_data.inc_val_ch08 = (*pscalerdata)[ind0 + 18];
            i_scaler_data.inc_raw_ch09 = (*pscalerdata)[ind0 + 19];
            i_scaler_data.inc_val_ch09 = (*pscalerdata)[ind0 + 20];
            i_scaler_data.inc_raw_ch10 = (*pscalerdata)[ind0 + 21];
            i_scaler_data.inc_val_ch10 = (*pscalerdata)[ind0 + 22];
            i_scaler_data.inc_raw_ch11 = (*pscalerdata)[ind0 + 23];
            i_scaler_data.inc_val_ch11 = (*pscalerdata)[ind0 + 24];
            i_scaler_data.inc_raw_ch12 = (*pscalerdata)[ind0 + 25];
            i_scaler_data.inc_val_ch12 = (*pscalerdata)[ind0 + 26];
            i_scaler_data.inc_raw_ch13 = (*pscalerdata)[ind0 + 27];
            i_scaler_data.inc_val_ch13 = (*pscalerdata)[ind0 + 28];
            i_scaler_data.inc_raw_ch14 = (*pscalerdata)[ind0 + 29];
            i_scaler_data.inc_val_ch14 = (*pscalerdata)[ind0 + 30];
            i_scaler_data.inc_raw_ch15 = (*pscalerdata)[ind0 + 31];
            i_scaler_data.inc_val_ch15 = (*pscalerdata)[ind0 + 32];
            ++slot_num;
            pscaler->push_back(i_scaler_data);
        }
    }

    struct sort_key
    {
        inline bool operator()(const DDASScalerData &left, const DDASScalerData &right)
        {
            if (left.slot_id < right.slot_id)
                return true;
            if (left.slot_id > right.slot_id)
                return false;
            if (left.timestamp < right.timestamp)
                return true;
            if (left.timestamp > right.timestamp)
                return false;
            return false;
        }
    };
    std::sort(pscaler->begin(), pscaler->end(), sort_key());

    // create a dataspace for scaler data
    hsize_t s_dim[] = {pscaler->size()};
    H5::DataSpace *s_dspace = new H5::DataSpace(SCALER_DATA_RANK, s_dim);

    // create a dataset under the new group
    H5::DataSet *s_dset = new H5::DataSet(group->createDataSet(SCALERS_DSET_NAME, scaler_dtype, *s_dspace));

    // write dataset: scalerdata
    s_dset->write(pscaler->data(), scaler_dtype);

    // clean up
    delete pscalerdata;
    delete pscalerts;
    delete pscalerlen;
    delete pscaler;

    delete s_dspace;
    delete s_dset;

    //
    return true;
}

bool _write_vme_scalerdata(std::vector<uint32_t> *pscalerdata, H5::Group *group,
                           std::vector<uint32_t> *pscalerlen, std::vector<time_t> *pscalerts)
{

    hsize_t dsize = pscalerts->size();
    if (dsize == 0)
    {
        return false;
    }

    const H5::CompType scaler_dtype(sizeof(ScalerData));
    scaler_dtype.insertMember(SCALER_TIMESTAMP, HOFFSET(ScalerData, timestamp), H5::PredType::NATIVE_LONG);
    scaler_dtype.insertMember(SCALER_DATETIME, HOFFSET(ScalerData, datetime), H5::StrType(H5::PredType::C_S1, 80));
    scaler_dtype.insertMember(SCALER_SLOT_ID, HOFFSET(ScalerData, slot_id), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH00, HOFFSET(ScalerData, inc_raw_ch00), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH00, HOFFSET(ScalerData, inc_val_ch00), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH01, HOFFSET(ScalerData, inc_raw_ch01), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH01, HOFFSET(ScalerData, inc_val_ch01), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH02, HOFFSET(ScalerData, inc_raw_ch02), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH02, HOFFSET(ScalerData, inc_val_ch02), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH03, HOFFSET(ScalerData, inc_raw_ch03), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH03, HOFFSET(ScalerData, inc_val_ch03), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH04, HOFFSET(ScalerData, inc_raw_ch04), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH04, HOFFSET(ScalerData, inc_val_ch04), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH05, HOFFSET(ScalerData, inc_raw_ch05), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH05, HOFFSET(ScalerData, inc_val_ch05), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH06, HOFFSET(ScalerData, inc_raw_ch06), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH06, HOFFSET(ScalerData, inc_val_ch06), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH07, HOFFSET(ScalerData, inc_raw_ch07), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH07, HOFFSET(ScalerData, inc_val_ch07), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH08, HOFFSET(ScalerData, inc_raw_ch08), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH08, HOFFSET(ScalerData, inc_val_ch08), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH09, HOFFSET(ScalerData, inc_raw_ch09), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH09, HOFFSET(ScalerData, inc_val_ch09), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH10, HOFFSET(ScalerData, inc_raw_ch10), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH10, HOFFSET(ScalerData, inc_val_ch10), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH11, HOFFSET(ScalerData, inc_raw_ch11), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH11, HOFFSET(ScalerData, inc_val_ch11), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH12, HOFFSET(ScalerData, inc_raw_ch12), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH12, HOFFSET(ScalerData, inc_val_ch12), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH13, HOFFSET(ScalerData, inc_raw_ch13), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH13, HOFFSET(ScalerData, inc_val_ch13), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH14, HOFFSET(ScalerData, inc_raw_ch14), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH14, HOFFSET(ScalerData, inc_val_ch14), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_RAW_CH15, HOFFSET(ScalerData, inc_raw_ch15), H5::PredType::NATIVE_INT);
    scaler_dtype.insertMember(SCALER_VAL_CH15, HOFFSET(ScalerData, inc_val_ch15), H5::PredType::NATIVE_INT);

    int seg_num = 32;
    // create 2d array for scaler data
    int s_dim0 = dsize, s_dim1 = (*pscalerlen)[0]; // assume all scaler of the same length
    uint32_t max_slot_num = 1 + s_dim1 / seg_num;  // max slot number (ddas only)
    uint32_t slot_num = 2;                         // first slot number (ddas only)
    time_t ts_tmp;
    char date_tmp[80];
    int ind0;

    // process scalerdata
    std::vector<ScalerData> *pscaler = new std::vector<ScalerData>();
    for (int i = 0; i < s_dim0; i++)
    {
        ts_tmp = (*pscalerts)[i];
        std::strftime(date_tmp, sizeof(date_tmp), "%c", std::localtime(&ts_tmp));
        slot_num = 2;
        while (slot_num <= max_slot_num)
        {
            ind0 = (slot_num - 2) * seg_num + i * s_dim1 - 1; // no crate_id
            ScalerData i_scaler_data;
            i_scaler_data.timestamp = ts_tmp;
            strcpy(i_scaler_data.datetime, date_tmp);
            i_scaler_data.slot_id = slot_num;
            i_scaler_data.inc_raw_ch00 = (*pscalerdata)[ind0 + 1];
            i_scaler_data.inc_val_ch00 = (*pscalerdata)[ind0 + 2];
            i_scaler_data.inc_raw_ch01 = (*pscalerdata)[ind0 + 3];
            i_scaler_data.inc_val_ch01 = (*pscalerdata)[ind0 + 4];
            i_scaler_data.inc_raw_ch02 = (*pscalerdata)[ind0 + 5];
            i_scaler_data.inc_val_ch02 = (*pscalerdata)[ind0 + 6];
            i_scaler_data.inc_raw_ch03 = (*pscalerdata)[ind0 + 7];
            i_scaler_data.inc_val_ch03 = (*pscalerdata)[ind0 + 8];
            i_scaler_data.inc_raw_ch04 = (*pscalerdata)[ind0 + 9];
            i_scaler_data.inc_val_ch04 = (*pscalerdata)[ind0 + 10];
            i_scaler_data.inc_raw_ch05 = (*pscalerdata)[ind0 + 11];
            i_scaler_data.inc_val_ch05 = (*pscalerdata)[ind0 + 12];
            i_scaler_data.inc_raw_ch06 = (*pscalerdata)[ind0 + 13];
            i_scaler_data.inc_val_ch06 = (*pscalerdata)[ind0 + 14];
            i_scaler_data.inc_raw_ch07 = (*pscalerdata)[ind0 + 15];
            i_scaler_data.inc_val_ch07 = (*pscalerdata)[ind0 + 16];
            i_scaler_data.inc_raw_ch08 = (*pscalerdata)[ind0 + 17];
            i_scaler_data.inc_val_ch08 = (*pscalerdata)[ind0 + 18];
            i_scaler_data.inc_raw_ch09 = (*pscalerdata)[ind0 + 19];
            i_scaler_data.inc_val_ch09 = (*pscalerdata)[ind0 + 20];
            i_scaler_data.inc_raw_ch10 = (*pscalerdata)[ind0 + 21];
            i_scaler_data.inc_val_ch10 = (*pscalerdata)[ind0 + 22];
            i_scaler_data.inc_raw_ch11 = (*pscalerdata)[ind0 + 23];
            i_scaler_data.inc_val_ch11 = (*pscalerdata)[ind0 + 24];
            i_scaler_data.inc_raw_ch12 = (*pscalerdata)[ind0 + 25];
            i_scaler_data.inc_val_ch12 = (*pscalerdata)[ind0 + 26];
            i_scaler_data.inc_raw_ch13 = (*pscalerdata)[ind0 + 27];
            i_scaler_data.inc_val_ch13 = (*pscalerdata)[ind0 + 28];
            i_scaler_data.inc_raw_ch14 = (*pscalerdata)[ind0 + 29];
            i_scaler_data.inc_val_ch14 = (*pscalerdata)[ind0 + 30];
            i_scaler_data.inc_raw_ch15 = (*pscalerdata)[ind0 + 31];
            i_scaler_data.inc_val_ch15 = (*pscalerdata)[ind0 + 32];
            ++slot_num;
            pscaler->push_back(i_scaler_data);
        }
    }
    struct sort_key
    {
        inline bool operator()(const ScalerData &left, const ScalerData &right)
        {
            if (left.slot_id < right.slot_id)
                return true;
            if (left.slot_id > right.slot_id)
                return false;
            if (left.timestamp < right.timestamp)
                return true;
            if (left.timestamp > right.timestamp)
                return false;
            return false;
        }
    };
    std::sort(pscaler->begin(), pscaler->end(), sort_key());

    // create a dataspace for scaler data
    hsize_t s_dim[] = {pscaler->size()};
    H5::DataSpace *s_dspace = new H5::DataSpace(SCALER_DATA_RANK, s_dim);

    // create a dataset under the new group
    H5::DataSet *s_dset = new H5::DataSet(group->createDataSet(SCALERS_DSET_NAME, scaler_dtype, *s_dspace));

    // write dataset: scalerdata
    s_dset->write(pscaler->data(), scaler_dtype);

    // clean up
    delete pscalerdata;
    delete pscalerts;
    delete pscalerlen;
    delete pscaler;

    delete s_dspace;
    delete s_dset;

    //
    return true;
}

bool write_fragdata_vme(std::vector<uint64_t> *pfragdata, H5::Group *group,
                        uint64_t &frag_cnt, std::string &comp_meth, int &gfactor)
{
    hsize_t dsize = pfragdata->size();
    if (dsize == 0)
    {
        return false;
    }

    // create 2d array for trace data
    // trace here means:
    //   event_id, geo, crate_id, n_channels, (all values ch0-32), (ov), (un)
    int dim0, dim1;
    int trace_length = dsize / frag_cnt; // frag_cnt (nrows), trace_length (ncolumns)
    dim0 = frag_cnt;
    dim1 = trace_length;
    int nrows_sub = T_NROWS_PER_WRITE; // how many rows write per time
    int current_row_id = 0;            // starting at the first fragment
    uint64_t tracedata_subarr[nrows_sub][dim1];

    hsize_t chunk_dims[2];
    chunk_dims[1] = trace_length;
    chunk_dims[0] = 1.0 * 1024 * 1024 / sizeof(uint64_t) / trace_length; // 1MB

    // create a new dataset under the defined group, TraceData
    H5::IntType trace_dtype(H5::PredType::NATIVE_LONG);
    trace_dtype.setOrder(H5T_ORDER_LE);

    // extensible dataset for Traces
    hsize_t trace_dims[2]; // initial dset shape
    trace_dims[0] = nrows_sub;
    trace_dims[1] = dim1;
    hsize_t max_trace_dims[2] = {(hsize_t)frag_cnt, (hsize_t)trace_length};
    H5::DataSpace trace_dspace(TRACE_DATA_RANK, trace_dims, max_trace_dims);

    // modify dataset creation properties, e.g. enable chunking
    H5::DSetCreatPropList cprops;
    cprops.setChunk(TRACE_DATA_RANK, chunk_dims);

    if (comp_meth == "szip")
    {
        cprops.setSzip(H5_SZIP_NN_OPTION_MASK, 16);
    }
    else if (comp_meth == "gzip")
    {
        cprops.setDeflate(gfactor);
    }

    // create Traces dataset
    H5::DataSet trace_dset = group->createDataSet(FRAGMENTS_DSET_NAME, trace_dtype, trace_dspace, cprops);

    H5::DataSpace fspace;
    hsize_t size[2];
    size[0] = 0;
    size[1] = dim1;
    hsize_t offset[2];
    offset[1] = 0;
    hsize_t trace_dims1[2];
    trace_dims1[1] = dim1;

    while (current_row_id < dim0)
    {
        // prep data
        for (int i = 0; i < nrows_sub; i++)
        {
            for (int j = 0; j < dim1; j++)
            {
                tracedata_subarr[i][j] = (*pfragdata)[(i + current_row_id) * dim1 + j];
            }
        }

        // extend size along dim0
        size[0] += nrows_sub;

        // extend the dataset
        trace_dset.extend(size);

        // select a hyperslab
        fspace = trace_dset.getSpace();
        offset[0] = current_row_id;
        trace_dims1[0] = nrows_sub;
        fspace.selectHyperslab(H5S_SELECT_SET, trace_dims1, offset);

        trace_dset.write(tracedata_subarr, trace_dtype, trace_dspace, fspace);

        current_row_id += nrows_sub;
    }

    // clean up
    delete pfragdata;

    // write attribute for column names
    // event_id, geo, crate_id, n_channels, (all values ch0-32), (ov), (un)
    int nch = (trace_length - 4) / 3;
    std::string column_names = "Event_ID,GEO,Crate_ID,N_Channels,";
    for (int i = 0; i < nch; ++i)
    {
        column_names.append("Value" + std::to_string(i) + ",");
    }
    for (int i = 0; i < nch; ++i)
    {
        column_names.append("Overflow_Code" + std::to_string(i) + ",");
    }
    for (int i = 0; i < nch - 1; ++i)
    {
        column_names.append("Threshold_Code" + std::to_string(i) + ",");
    }
    column_names.append("Threshold_Code" + std::to_string(nch - 1));
    std::cout << column_names << std::endl;
    H5::StrType stype(H5::PredType::C_S1, column_names.length() + 1);
    H5::Attribute *attr = new H5::Attribute(
        trace_dset.createAttribute("Column Names", stype,
                                   H5::DataSpace(H5S_SCALAR)));
    attr->write(stype, column_names);

    delete attr;

    // // create mem data type for FragmentData_V785
    // const H5::CompType frag_dtype(sizeof(FragmentData_V785));
    // auto int16_type = H5::PredType::NATIVE_INT16;
    // auto vint16_type = H5::VarLenType(&int16_type);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_EVENT_ID, HOFFSET(FragmentData_V785, event_id), H5::PredType::NATIVE_LONG);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_GEO, HOFFSET(FragmentData_V785, geo), H5::PredType::NATIVE_INT16);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_CRATE_ID, HOFFSET(FragmentData_V785, crate_id), H5::PredType::NATIVE_INT16);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_TOTAL_CHANNELS, HOFFSET(FragmentData_V785, n_channels), H5::PredType::NATIVE_INT16);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_VALUE, HOFFSET(FragmentData_V785, val_channel_handle), vint16_type);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_OVERFLOW, HOFFSET(FragmentData_V785, ov_channel_handle), vint16_type);
    // frag_dtype.insertMember(V785_FRAGMENT_DATA_THRESHOLD, HOFFSET(FragmentData_V785, un_channel_handle), vint16_type);

    // hsize_t dim[] = {dsize};
    // H5::DataSpace *dspace = new H5::DataSpace(FRAGMENT_DATA_RANK, dim);

    // // create a dataset under the new group
    // H5::DataSet *dset = new H5::DataSet(group->createDataSet(FRAGMENTS_DSET_NAME, frag_dtype, *dspace));

    // for (auto it = pfragdata->begin(); it != pfragdata->end(); ++it)
    // {
    //     (*it).val_channel.assign(static_cast<uint16_t *>((*it).val_channel_handle.p),
    //                              static_cast<uint16_t *>((*it).val_channel_handle.p) + (*it).val_channel_handle.len);

    //     (*it).ov_channel.assign(static_cast<uint16_t *>((*it).ov_channel_handle.p),
    //                             static_cast<uint16_t *>((*it).ov_channel_handle.p) + (*it).ov_channel_handle.len);

    //     (*it).un_channel.assign(static_cast<uint16_t *>((*it).un_channel_handle.p),
    //                             static_cast<uint16_t *>((*it).un_channel_handle.p) + (*it).un_channel_handle.len);
    // }

    // write dataset: fragments
    // dset->write(pfragdata->data(), frag_dtype);

    // clean up
    // delete pfragdata;
    // delete dset;
    // delete dspace;
    //

    return true;
}
