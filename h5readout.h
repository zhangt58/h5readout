#pragma once

#include <cstdint>
#include <linux/limits.h> // PATH_MAX

#include <vector>
#include <string>
#include <iostream>

#include "argh.h"    // args parser
#include "H5Cpp.h"   // HDF5
#include "modules.h" // VME modules

// version
const std::string APP_VERSION = "3.0";

// rows of trace data written at once.
const int T_NROWS_PER_WRITE = 1;

// type for run meta
typedef struct RunMetaData
{
  uint32_t number;   // run number
  char title[80];    // title
  time_t ts0;        // begin timestamp
  time_t ts1;        // end timestamp
  char date0[80];    // datetime from ts0
  char date1[80];    // datetime from ts1
  uint32_t dt;       // elapsed time in sec
  char fmt[6];       // ring format
  uint64_t n_events; // total events
  uint64_t n_frags;  // total fragments
  char ctrl_type[6]; // device controller type
} RunMetaData;

// runmetadata compound data type
const std::string META_DATA_RUN_NUMBER("Run_Number");           // uint32_t
const std::string META_DATA_TITLE("Title");                     // char[80]
const std::string META_DATA_TIMESTAMP_0("Begin_Timestamp");     // time_t
const std::string META_DATA_TIMESTAMP_1("End_Timestamp");       // time_t
const std::string META_DATA_ELAPSEDTIME("Elapsed_Seconds");     // uint32_t
const std::string META_DATA_DATETIME_0("Begin_Datetime");       // char[80]
const std::string META_DATA_DATETIME_1("End_Datetime");         // char[80]
const std::string META_DATA_TOTAL_EVENTS("Total_Events");       // uint64_t
const std::string META_DATA_TOTAL_FRAGMENTS("Total_Fragments"); // uint64_t
const std::string META_DATA_RING_FORMAT("Ring_Format");         // char[6]
const std::string META_DATA_CTRL_TYPE("Controller");            // char[6]

// type for physics fragment data (DDAS)
typedef struct FragmentData
{
  uint64_t fragment_id;
  uint64_t event_id;
  uint32_t crate_id;
  uint32_t slot_id;
  uint32_t channel_id;
  double timestamp;
  uint64_t coarse_time;
  uint32_t energy;
  uint32_t finish_code;
  uint32_t overflow_code;
  uint32_t cfd_fail_bit;
  uint32_t trace_length;
  uint32_t adc_frequency;
  uint32_t adc_resolution;
  uint16_t adc_over_underflow;
} FragmentData;

// fragmentdata H5 compound data type (DDAS)
// See Also: http://docs.nscl.msu.edu/daq/newsite/ddas-1.1/DDASHit_8h_source.html
// directly readout parameters
const std::string FRAGMENT_DATA_CRATE_ID("Crate_ID");                      // uint32_t
const std::string FRAGMENT_DATA_SLOT_ID("Slot_ID");                        // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_ID("Channel_ID");                  // uint32_t
const std::string FRAGMENT_DATA_TIMESTAMP("Time");                         // double
const std::string FRAGMENT_DATA_COARSE_TIME("Coarse_Time");                // uint64_t
const std::string FRAGMENT_DATA_ENERGY("Energy");                          // uint32_t
const std::string FRAGMENT_DATA_FINISH_CODE("Err_Pileup");                 // uint32_t
const std::string FRAGMENT_DATA_OVERFLOW_CODE("Overflow_Code");            // uint32_t
const std::string FRAGMENT_DATA_CFD_FAIL_BIT("Err_CFD_Fail");              // uint32_t
const std::string FRAGMENT_DATA_TRACE_LENGTH("Trace_Length");              // uint32_t
const std::string FRAGMENT_DATA_MODMSPS("ADC_Frequency");                  // uint32_t
const std::string FRAGMENT_DATA_ADC_RESOLUTION("ADC_Resolution");          // uint32_t
const std::string FRAGMENT_DATA_ADC_OVER_UNDER_FLOW("Err_ADC_Saturation"); // uint16_t
// others, including post-processed parameters
const std::string FRAGMENT_DATA_FRAGMENT_ID("Fragment_ID"); // uint64_t
const std::string FRAGMENT_DATA_EVENT_ID("Event_ID");       // uint64_t

// group name for physics event data
const std::string PHYSICS_EVENT_GROUP_NAME("Events");

// dataset names under physics event data group (ddas)
const std::string FRAGMENTS_DSET_NAME("FragmentData");
const std::string TRACES_DSET_NAME("TraceData");

// group name for scaler data
const std::string SCALER_GROUP_NAME("Scalers");
const std::string SCALERS_DSET_NAME("ScalerData");

/*
H5 data structure:
 / (root group)
 /attr: metadata
 /Events # physics events
  - FragmentData (dataset)
  - TraceData (dataset)
/Scalers # scalers
  - ScalerData (dataset)
*/

// fragment data rank
const int FRAGMENT_DATA_RANK = 1;

// trace data rank
const int TRACE_DATA_RANK = 2;

// scaler data rank
const int SCALER_DATA_RANK = 1;

// type for DDAS scaler data
typedef struct DDASScalerData
{
  time_t timestamp;
  uint32_t slot_id; // beginning at 2
  char datetime[80];
  uint32_t crate_id;
  uint32_t inc_raw_ch00; // raw trigger incremental count at slot-id.ch0
  uint32_t inc_val_ch00; // validated trigger incremental count
  uint32_t inc_raw_ch01;
  uint32_t inc_val_ch01;
  uint32_t inc_raw_ch02;
  uint32_t inc_val_ch02;
  uint32_t inc_raw_ch03;
  uint32_t inc_val_ch03;
  uint32_t inc_raw_ch04;
  uint32_t inc_val_ch04;
  uint32_t inc_raw_ch05;
  uint32_t inc_val_ch05;
  uint32_t inc_raw_ch06;
  uint32_t inc_val_ch06;
  uint32_t inc_raw_ch07;
  uint32_t inc_val_ch07;
  uint32_t inc_raw_ch08;
  uint32_t inc_val_ch08;
  uint32_t inc_raw_ch09;
  uint32_t inc_val_ch09;
  uint32_t inc_raw_ch10;
  uint32_t inc_val_ch10;
  uint32_t inc_raw_ch11;
  uint32_t inc_val_ch11;
  uint32_t inc_raw_ch12;
  uint32_t inc_val_ch12;
  uint32_t inc_raw_ch13;
  uint32_t inc_val_ch13;
  uint32_t inc_raw_ch14;
  uint32_t inc_val_ch14;
  uint32_t inc_raw_ch15;
  uint32_t inc_val_ch15;
} DDASScalerData;

// type for non-DDAS (e.g. VME) scaler data
typedef struct ScalerData
{
  time_t timestamp;
  char datetime[80];
  uint32_t slot_id;      // beginning at 2 (?assume the same as ddas)
  uint32_t inc_raw_ch00; // raw trigger incremental count at ch0
  uint32_t inc_val_ch00; // validated trigger incremental count
  uint32_t inc_raw_ch01;
  uint32_t inc_val_ch01;
  uint32_t inc_raw_ch02;
  uint32_t inc_val_ch02;
  uint32_t inc_raw_ch03;
  uint32_t inc_val_ch03;
  uint32_t inc_raw_ch04;
  uint32_t inc_val_ch04;
  uint32_t inc_raw_ch05;
  uint32_t inc_val_ch05;
  uint32_t inc_raw_ch06;
  uint32_t inc_val_ch06;
  uint32_t inc_raw_ch07;
  uint32_t inc_val_ch07;
  uint32_t inc_raw_ch08;
  uint32_t inc_val_ch08;
  uint32_t inc_raw_ch09;
  uint32_t inc_val_ch09;
  uint32_t inc_raw_ch10;
  uint32_t inc_val_ch10;
  uint32_t inc_raw_ch11;
  uint32_t inc_val_ch11;
  uint32_t inc_raw_ch12;
  uint32_t inc_val_ch12;
  uint32_t inc_raw_ch13;
  uint32_t inc_val_ch13;
  uint32_t inc_raw_ch14;
  uint32_t inc_val_ch14;
  uint32_t inc_raw_ch15;
  uint32_t inc_val_ch15;
} ScalerData;

// for H5 compound data type
const std::string SCALER_TIMESTAMP("Timestamp"); // time_t
const std::string SCALER_DATETIME("Datetime");   // char[80]
const std::string SCALER_CRATE_ID("Crate_ID");   //uint32_t
const std::string SCALER_SLOT_ID("Slot_ID");     //uint32_t
const std::string SCALER_RAW_CH00("Raw_Ch00");   //uint32_t
const std::string SCALER_VAL_CH00("Val_Ch00");   //uint32_t
const std::string SCALER_RAW_CH01("Raw_Ch01");   //uint32_t
const std::string SCALER_VAL_CH01("Val_Ch01");   //uint32_t
const std::string SCALER_RAW_CH02("Raw_Ch02");   //uint32_t
const std::string SCALER_VAL_CH02("Val_Ch02");   //uint32_t
const std::string SCALER_RAW_CH03("Raw_Ch03");   //uint32_t
const std::string SCALER_VAL_CH03("Val_Ch03");   //uint32_t
const std::string SCALER_RAW_CH04("Raw_Ch04");   //uint32_t
const std::string SCALER_VAL_CH04("Val_Ch04");   //uint32_t
const std::string SCALER_RAW_CH05("Raw_Ch05");   //uint32_t
const std::string SCALER_VAL_CH05("Val_Ch05");   //uint32_t
const std::string SCALER_RAW_CH06("Raw_Ch06");   //uint32_t
const std::string SCALER_VAL_CH06("Val_Ch06");   //uint32_t
const std::string SCALER_RAW_CH07("Raw_Ch07");   //uint32_t
const std::string SCALER_VAL_CH07("Val_Ch07");   //uint32_t
const std::string SCALER_RAW_CH08("Raw_Ch08");   //uint32_t
const std::string SCALER_VAL_CH08("Val_Ch08");   //uint32_t
const std::string SCALER_RAW_CH09("Raw_Ch09");   //uint32_t
const std::string SCALER_VAL_CH09("Val_Ch09");   //uint32_t
const std::string SCALER_RAW_CH10("Raw_Ch10");   //uint32_t
const std::string SCALER_VAL_CH10("Val_Ch10");   //uint32_t
const std::string SCALER_RAW_CH11("Raw_Ch11");   //uint32_t
const std::string SCALER_VAL_CH11("Val_Ch11");   //uint32_t
const std::string SCALER_RAW_CH12("Raw_Ch12");   //uint32_t
const std::string SCALER_VAL_CH12("Val_Ch12");   //uint32_t
const std::string SCALER_RAW_CH13("Raw_Ch13");   //uint32_t
const std::string SCALER_VAL_CH13("Val_Ch13");   //uint32_t
const std::string SCALER_RAW_CH14("Raw_Ch14");   //uint32_t
const std::string SCALER_VAL_CH14("Val_Ch14");   //uint32_t
const std::string SCALER_RAW_CH15("Raw_Ch15");   //uint32_t
const std::string SCALER_VAL_CH15("Val_Ch15");   //uint32_t

/**
 * Argument parser
 *
 *
 */
class ArgumentParser
{

public:
  ArgumentParser();
  ~ArgumentParser();

  // default param value
  std::string get_default(std::string name);

  // parse args
  void parse(int argc, char **argv);

  // has option?
  bool has_option(std::string opt);

  // print help, if -h or missing mandatory: -i, return True if printed
  bool print_help_if_possible();

  // print version info, if -v, return True if printed
  bool print_version_if_possible();

  void print_help();
  void print_examples();

  // version, if --version
  std::string get_version();

  // gzip compress level
  int get_gzip_compress_level();

  // compress method
  std::string get_compress_method();

  // chunk size, pointer
  hsize_t *get_chunk_dims();

  // max evt number
  uint64_t get_max_evt();

  // verbosity level, 0 (default), 1, 2
  // 0: no messages
  // 1: scaler
  // 2: scaler, events
  // 3: all
  int get_verbosity();

  // input data source (original arg)
  std::string get_input_data_source();

  // valida uri input data source
  std::string get_input_data_source_as_uri();

  // validate input data source
  bool validate_input_data_source();

  // full path for output data file
  std::string get_output_filepath();

  // print all args
  void print_all_args();

  // prog name and desc
  std::string get_progname();
  std::string get_progdesc();

  // controller type, DDAS or VME
  std::string get_ctrl_type();

  // module list for VME
  std::vector<std::string> *get_module_list();

  // exclude item type(s)
  std::vector<uint16_t> *get_exclude_types();

private:
  std::string m_prog;
  std::string m_progname;
  std::string m_progdesc;
  argh::parser m_parser;

  std::map<std::string, std::string> m_default_params;

  int m_gzip_comp_level;
  std::string m_comp_method;
  hsize_t *m_chunk_dims;
  uint64_t m_max_evt;

  std::vector<uint16_t> *m_exclude_types = new std::vector<uint16_t>();
  std::vector<std::string> *m_module_list = new std::vector<std::string>();

  int m_verbosity;

  std::string m_ctrl_type;
  std::string m_ifname;
  std::string m_ofname;
  std::string m_ifname_full;
  std::string m_ofname_full;

  void init_options();
};

/**
 * write_metadata: Write metadata.
 * 
 */
bool write_metadata(RunMetaData &metadata, H5::H5File *group);

/**
 * write_fragdata: Write fragment data
 *  
 */
bool write_fragdata(std::vector<FragmentData> *pfragdata, H5::Group *group);

/**
 * write_fragdata_vme: Write fragment data
 *  
 */
bool write_fragdata_vme(std::vector<uint64_t> *pfragdata, H5::Group *group,
                        uint64_t &frag_cnt, std::string &comp_meth, int &gfactor);

/**
 * write_tracedata: Write trace data
 * 
 */
bool write_tracedata(std::vector<uint16_t> *ptracedata, H5::Group *group, uint64_t &frag_cnt,
                     hsize_t *chunk_dims, std::string &comp_meth, int &gfactor);

/**
 * write_scalerdata: Write scaler data
 *  
 */
bool write_scalerdata(std::vector<uint32_t> *pscalerdata, H5::Group *group,
                      std::vector<uint32_t> *pscalerlen, std::vector<time_t> *pscalerts,
                      std::string &stype);

bool _write_ddas_scalerdata(std::vector<uint32_t> *pscalerdata, H5::Group *group,
                            std::vector<uint32_t> *pscalerlen, std::vector<time_t> *pscalerts);

bool _write_vme_scalerdata(std::vector<uint32_t> *pscalerdata, H5::Group *group,
                           std::vector<uint32_t> *pscalerlen, std::vector<time_t> *pscalerts);

// VME V785 fragmentdata
const std::string V785_FRAGMENT_DATA_FRAG_ID("Fragment_ID");     // uint64_t
const std::string V785_FRAGMENT_DATA_EVENT_ID("Event_ID");       // uint64_t
const std::string V785_FRAGMENT_DATA_GEO("GEO");                 // uint16_t
const std::string V785_FRAGMENT_DATA_CRATE_ID("Crate_ID");       // uint16_t
const std::string V785_FRAGMENT_DATA_TOTAL_CHANNELS("Channels"); // uint16_t
const std::string V785_FRAGMENT_DATA_CHANNEL_ID("Channel_ID");   // uint16_t
const std::string V785_FRAGMENT_DATA_VALUE("Value");             // uint16_t
const std::string V785_FRAGMENT_DATA_OVERFLOW("Overflow");       // uint16_t
const std::string V785_FRAGMENT_DATA_THRESHOLD("Threshold");     // uint16_t