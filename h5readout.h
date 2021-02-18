#pragma once

#include <cstdint>
#include <vector>

#include "H5Cpp.h"

// type for physics fragment data
typedef struct FragmentData {
    uint64_t event_id;
    uint32_t crate_id;
    uint32_t slot_id;
    uint32_t channel_id;
    double timestamp;
    uint64_t coarse_time;
    uint32_t energy;
    uint32_t time_high;
    uint32_t time_low;
    uint32_t time_cfd;
    uint32_t finish_code;
    uint32_t channel_length;
    uint32_t channel_header_length;
    uint32_t overflow_code;
    uint32_t cfd_trig_source_bit;
    uint32_t cfd_fail_bit;
    uint32_t trace_length;
    uint32_t adc_frequency;
    uint32_t hardware_revision;
    uint32_t adc_resolution;
    uint16_t adc_over_underflow;
} FragmentData;

// type for physics event data
typedef struct EventData {
   uint16_t event_id;
   uint16_t frag_cnt;
   std::vector<FragmentData> *fragdata = new std::vector<FragmentData>();
} EventData;

// for H5 compound data type
// See Also: http://docs.nscl.msu.edu/daq/newsite/ddas-1.1/DDASHit_8h_source.html
const std::string FRAGMENT_DATA_CRATE_ID             ( "Crate ID"             ); // uint32_t
const std::string FRAGMENT_DATA_SLOT_ID              ( "Slot ID"              ); // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_ID           ( "Channel ID"           ); // uint32_t
const std::string FRAGMENT_DATA_TIMESTAMP            ( "Timestamp"            ); // double
const std::string FRAGMENT_DATA_COARSE_TIME          ( "Coarse Time"          ); // uint64_t
const std::string FRAGMENT_DATA_ENERGY               ( "Energy"               ); // uint32_t
const std::string FRAGMENT_DATA_TIME_HIGH            ( "Time High"            ); // uint32_t
const std::string FRAGMENT_DATA_TIME_LOW             ( "Time Low"             ); // uint32_t
const std::string FRAGMENT_DATA_TIME_CFD             ( "Time CFD"             ); // uint32_t
const std::string FRAGMENT_DATA_FINISH_CODE          ( "Finish Code"          ); // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_LENGTH       ( "Channel Length"       ); // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_HEADER_LENGTH( "Channel Header Length"); // uint32_t
const std::string FRAGMENT_DATA_OVERFLOW_CODE        ( "Overflow Code"        ); // uint32_t
const std::string FRAGMENT_DATA_CFD_TRIG_SOURCE_BIT  ( "CFD Trig Source Bit"  ); // uint32_t
const std::string FRAGMENT_DATA_CFD_FAIL_BIT         ( "CFD Fail Bit"         ); // uint32_t
const std::string FRAGMENT_DATA_TRACE_LENGTH         ( "Trace Length"         ); // uint32_t
const std::string FRAGMENT_DATA_MODMSPS              ( "ADC Frequency"        ); // uint32_t
const std::string FRAGMENT_DATA_HW_REVISION          ( "Hardware Revision"    ); // uint32_t
const std::string FRAGMENT_DATA_ADC_RESOLUTION       ( "ADC Resolution"       ); // uint32_t
const std::string FRAGMENT_DATA_ADC_OVER_UNDER_FLOW  ( "ADC Over/Underflow"   ); // uint16_t
// others
const std::string FRAGMENT_DATA_EVENT_ID             ( "Event ID"             ); // uint64_t
// dataset name
const std::string FRAGMENTS_DSET_NAME( "Fragments" );
// fragment data rank
const int FRAGMENT_DATA_RANK = 1;

/*
data structure:
 / (root group)
   Fragments (dataset)
   Traces (dataset)
*/
