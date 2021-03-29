#pragma once

#include <cstdint>
#include <vector>
#include "H5Cpp.h"

// type for physics fragment data (V785)
typedef struct FragmentData_V785
{
    uint64_t event_id;
    uint16_t geo;
    uint16_t crate_id;
    uint16_t n_channels;
    std::vector<uint16_t> val_channel; // values
    hvl_t val_channel_handle;
    std::vector<uint16_t> ov_channel; // overflow_code;
    hvl_t ov_channel_handle;
    std::vector<uint16_t> un_channel; // threshold_code;
    hvl_t un_channel_handle;
} FragmentData_V785;

namespace VME
{
    class V785Hit
    {
    private:
        uint16_t m_geo;
        uint16_t m_channel_id;
        uint16_t m_crate_id;
        uint16_t m_value;
        uint16_t m_channels;
        uint16_t m_overflow_code;
        uint16_t m_threshold_code;
        uint64_t m_event_id;
        uint16_t m_type;

    public:
        V785Hit();
        V785Hit(uint32_t &);
        ~V785Hit();

        void reset();
        uint16_t get_geo();
        uint16_t get_channel_id();
        uint16_t get_crate_id();
        uint16_t get_value();
        uint16_t get_channels();
        uint16_t get_overflow_code();
        uint16_t get_threshold_code();
        uint64_t get_event_id();
        uint16_t get_type();

        bool is_header();
        bool is_datum();
        bool is_eob();
    };
};
