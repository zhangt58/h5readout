#include <stdexcept>
#include "modules.h"

namespace VME
{

    V785Hit::V785Hit() : m_geo(0),
                         m_channel_id(0),
                         m_channels(0),
                         m_crate_id(0),
                         m_value(0),
                         m_overflow_code(0),
                         m_threshold_code(0),
                         m_event_id(0),
                         m_type(0)
    {
    }

    V785Hit::V785Hit(uint32_t &record)
    {
        reset();
        // GEO
        m_geo = (record >> 27) & 0x1F;
        // type, b010: header, b000: datum, b100: EOB
        m_type = (record >> 24) & 0x07;

        // only for header
        if (is_header())
        {
            // crate_id
            m_crate_id = (record >> 16) & 0xFF;
            // channels
            m_channels = (record >> 8) & 0x3F;
        }
        else if (is_datum()) // datum
        {
            // channel_id (V785), for V785N (record >> 17) & 0xF;
            m_channel_id = (record >> 16) & 0x1F;
            // threshold code
            m_threshold_code = (record >> 13) & 1;
            // overflow code
            m_overflow_code = (record >> 12) & 1;
            // value
            m_value = record & 0xFFF;
        }
        else if (is_eob()) // EOB
        {
            // event id
            m_event_id = record & 0xFFFFFF;
        }
    }

    bool V785Hit::is_header()
    {
        return m_type == 2;
    }
    bool V785Hit::is_datum()
    {
        return m_type == 0;
    }
    bool V785Hit::is_eob()
    {
        return m_type == 4;
    }

    void V785Hit::reset()
    {
        m_geo = 0;
        m_crate_id = 0;
        m_type = 0;
        m_channel_id = 0;
        m_value = 0;
        m_overflow_code = 0;
        m_threshold_code = 0;
        m_event_id = 0;
        m_channels = 0;
    }

    V785Hit::~V785Hit() {}

    uint16_t V785Hit::get_geo() { return m_geo; }
    uint16_t V785Hit::get_type() { return m_type; }
    uint16_t V785Hit::get_channel_id() { return m_channel_id; }
    uint16_t V785Hit::get_crate_id() { return m_crate_id; }
    uint16_t V785Hit::get_value() { return m_value; }
    uint16_t V785Hit::get_channels() { return m_channels; }
    uint16_t V785Hit::get_overflow_code() { return m_overflow_code; }
    uint16_t V785Hit::get_threshold_code() { return m_threshold_code; }
    uint64_t V785Hit::get_event_id() { return m_event_id; }
}