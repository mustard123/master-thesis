#include <stdint.h>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <unordered_map>
#include <liquid/liquid.h>
#include <vector>
#include <fstream>
#include <volk/volk.h>
#include <boost/circular_buffer.hpp>
#include <gnuradio/sincos.h>
#include <numeric>
#include <pmt/pmt.h>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/parsers.hpp>
#include <experimental/filesystem>
#include <bitset>

#define MAC_CRC_SIZE 2u
#define MAX_PWR_QUEUE_SIZE 4

typedef std::complex<float> gr_complex;

typedef enum sf
{
    SF7 = 7,
    SF8,
    SF9,
    SF10,
    SF11,
    SF12
} sf_t;

typedef struct __attribute__((__packed__)) loratap_channel
{
    uint32_t frequency; /* LoRa frequency (Hz) */
    uint8_t bandwidth;  /* Channel bandwidth (KHz) in 125 KHz steps */
    uint8_t sf;         /* LoRa SF (sf_t) [7, 8, 9, 10, 11, 12] */
} loratap_channel_t;

typedef struct __attribute__((__packed__)) loratap_rssi
{
    uint8_t packet_rssi;  /* LoRa packet RSSI, if snr >= 0 then dBm value is -139 + packet_rssi, otherwise dBm value is -139 + packet_rssi * .25 */
    uint8_t max_rssi;     /* LoRa receiver max RSSI (dBm value is -139 + rssi) */
    uint8_t current_rssi; /* LoRa receiver current RSSI (dBm value is -139 + rssi) */
    uint8_t snr;          /* LoRa SNR (dB value is (snr[two's complement])/4) */
} loratap_rssi_t;

typedef struct __attribute__((__packed__)) loratap_header
{
    uint8_t lt_version; /* LoRatap header version */
    uint8_t lt_padding;
    uint16_t lt_length; /* LoRatap header length */
    loratap_channel_t channel;
    loratap_rssi_t rssi;
    uint8_t sync_word; /* LoRa radio sync word [0x34 = LoRaWAN] */
} loratap_header_t;

const uint8_t prng_header[] = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
/**
         *  Whitening sequence
         */
const uint8_t prng_payload[] = {
    0xff,
    0xff,
    0x2d,
    0xff,
    0x78,
    0xff,
    0xe1,
    0xff,
    0x00,
    0xff,
    0xd2,
    0x2d,
    0x55,
    0x78,
    0x4b,
    0xe1,
    0x66,
    0x00,
    0x1e,
    0xd2,
    0xff,
    0x55,
    0x2d,
    0x4b,
    0x78,
    0x66,
    0xe1,
    0x1e,
    0xd2,
    0xff,
    0x87,
    0x2d,
    0xcc,
    0x78,
    0xaa,
    0xe1,
    0xb4,
    0xd2,
    0x99,
    0x87,
    0xe1,
    0xcc,
    0x00,
    0xaa,
    0x00,
    0xb4,
    0x00,
    0x99,
    0x00,
    0xe1,
    0xd2,
    0x00,
    0x55,
    0x00,
    0x99,
    0x00,
    0xe1,
    0x00,
    0xd2,
    0xd2,
    0x87,
    0x55,
    0x1e,
    0x99,
    0x2d,
    0xe1,
    0x78,
    0xd2,
    0xe1,
    0x87,
    0xd2,
    0x1e,
    0x55,
    0x2d,
    0x99,
    0x78,
    0x33,
    0xe1,
    0x55,
    0xd2,
    0x4b,
    0x55,
    0x66,
    0x99,
    0x1e,
    0x33,
    0x2d,
    0x55,
    0x78,
    0x4b,
    0xe1,
    0x66,
    0x00,
    0x1e,
    0x00,
    0x2d,
    0x00,
    0x78,
    0xd2,
    0xe1,
    0x87,
    0x00,
    0xcc,
    0x00,
    0x78,
    0x00,
    0x33,
    0xd2,
    0x55,
    0x87,
    0x99,
    0xcc,
    0x33,
    0x78,
    0x55,
    0x33,
    0x99,
    0x55,
    0x33,
    0x99,
    0x87,
    0x33,
    0xcc,
    0x55,
    0xaa,
    0x99,
    0x66,
    0x33,
    0x1e,
    0x87,
    0x2d,
    0xcc,
    0x78,
    0xaa,
    0x33,
    0x66,
    0x55,
    0x1e,
    0x99,
    0x2d,
    0xe1,
    0x78,
    0x00,
    0x33,
    0x00,
    0x55,
    0xd2,
    0x99,
    0x55,
    0xe1,
    0x4b,
    0x00,
    0xb4,
    0x00,
    0x4b,
    0xd2,
    0x66,
    0x55,
    0xcc,
    0x4b,
    0xaa,
    0xb4,
    0x66,
    0x4b,
    0xcc,
    0x66,
    0xaa,
    0xcc,
    0xb4,
    0xaa,
    0x4b,
    0x66,
    0x66,
    0xcc,
    0xcc,
    0xaa,
    0x78,
    0xb4,
    0x33,
    0x4b,
    0x55,
    0x66,
    0x4b,
    0xcc,
    0x66,
    0x78,
    0xcc,
    0x33,
    0x78,
    0x55,
    0xe1,
    0x4b,
    0x00,
    0x66,
    0xd2,
    0xcc,
    0x87,
    0x78,
    0x1e,
    0xe1,
    0xff,
    0x00,
    0xff,
    0xd2,
    0x2d,
    0x87,
    0xaa,
    0x1e,
    0x66,
    0xff,
    0xcc,
    0xff,
    0xaa,
    0x2d,
    0x66,
    0xaa,
    0x1e,
    0x66,
    0xff,
    0xcc,
    0x2d,
    0xaa,
    0xaa,
    0x66,
    0xb4,
    0x1e,
    0x4b,
    0xff,
    0x66,
    0x2d,
    0x1e,
    0xaa,
    0x2d,
    0xb4,
    0xaa,
    0x4b,
    0xb4,
    0x66,
    0x99,
    0x1e,
    0xe1,
    0x2d,
    0xd2,
    0xaa,
    0x55,
    0xb4,
    0x99,
    0x99,
    0xe1,
    0xe1,
    0x00,
    0xd2,
    0xd2,
    0x55,
    0x87,
    0x99,
    0xcc,
    0xe1,
    0xaa,
    0x00,
    0x66,
    0xd2,
    0xcc,
    0x87,
    0x78,
    0xcc,
    0xe1,
    0xaa,
    0xd2,
    0x66,
    0x87,
    0xcc,
    0x1e,
    0x78,
    0xff,
    0xe1,
    0x2d,
    0xd2,
    0x78,
    0x87,
    0x33,
    0x1e,
    0x87,
    0xff,
    0x1e,
    0x2d,
    0x2d,
    0x78,
    0x78,
    0x33,
    0x33,
    0x87,
    0x87,
    0x1e,
    0xcc,
    0x2d,
    0x78,
    0x78,
    0xe1,
    0x33,
    0xd2,
    0x87,
    0x55,
    0xcc,
    0x4b,
    0x78,
    0x66,
    0xe1,
    0xcc,
    0xd2,
    0xaa,
    0x55,
    0xb4,
    0x4b,
    0x99,
    0x66,
    0x33,
    0xcc,
    0x55,
    0xaa,
    0x99,
    0xb4,
    0xe1,
    0x99,
    0xd2,
    0x33,
    0x55,
    0x55,
    0x4b,
    0x99,
    0xb4,
    0xe1,
    0x99,
    0xd2,
    0x33,
    0x55,
    0x55,
    0x4b,
    0x4b,
    0xb4,
    0xb4,
    0x99,
    0x4b,
    0x33,
    0xb4,
    0x55,
    0x99,
    0x4b,
    0x33,
    0xb4,
    0x87,
    0x4b,
    0x1e,
    0xb4,
    0x2d,
    0x99,
    0xaa,
    0x33,
    0x66,
    0xc7,
    0x1e,
    0x1e,
    0x2d,
    0x2d,
    0xaa,
    0xaa,
    0x66,
    0x66,
    0xcc,
    0x1e,
    0x78,
    0x2d,
    0x33,
    0xaa,
    0x87,
    0x66,
    0x1e,
    0xcc,
    0xff,
    0x78,
    0x2d,
    0x33,
    0xaa,
    0x87,
    0x66,
    0x1e,
    0x1e,
    0xff,
    0xff,
    0x2d,
    0xff,
    0xaa,
    0xff,
    0x66,
    0x2d,
    0x1e,
    0xaa,
    0xff,
    0xb4,
    0xff,
    0x99,
    0xff,
    0x33,
    0x2d,
    0x87,
    0xaa,
    0xcc,
    0xb4,
    0x78,
    0x99,
    0x33,
    0x33,
    0x87,
    0x87,
    0xcc,
    0xcc,
    0xaa,
    0x78,
    0xb4,
    0x33,
    0x4b,
    0x87,
    0xb4,
    0xcc,
    0x99,
    0xaa,
    0xe1,
    0xb4,
    0xd2,
    0x4b,
    0x87,
    0xb4,
    0xcc,
    0x99,
    0x78,
    0xe1,
    0xe1,
    0xd2,
    0x00,
    0x87,
    0x00,
    0xcc,
    0xd2,
    0x78,
    0x87,
    0xe1,
    0x1e,
    0x00,
    0x2d,
    0x00,
    0xaa,
    0xd2,
    0xb4,
    0x87,
    0x4b,
    0x1e,
    0xb4,
    0x2d,
    0x4b,
    0xaa,
    0xb4,
    0xb4,
    0x4b,
    0x4b,
    0x66,
    0xb4,
    0x1e,
    0x4b,
    0xff,
    0xb4,
    0xff,
    0x4b,
    0x2d,
    0x66,
    0x78,
    0x1e,
    0x33,
    0xff,
    0x55,
    0xff,
    0x4b,
    0x2d,
    0xb4,
    0x78,
    0x99,
    0x33,
    0xe1,
    0x55,
    0x00,
    0x4b,
    0xd2,
    0xb4,
    0x55,
    0x99,
    0x99,
    0xe1,
    0x33,
    0x00,
    0x87,
    0xd2,
    0x1e,
    0x55,
    0xff,
    0x99,
    0xff,
    0x33,
    0xff,
    0x87,
    0xff,
    0x1e,
    0x00,
    0x00,
    0x00,
    0x00,
    0x87,
    0xe1,
    0xaa,
    0xcc,
};

static std::vector<std::string> term_colors = {
    "\e[0m",
    "\e[1m\e[91m",
    "\e[1m\e[92m",
    "\e[1m\e[93m",
    "\e[1m\e[94m",
    "\e[1m\e[95m",
    "\e[1m\e[96m",
    "\e[1m\e[97m",
    "\e[1m\e[31m",
    "\e[1m\e[32m",
    "\e[1m\e[33m",
    "\e[1m\e[34m",
    "\e[1m\e[35m",
    "\e[1m\e[36m"};

inline int32_t wrap_index(int32_t i, int32_t n)
{
    return ((i % n) + n) % n;
}

template <typename T>
inline void print_vector_hex(std::ostream& out, const T* v, const uint32_t size, bool endline) {
    for (uint32_t i = 0u; i < size; i++) {
        out << " " << std::hex << std::setw(2) << std::setfill('0') << (int)v[i];
    }

    if(endline)
        out << std::endl;

    out << std::flush;
}

template <typename T>
inline std::string to_bin(const T v, const uint32_t bitwidth)
{
#ifdef LSB_FIRST
    const uint64_t maxpow = bitwidth ? (1ull << (bitwidth - 1)) : 0;
    uint64_t mask;

    std::string result = "";

    for (mask = 0x1; mask <= maxpow; mask <<= 1)
    {
        result += (v & mask) ? "1" : "0";
    }
#else
    uint64_t mask = bitwidth ? (1ull << bitwidth) : 0;

    std::string result = "";

    while (mask >>= 1)
    {
        result += (v & mask) ? "1" : "0";
    }
#endif

    return result;
}

template <typename T>
inline void print_vector_bin(std::ostream &out, const std::vector<T> &v, const std::string &prefix, const int element_len_bits)
{
    out << prefix << ": ";

    for (const T &x : v)
        out << to_bin(x, element_len_bits) << ", ";

    out << std::endl
        << std::flush;
}
template <typename T>
inline void print_interleave_matrix(std::ostream &out, const std::vector<T> &v, const uint32_t sf)
{
    uint32_t cr = v.size();

    for (uint32_t i = 0; i < cr; i++)
        out << "-";
    out << std::endl;

    out << "LSB" << std::endl;

    for (int32_t i = sf - 1; i >= 0; i--)
    {
        for (int32_t j = 0; j < (int32_t)cr; j++)
        {
            out << term_colors[wrap_index(j - i, (int32_t)sf) + 1] << to_bin(v[j], sf)[i] << term_colors[0];
        }
        out << std::endl;
    }

    out << "MSB" << std::endl;

    for (uint32_t i = 0; i < cr; i++)
        out << "-";
    out << std::endl;

    out << std::flush;
}
typedef struct __attribute__((__packed__)) loraphy_header
{
    uint8_t length;
    uint8_t crc_msn : 4;
    uint8_t has_mac_crc : 1;
    uint8_t cr : 3;
    uint8_t crc_lsn : 4;
    uint8_t reserved : 4;
} loraphy_header_t;

/**
         *  \brief  **DecoderState** : Each state the LoRa decoder can be in.
         */
enum class DecoderState
{
    DETECT,
    SYNC,
    FIND_SFD,
    PAUSE,
    DECODE_HEADER,
    DECODE_PAYLOAD,
    STOP
};

DecoderState d_state; ///< Holds the current state of the decoder (state machine).

std::vector<gr_complex> d_downchirp;  ///< The complex ideal downchirp.
std::vector<float> d_downchirp_ifreq; ///< The instantaneous frequency of the ideal downchirp.

std::vector<gr_complex> d_upchirp;    ///< The complex ideal upchirp.
std::vector<float> d_upchirp_ifreq;   ///< The instantaneous frequency of the ideal upchirp.
std::vector<float> d_upchirp_ifreq_v; ///< The instantaneous frequency of the ideal upchirp.

std::vector<gr_complex> d_fft;     ///< Vector containing the FFT resuls.
std::vector<gr_complex> d_mult_hf; ///< Vector containing the FFT decimation.
std::vector<gr_complex> d_tmp;     ///< Vector containing the FFT decimation.

bool d_implicit;                                               ///< Implicit header mode.
bool d_reduced_rate;                                           ///< Use reduced rate (only configurable in implicit header mode).
uint8_t d_sf;                                                  ///< The Spreading Factor.
uint32_t d_bw;                                                 ///< The receiver bandwidth (fixed to `125kHz`).
loraphy_header_t d_phdr;                                       ///< LoRa PHY header.
uint16_t d_mac_crc;                                            ///< The MAC CRC.
double d_bits_per_second;                                      ///< Indicator of how many bits are transferred each second.
uint32_t d_delay_after_sync;                                   ///< The number of samples to skip in `DecoderState::PAUSE`.
uint32_t d_samples_per_second;                                 ///< The number of samples taken per second by GNU Radio.
double d_symbols_per_second;                                   ///< Indicator of how many symbols (read: chirps) are transferred each second.
double d_bits_per_symbol;                                      ///< The number of bits each of the symbols contain.
uint32_t d_samples_per_symbol;                                 ///< The number of samples in one symbol.
double d_period;                                               ///< Period of the symbol.
uint32_t d_number_of_bins;                                     ///< Indicates in how many parts or bins a symbol is decimated, i.e. the max value to decode out of one payload symbol.
uint32_t d_number_of_bins_hdr;                                 ///< Indicates in how many parts or bins a HDR symbol is decimated, i.e. the max value to decode out of one HDR symbol.
int32_t d_payload_symbols;                                     ///< The number of symbols needed to decode the payload. Calculated from an indicator in the HDR.
uint32_t d_payload_length;                                     ///< The number of words after decoding the HDR or payload. Calculated from an indicator in the HDR.
uint32_t d_corr_fails;                                         ///< Indicates how many times the correlation failed. After some tries, the state will revert to `DecoderState::DETECT`.
float d_energy_threshold;                                      ///< The absolute threshold to distinguish signal from noise.
const uint8_t *d_whitening_sequence;                           ///< A pointer to the whitening sequence to be used in decoding. Determined by the SF in the ctor.
float d_snr;                                                   ///< Signal to noise ratio
boost::circular_buffer<float> d_pwr_queue(MAX_PWR_QUEUE_SIZE); ///< Queue holding symbol power values

int num_of_detects = 0;

std::vector<uint32_t> d_words;           ///< Vector containing the demodulated words.
std::vector<uint8_t> d_demodulated;      ///< Vector containing the words after deinterleaving.
std::vector<uint8_t> d_words_deshuffled; ///< Vector containing the words after deshuffling.
std::vector<uint8_t> d_words_dewhitened; ///< Vector containing the words after dewhitening.
std::vector<uint8_t> d_decoded;          ///< Vector containing the words after Hamming decode or the final decoded words.

std::ofstream d_debug_samples; ///< Debug utputstream for complex values.
std::ofstream d_debug;         ///< Outputstream for the debug log.

fftplan d_q;   ///< The LiquidDSP::FFT_Plan.
fftplan d_qr;  ///< The LiquidDSP::FFT_Plan in reverse.
fec d_h48_fec; ///< LiquidDSP Hamming 4/8 FEC.

uint32_t d_decim_factor; ///< The number of samples (data points) in each bin.
float d_cfo_estimation;  ///< An estimation for the current Center Frequency Offset.
double d_dt;             ///< Indicates how fast the frequency changes in a symbol (chirp).
bool d_enable_fine_sync; ///< Enable drift correction
int32_t d_fine_sync;     ///< Amount of drift correction to apply for next symbol

gr_complex *input;
uintmax_t file_size;
uintmax_t read_samples = 0;
std::ofstream outfile;

static inline gr_complex
gr_expj(float phase)
{
    float t_imag, t_real;
    gr::sincosf(phase, &t_imag, &t_real);
    return gr_complex(t_real, t_imag);
}

inline uint32_t build_packet(uint8_t *buffer, uint32_t offset, const void *header, uint32_t header_size)
{
    memcpy(buffer + offset, header, header_size);
    offset += header_size;

    return offset;
}

inline uint32_t rotl(uint32_t bits, uint32_t count = 1u, const uint32_t size = 8u)
{
    const uint32_t len_mask = (1u << size) - 1u;

    count %= size;    // Limit bit rotate count to size
    bits &= len_mask; // Limit given bits to size

    return ((bits << count) & len_mask) | (bits >> (size - count));
}

inline void swap_nibbles(uint8_t *array, uint32_t length)
{
    for (uint32_t i = 0; i < length; i++)
    {
        array[i] = ((array[i] & 0x0f) << 4) | ((array[i] & 0xf0) >> 4);
    }
}

inline void instantaneous_frequency(const gr_complex *in_samples, float *out_ifreq, const uint32_t window)
{
    if (window < 2u)
    {
        std::cerr << "[LoRa Decoder] WARNING : window size < 2 !" << std::endl;
        return;
    }

    /* instantaneous_phase */
    for (uint32_t i = 1u; i < window; i++)
    {
        const float iphase_1 = std::arg(in_samples[i - 1]);
        float iphase_2 = std::arg(in_samples[i]);

        // Unwrapped loops from liquid_unwrap_phase
        while ((iphase_2 - iphase_1) > M_PI)
            iphase_2 -= 2.0f * M_PI;
        while ((iphase_2 - iphase_1) < -M_PI)
            iphase_2 += 2.0f * M_PI;

        out_ifreq[i - 1] = iphase_2 - iphase_1;
    }

    // Make sure there is no strong gradient if this value is accessed by mistake
    out_ifreq[window - 1] = out_ifreq[window - 2];
}

inline void instantaneous_phase(const gr_complex *in_samples, float *out_iphase, const uint32_t window)
{
    out_iphase[0] = std::arg(in_samples[0]);

    for (uint32_t i = 1u; i < window; i++)
    {
        out_iphase[i] = std::arg(in_samples[i]);
        // = the same as atan2(imag(in_samples[i]),real(in_samples[i]));

        // Unwrapped loops from liquid_unwrap_phase
        while ((out_iphase[i] - out_iphase[i - 1]) > M_PI)
            out_iphase[i] -= 2.0f * M_PI;
        while ((out_iphase[i] - out_iphase[i - 1]) < -M_PI)
            out_iphase[i] += 2.0f * M_PI;
    }
}

inline uint32_t select_bits(const uint32_t data, const uint8_t *indices, const uint8_t n)
{
    uint32_t r = 0u;

    for (uint8_t i = 0u; i < n; ++i)
        r |= (data & (1u << indices[i])) ? (1u << i) : 0u;

    return r;
}

float detect_preamble_autocorr(const gr_complex *samples, const uint32_t window)
{
    const gr_complex *chirp1 = samples;
    const gr_complex *chirp2 = samples + d_samples_per_symbol;
    float magsq_chirp1[window];
    float magsq_chirp2[window];
    float energy_chirp1 = 0;
    float energy_chirp2 = 0;
    float autocorr = 0;
    gr_complex dot_product;

    volk_32fc_x2_conjugate_dot_prod_32fc(&dot_product, chirp1, chirp2, window);
    volk_32fc_magnitude_squared_32f(magsq_chirp1, chirp1, window);
    volk_32fc_magnitude_squared_32f(magsq_chirp2, chirp2, window);
    volk_32f_accumulator_s32f(&energy_chirp1, magsq_chirp1, window);
    volk_32f_accumulator_s32f(&energy_chirp2, magsq_chirp2, window);

    // When using implicit mode, stop when energy is halved.
    d_energy_threshold = energy_chirp2 / 2.0f;

    // For calculating the SNR later on
    d_pwr_queue.push_back(energy_chirp1 / d_samples_per_symbol);

    // Autocorr value
    autocorr = abs(dot_product / gr_complex(sqrt(energy_chirp1 * energy_chirp2), 0));

    return autocorr;
}

float experimental_determine_cfo(const gr_complex *samples, uint32_t window)
{
    gr_complex mult[window];
    float mult_ifreq[window];

    volk_32fc_x2_multiply_32fc(mult, samples, &d_downchirp[0], window);
    instantaneous_frequency(mult, mult_ifreq, window);

    return mult_ifreq[256] / (2.0 * M_PI) * d_samples_per_second;
}

void extract_data_only(bool is_header)
{
    static const uint8_t data_indices[4] = {1, 2, 3, 5};
    uint32_t len = d_words_dewhitened.size();

    for (uint32_t i = 0u; i < len; i += 2u)
    {
        const uint8_t d2 = (i + 1u < len) ? select_bits(d_words_dewhitened[i + 1u], data_indices, 4u) & 0xFF : 0u;
        const uint8_t d1 = (select_bits(d_words_dewhitened[i], data_indices, 4u) & 0xFF);

        if (is_header)
            d_decoded.push_back((d1 << 4u) | d2);
        else
            d_decoded.push_back((d2 << 4u) | d1);
    }
}

float cross_correlate_ifreq_fast(const float *samples_ifreq, const float *ideal_chirp, const uint32_t window)
{
    float result = 0;
    volk_32f_x2_dot_prod_32f(&result, samples_ifreq, ideal_chirp, window);
    return result;
}

void fine_sync(const gr_complex *in_samples, uint32_t bin_idx, int32_t search_space)
{
    int32_t shift_ref = (bin_idx + 1) * d_decim_factor;
    float samples_ifreq[d_samples_per_symbol];
    float max_correlation = 0.0f;
    int32_t lag = 0;

    instantaneous_frequency(in_samples, samples_ifreq, d_samples_per_symbol);

    for (int32_t i = -search_space + 1; i < search_space; i++)
    {
        //float c = cross_correlate_fast(in_samples, &d_upchirp_v[shift_ref+i+d_samples_per_symbol], d_samples_per_symbol);
        float c = cross_correlate_ifreq_fast(samples_ifreq, &d_upchirp_ifreq_v[shift_ref + i + d_samples_per_symbol], d_samples_per_symbol);
        if (c > max_correlation)
        {
            max_correlation = c;
            lag = i;
        }
    }

#ifdef DEBUG
    //d_debug << "FINE: " << -lag << std::endl;
#endif

    d_fine_sync = -lag;

    //if(abs(d_fine_sync) >= d_decim_factor / 2)
    //    d_fine_sync = 0;
    //d_fine_sync = 0;
}

typedef std::complex<float> gr_complex;

void deshuffle(const uint8_t *shuffle_pattern, const bool is_header)
{
    const uint32_t to_decode = is_header ? 5u : d_demodulated.size();
    const uint32_t len = sizeof(shuffle_pattern) / sizeof(uint8_t);
    uint8_t result;

    for (uint32_t i = 0u; i < to_decode; i++)
    {
        result = 0u;

        for (uint32_t j = 0u; j < len; j++)
        {
            result |= !!(d_demodulated[i] & (1u << shuffle_pattern[j])) << j;
        }

        d_words_deshuffled.push_back(result);
    }

    // print_vector_bin(d_debug, d_words_deshuffled, "S", sizeof(uint8_t) * 8);

    // We're done with these words
    if (is_header)
    {
        d_demodulated.erase(d_demodulated.begin(), d_demodulated.begin() + 5u);
        d_words_deshuffled.push_back(0);
    }
    else
    {
        d_demodulated.clear();
    }
}

void dewhiten(const uint8_t *prng)
{
    const uint32_t len = d_words_deshuffled.size();

    for (uint32_t i = 0u; i < len; i++)
    {
        uint8_t xor_b = d_words_deshuffled[i] ^ prng[i];
        d_words_dewhitened.push_back(xor_b);
    }

    print_vector_bin(d_debug, d_words_dewhitened, "W", sizeof(uint8_t) * 8);

    d_words_deshuffled.clear();
}

void hamming_decode(bool is_header)
{
    switch (d_phdr.cr)
    {
    case 4:
    case 3:
    { // Hamming(8,4) or Hamming(7,4)
        //hamming_decode_soft(is_header);
        uint32_t n = ceil(d_words_dewhitened.size() * 4.0f / (4.0f + d_phdr.cr));
        uint8_t decoded[n];

        fec_decode(d_h48_fec, n, &d_words_dewhitened[0], decoded);
        if (!is_header)
            swap_nibbles(decoded, n);
        d_decoded.assign(decoded, decoded + n);
        break;
    }
    case 2:
    case 1:
    { // Hamming(6,4) or Hamming(5,4)
        // TODO: Report parity error to the user
        extract_data_only(is_header);
        break;
    }
    }

    d_words_dewhitened.clear();
}

void samples_to_file(const std::string path, const gr_complex *v, const uint32_t length, const uint32_t elem_size)
{
#ifdef DEBUG
    std::ofstream out_file;
    out_file.open(path.c_str(), std::ios::out | std::ios::binary);

    //for(std::vector<gr_complex>::const_iterator it = v.begin(); it != v.end(); ++it) {
    for (uint32_t i = 0u; i < length; i++)
    {
        out_file.write(reinterpret_cast<const char *>(&v[i]), elem_size);
    }

    out_file.close();
#else
    (void)path;
    (void)v;
    (void)length;
    (void)elem_size;
#endif
}

void build_ideal_chirps(void)
{
    d_downchirp.resize(d_samples_per_symbol);
    d_upchirp.resize(d_samples_per_symbol);
    d_downchirp_ifreq.resize(d_samples_per_symbol);
    d_upchirp_ifreq.resize(d_samples_per_symbol);
    d_upchirp_ifreq_v.resize(d_samples_per_symbol * 3);
    gr_complex tmp[d_samples_per_symbol * 3];

    const double T = -0.5 * d_bw * d_symbols_per_second;
    const double f0 = (d_bw / 2.0);
    const double pre_dir = 2.0 * M_PI;
    double t;
    gr_complex cmx = gr_complex(1.0f, 1.0f);

    for (uint32_t i = 0u; i < d_samples_per_symbol; i++)
    {
        // Width in number of samples = samples_per_symbol
        // See https://en.wikipedia.org/wiki/Chirp#Linear
        t = d_dt * i;
        d_downchirp[i] = cmx * gr_expj(pre_dir * t * (f0 + T * t));
        d_upchirp[i] = cmx * gr_expj(pre_dir * t * (f0 + T * t) * -1.0f);
    }

    // Store instantaneous frequency
    instantaneous_frequency(&d_downchirp[0], &d_downchirp_ifreq[0], d_samples_per_symbol);
    instantaneous_frequency(&d_upchirp[0], &d_upchirp_ifreq[0], d_samples_per_symbol);

    samples_to_file("/tmp/downchirp", &d_downchirp[0], d_downchirp.size(), sizeof(gr_complex));
    samples_to_file("/tmp/upchirp", &d_upchirp[0], d_upchirp.size(), sizeof(gr_complex));

    // Upchirp sequence
    memcpy(tmp, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
    memcpy(tmp + d_samples_per_symbol, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
    memcpy(tmp + d_samples_per_symbol * 2, &d_upchirp[0], sizeof(gr_complex) * d_samples_per_symbol);
    instantaneous_frequency(tmp, &d_upchirp_ifreq_v[0], d_samples_per_symbol * 3);
}

float sliding_norm_cross_correlate_upchirp(const float *samples_ifreq, const uint32_t window, int32_t *index)
{
    float max_correlation = 0;

    // Cross correlate
    for (uint32_t i = 0; i < window; i++)
    {
        const float max_corr = cross_correlate_ifreq_fast(samples_ifreq + i, &d_upchirp_ifreq[0], window - 1u);

        if (max_corr > max_correlation)
        {
            *index = i;
            max_correlation = max_corr;
        }
    }

    return max_correlation;
}

float detect_upchirp(const gr_complex *samples, const uint32_t window, int32_t *index)
{
    float samples_ifreq[window * 2];
    instantaneous_frequency(samples, samples_ifreq, window * 2);

    return sliding_norm_cross_correlate_upchirp(samples_ifreq, window, index);
}

void deinterleave(const uint32_t ppm)
{
    const uint32_t bits_per_word = d_words.size();
    const uint32_t offset_start = ppm - 1u;

    std::vector<uint8_t> words_deinterleaved(ppm, 0u);

    if (bits_per_word > 8u)
    {
        // Not sure if this can ever occur. It would imply coding rate high than 4/8 e.g. 4/9.
        std::cerr << "[LoRa Decoder] WARNING : Deinterleaver: More than 8 bits per word. uint8_t will not be sufficient!\nBytes need to be stored in intermediate array and then packed into words_deinterleaved!" << std::endl;
        exit(1);
    }

    for (uint32_t i = 0u; i < bits_per_word; i++)
    {
        const uint32_t word = rotl(d_words[i], i, ppm);

        for (uint32_t j = (1u << offset_start), x = offset_start; j; j >>= 1u, x--)
        {
            words_deinterleaved[x] |= !!(word & j) << i;
        }
    }

    print_interleave_matrix(std::cout, d_words, ppm);
    print_vector_bin(std::cout, words_deinterleaved, "D", sizeof(uint8_t) * 8u);

    // Add to demodulated data
    d_demodulated.insert(d_demodulated.end(), words_deinterleaved.begin(), words_deinterleaved.end());

    // Cleanup
    d_words.clear();
}

void determine_snr()
{
    if (d_pwr_queue.size() >= 2)
    {
        float pwr_noise = d_pwr_queue[0];
        float pwr_signal = d_pwr_queue[d_pwr_queue.size() - 1];
        d_snr = pwr_signal / pwr_noise;
    }
}

float stddev(const float *values, const uint32_t len, const float mean)
{
    float variance = 0.0f;

    for (uint32_t i = 0u; i < len; i++)
    {
        const float temp = values[i] - mean;
        variance += temp * temp;
    }

    variance /= (float)len;
    return std::sqrt(variance);
}

float cross_correlate_ifreq(const float *samples_ifreq, const std::vector<float> &ideal_chirp, const uint32_t to_idx)
{
    float result = 0.0f;

    const float average = std::accumulate(samples_ifreq, samples_ifreq + to_idx, 0.0f) / (float)(to_idx);
    const float chirp_avg = std::accumulate(&ideal_chirp[0], &ideal_chirp[to_idx], 0.0f) / (float)(to_idx);
    const float sd = stddev(samples_ifreq, to_idx, average) * stddev(&ideal_chirp[0], to_idx, chirp_avg);

    for (uint32_t i = 0u; i < to_idx; i++)
    {
        result += (samples_ifreq[i] - average) * (ideal_chirp[i] - chirp_avg) / sd;
    }

    result /= (float)(to_idx);

    return result;
}

float detect_downchirp(const gr_complex *samples, const uint32_t window)
{
    float samples_ifreq[window];
    instantaneous_frequency(samples, samples_ifreq, window);

    return cross_correlate_ifreq(samples_ifreq, d_downchirp_ifreq, window - 1u);
}

uint32_t max_frequency_gradient_idx(const gr_complex *samples)
{
    float samples_ifreq[d_samples_per_symbol];
    float samples_ifreq_avg[d_number_of_bins];

    samples_to_file("/tmp/data", &samples[0], d_samples_per_symbol, sizeof(gr_complex));

    instantaneous_frequency(samples, samples_ifreq, d_samples_per_symbol);

    for (uint32_t i = 0; i < d_number_of_bins; i++)
    {
        volk_32f_accumulator_s32f(&samples_ifreq_avg[i], &samples_ifreq[i * d_decim_factor], d_decim_factor);
        samples_ifreq_avg[i] /= d_decim_factor;
    }

    float max_gradient = 0.1f;
    float gradient = 0.0f;
    uint32_t max_index = 0;
    for (uint32_t i = 1u; i < d_number_of_bins; i++)
    {
        gradient = samples_ifreq_avg[i - 1] - samples_ifreq_avg[i];
        if (gradient > max_gradient)
        {
            max_gradient = gradient;
            max_index = i + 1;
        }
    }

    return (d_number_of_bins - max_index) % d_number_of_bins;
}

float determine_energy(const gr_complex *samples)
{
    float magsq_chirp[d_samples_per_symbol];
    float energy_chirp = 0;
    volk_32fc_magnitude_squared_32f(magsq_chirp, samples, d_samples_per_symbol);
    volk_32f_accumulator_s32f(&energy_chirp, magsq_chirp, d_samples_per_symbol);

    return energy_chirp;
}

bool demodulate(const gr_complex *samples, const bool reduced_rate)
{
    // DBGR_TIME_MEASUREMENT_TO_FILE("SFxx_method");

    // DBGR_START_TIME_MEASUREMENT(false, "only");

    uint32_t bin_idx = max_frequency_gradient_idx(samples);
    //uint32_t bin_idx = get_shift_fft(samples);
    if (d_enable_fine_sync)
        fine_sync(samples, bin_idx, std::max(d_decim_factor / 4u, 2u));

    // DBGR_INTERMEDIATE_TIME_MEASUREMENT();

    // Header has additional redundancy
    if (reduced_rate || d_sf > 10)
    {
        bin_idx = std::lround(bin_idx / 4.0f) % d_number_of_bins_hdr;
    }

    // Decode (actually gray encode) the bin to get the symbol value
    const uint32_t word = bin_idx ^ (bin_idx >> 1u);

#ifdef DEBUG
    d_debug << gr::lora::to_bin(word, reduced_rate ? d_sf - 2u : d_sf) << " " << word << " (bin " << bin_idx << ")" << std::endl;
#endif
    d_words.push_back(word);

    std::cout << "Word is:" << word << " Consumed Samples:" << read_samples << std::endl;
    //outfile << std::bitset<32>(word) << ',' << read_samples << std::endl;
    outfile << std::hex << word << ',' << std::dec << read_samples << std::endl;
    // Look for 4+cr symbols and stop
    if (d_words.size() == (4u + d_phdr.cr))
    {
        // Deinterleave
        deinterleave((reduced_rate || d_sf > 10) ? d_sf - 2u : d_sf);

        return true; // Signal that a block is ready for decoding
    }

    return false; // We need more words in order to decode a block
}
void decode(const bool is_header)
{
    static const uint8_t shuffle_pattern[] = {5, 0, 1, 2, 4, 3, 6, 7};

    // For determining shuffle pattern
    //if (!is_header)
    //    values_to_file("/tmp/before_deshuffle", &d_demodulated[0], d_demodulated.size(), 8);

    deshuffle(shuffle_pattern, is_header);

    // For determining whitening sequence
    //if (!is_header)
    //    values_to_file("/tmp/after_deshuffle", &d_words_deshuffled[0], d_words_deshuffled.size(), 8);

    dewhiten(is_header ? prng_header : d_whitening_sequence);

    //if (!is_header)
    //    values_to_file("/tmp/after_dewhiten", &d_words_dewhitened[0], d_words_dewhitened.size(), 8);

    hamming_decode(is_header);
}

void msg_lora_frame(void)
{
    uint32_t len = sizeof(loratap_header_t) + sizeof(loraphy_header_t) + d_payload_length;
    uint32_t offset = 0;
    uint8_t buffer[len];
    loratap_header_t loratap_header;

    memset(buffer, 0, sizeof(uint8_t) * len);
    memset(&loratap_header, 0, sizeof(loratap_header));

    loratap_header.rssi.snr = (uint8_t)(10.0f * log10(d_snr) + 0.5);

    offset = build_packet(buffer, offset, &loratap_header, sizeof(loratap_header_t));
    uint32_t loratap_offset = offset;
    offset = build_packet(buffer, offset, &d_phdr, sizeof(loraphy_header_t));
    offset = build_packet(buffer, offset, &d_decoded[0], d_payload_length);
    if (offset != len)
    {
        std::cerr << "decoder_impl::msg_lora_frame: invalid write" << std::endl;
        exit(1);
    }

    pmt::pmt_t payload_blob = pmt::make_blob(buffer, sizeof(uint8_t) * len);
    std::cout << "The PMT contains (with lora_tap header) " << payload_blob << std::endl;
    std::cout << "Raw LoRa as HEX (without lora_tap header): ";
    print_vector_hex(std::cout, &buffer[loratap_offset], len-loratap_offset, true);
    
}

void conjugate(gr_complex *&input, int number_of_samples)
{
    volk_32fc_conjugate_32fc(input, input, number_of_samples);
}

int work(gr_complex *&input)
{

    d_fine_sync = 0; // Always reset fine sync

    switch (d_state)
    {
    case DecoderState::DETECT:
    {
        num_of_detects++;
        float correlation = detect_preamble_autocorr(input, d_samples_per_symbol);

        if (correlation >= 0.90f)
        {
            determine_snr();
#ifdef DEBUG
            d_debug << "Ca: " << correlation << std::endl;
#endif
            d_corr_fails = 0u;
            d_state = DecoderState::SYNC;
            break;
        }

        // consume_each(d_samples_per_symbol);
        read_samples = read_samples + d_samples_per_symbol;
        input = input + d_samples_per_symbol;
        break;
    }

    case DecoderState::SYNC:
    {
        int i = 0;
        detect_upchirp(input, d_samples_per_symbol, &i);

        //float cfo = experimental_determine_cfo(&input[i], d_samples_per_symbol);
        //pmt::pmt_t kv = pmt::cons(pmt::intern(std::string("cfo")), pmt::from_double(cfo));
        //message_port_pub(pmt::mp("control"), kv);

        samples_to_file("/tmp/detect", &input[i], d_samples_per_symbol, sizeof(gr_complex));

        // consume_each(i);
        read_samples = read_samples + i;
        input = input + i;

        d_state = DecoderState::FIND_SFD;
        break;
    }

    case DecoderState::FIND_SFD:
    {
        const float c = detect_downchirp(input, d_samples_per_symbol);

#ifdef DEBUG
        d_debug << "Cd: " << c << std::endl;
#endif

        if (c > 0.96f)
        {
#ifdef DEBUG
            d_debug << "SYNC: " << c << std::endl;
#endif
            // Debug stuff
            samples_to_file("/tmp/sync", input, d_samples_per_symbol, sizeof(gr_complex));

            d_state = DecoderState::PAUSE;
        }
        else
        {
            if (c < -0.97f)
            {
                fine_sync(input, d_number_of_bins - 1, d_decim_factor * 4);
            }
            else
            {
                d_corr_fails++;
            }

            if (d_corr_fails > 4u)
            {
                d_state = DecoderState::DETECT;
#ifdef DEBUG
                d_debug << "Lost sync" << std::endl;
#endif
            }
        }

        // consume_each((int32_t)d_samples_per_symbol + d_fine_sync);
        read_samples = read_samples + d_samples_per_symbol + d_fine_sync;
        input = input + d_samples_per_symbol + d_fine_sync;

        break;
    }

    case DecoderState::PAUSE:
    {
        if (d_implicit)
        {
            d_state = DecoderState::DECODE_PAYLOAD;
            d_payload_symbols = 1;
        }
        else
        {
            d_state = DecoderState::DECODE_HEADER;
        }
        // consume_each(d_samples_per_symbol + d_delay_after_sync);
        read_samples = read_samples + d_samples_per_symbol + d_delay_after_sync;
        input = input + d_samples_per_symbol + d_delay_after_sync;

        break;
    }

    case DecoderState::DECODE_HEADER:
    {
        d_phdr.cr = 4u;

        if (demodulate(input, true))
        {
            decode(true);
            //print_vector_hex(std::cout, &d_decoded[0], d_decoded.size(), false);
            memcpy(&d_phdr, &d_decoded[0], sizeof(loraphy_header_t));
            if (d_phdr.cr > 4)
                d_phdr.cr = 4;
            d_decoded.clear();

            d_payload_length = d_phdr.length + MAC_CRC_SIZE * d_phdr.has_mac_crc;
            //d_phy_crc = SM(decoded[1], 4, 0xf0) | MS(decoded[2], 0xf0, 4);

            // Calculate number of payload symbols needed
            uint8_t redundancy = (d_sf > 10 ? 2 : 0);
            const int symbols_per_block = d_phdr.cr + 4u;
            const float bits_needed = float(d_payload_length) * 8.0f;
            const float symbols_needed = bits_needed * (symbols_per_block / 4.0f) / float(d_sf - redundancy);
            const int blocks_needed = (int)std::ceil(symbols_needed / symbols_per_block);
            d_payload_symbols = blocks_needed * symbols_per_block;

#ifdef DEBUG
            d_debug << "LEN: " << d_payload_length << " (" << d_payload_symbols << " symbols)" << std::endl;
#endif

            d_state = DecoderState::DECODE_PAYLOAD;
        }

        // consume_each((int32_t)d_samples_per_symbol + d_fine_sync);
        read_samples = read_samples + d_samples_per_symbol + d_fine_sync;
        input = input + d_samples_per_symbol + d_fine_sync;

        break;
    }

    case DecoderState::DECODE_PAYLOAD:
    {
        if (d_implicit && determine_energy(input) < d_energy_threshold)
        {
            d_payload_symbols = 0;
            //d_demodulated.erase(d_demodulated.begin(), d_demodulated.begin() + 7u); // Test for SF 8 with header
            d_payload_length = (int32_t)(d_demodulated.size() / 2);
        }
        else if (demodulate(input, d_implicit && d_reduced_rate))
        {
            if (!d_implicit)
                d_payload_symbols -= (4u + d_phdr.cr);
        }

        if (d_payload_symbols <= 0)
        {
            decode(false);
            //print_vector_hex(std::cout, &d_decoded[0], d_payload_length, true);
            msg_lora_frame();

            d_state = DecoderState::DETECT;
            d_decoded.clear();
            d_words.clear();
            d_words_dewhitened.clear();
            d_words_deshuffled.clear();
            d_demodulated.clear();
        }

        // consume_each((int32_t)d_samples_per_symbol + d_fine_sync);
        read_samples = read_samples + d_samples_per_symbol + d_fine_sync;
        input = input + d_samples_per_symbol + d_fine_sync;

        break;
    }

    case DecoderState::STOP:
    {
        // consume_each(d_samples_per_symbol);
        read_samples = read_samples + d_samples_per_symbol;
        input = input + d_samples_per_symbol;

        break;
    }

    default:
    {
        std::cerr << "[LoRa Decoder] WARNING : No state! Shouldn't happen\n";
        break;
    }
    }

    // DBGR_INTERMEDIATE_TIME_MEASUREMENT();

    // Tell runtime system how many output items we produced.
    return 0;
}

int main(int argc, char **argv)
{

    outfile.open("words.csv", std::ios::out | std::ios::trunc);

    namespace po = boost::program_options;
    namespace pa = boost::filesystem;

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")
    ("file,f", po::value<std::string>(), "path to file")
    ("sf,s", po::value<int>(), "spreading factor, default is 9")
    ("sr,r", po::value<int>(), "sampling rate, how many samples per second, default is 1000000");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("sf"))
    {
        d_sf = (uint8_t)vm["sf"].as<int>();
    }
    else
    {
        d_sf = 9;
    }
    if (vm.count("bw"))
    {
        d_bw = (uint8_t)vm["bw"].as<int>();
    }
    else
    {
        d_bw = 125000;
    }
    if (vm.count("sr"))
    {
        d_samples_per_second = (uint8_t)vm["sr"].as<int>();
    }
    else
    {
        d_samples_per_second = 1000000;
    }
    

#ifdef DEBUG
    d_debug_samples.open("/tmp/grlora_debug", std::ios::out | std::ios::binary);
    d_debug.open("/tmp/grlora_debug_txt", std::ios::out);
    d_dbg.attach();
#endif

    //d_bw = 125000; via program args
    d_implicit = false;
    // d_reduced_rate = reduced_rate;
    d_phdr.cr = 4; // header cr is always 4, phy cr is set in header which gets decoded first to get the value
    d_phdr.has_mac_crc = 1;
    //d_samples_per_second = 1000000; via program args
    d_payload_symbols = 0;
    d_cfo_estimation = 0.0f;
    d_dt = 1.0f / d_samples_per_second;
    // d_sf = 9; via program args
    d_bits_per_second = (double)d_sf * (double)(4.0 / (4.0 + d_phdr.cr)) / (1u << d_sf) * d_bw;
    d_symbols_per_second = (double)d_bw / (1u << d_sf);
    d_period = 1.0f / (double)d_symbols_per_second;
    d_bits_per_symbol = (double)(d_bits_per_second / d_symbols_per_second);
    d_samples_per_symbol = (uint32_t)(d_samples_per_second / d_symbols_per_second);
    d_delay_after_sync = d_samples_per_symbol / 4u;
    d_number_of_bins = (uint32_t)(1u << d_sf);
    d_number_of_bins_hdr = (uint32_t)(1u << (d_sf - 2));
    d_decim_factor = d_samples_per_symbol / d_number_of_bins;
    d_energy_threshold = 0.0f;
    d_whitening_sequence = prng_payload;
    d_fine_sync = 0;
    d_enable_fine_sync = true;
    // set_output_multiple(2 * d_samples_per_symbol);

    // Radio config
    d_state = DecoderState::DETECT;


    if (d_sf < 6 || d_sf > 13)
    {
        std::cout << d_sf << std::endl;
        std::cerr << "[LoRa Decoder] ERROR : Spreading factor should be between 6 and 12 (inclusive)!" << std::endl
                  << "                       Other values are currently not supported." << std::endl;
        exit(1);
    }

    std::cout << "Bits (nominal) per symbol: \t" << d_bits_per_symbol << std::endl;
    std::cout << "Bins per symbol: \t" << d_number_of_bins << std::endl;
    std::cout << "Samples per symbol: \t" << d_samples_per_symbol << std::endl;
    std::cout << "Decimation: \t\t" << d_decim_factor << std::endl;
    std::cout << "Symbols per second: \t\t" << d_symbols_per_second << std::endl;
    std::cout << "Bits per second: \t\t" << d_bits_per_second << std::endl;
    if (!d_enable_fine_sync)
    {
        std::cout << "Warning: clock drift correction disabled" << std::endl;
    }
    if (d_implicit)
    {
        std::cout << "CR: \t\t" << (int)d_phdr.cr << std::endl;
        std::cout << "CRC: \t\t" << (int)d_phdr.has_mac_crc << std::endl;
    }

    // Locally generated chirps
    build_ideal_chirps();

    // FFT decoding preparations
    d_fft.resize(d_samples_per_symbol);
    d_mult_hf.resize(d_samples_per_symbol);
    d_tmp.resize(d_number_of_bins);
    d_q = fft_create_plan(d_samples_per_symbol, &d_mult_hf[0], &d_fft[0], LIQUID_FFT_FORWARD, 0);
    d_qr = fft_create_plan(d_number_of_bins, &d_tmp[0], &d_mult_hf[0], LIQUID_FFT_BACKWARD, 0);

    // Hamming coding
    fec_scheme fs = LIQUID_FEC_HAMMING84;
    d_h48_fec = fec_create(fs, NULL);

    // Register gnuradio ports
    // message_port_register_out(pmt::mp("frames"));
    // message_port_register_out(pmt::mp("control"));

    //const gr_complex *input = (gr_complex *)input_items[0];
    //const gr_complex *raw_input = (gr_complex *) input_items[1]; // Input bypassed by low pass filter

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm.count("file"))
    {
        std::string path = vm["file"].as<std::string>();
        std::string abs_path = boost::filesystem::canonical(path).string();
        std::cout << "File path set to "
                  << abs_path << ".\n";
        try
        {
            file_size = std::experimental::filesystem::file_size(abs_path);
            uintmax_t total_samples = file_size / sizeof(gr_complex);

            std::cout << "File size = " << file_size << '\n';
            std::cout << "Total samples = " << total_samples << '\n';

            int BUFFER_SIZE = file_size; //4* d_samples_per_symbol * sizeof(gr_complex); // allocate space for 4 symbols
            char buffer[BUFFER_SIZE];
            std::ifstream fin(abs_path, std::ios::in | std::ios::binary);
            fin.read(buffer, BUFFER_SIZE);
            input = (gr_complex *)buffer;
            outfile << "word,sample" << std::endl;
            while (read_samples < total_samples - 2 * d_samples_per_symbol) // remove (-2 * d_samples_per_symbol) if decoding signal generate with encoder for getting all the symbols in the csv file. The preamble detection looks two symbols ahead, prevent overflow in detect mode if signal is never detected.
            {
                work(input);
            }

            outfile.close();
        }
        catch (std::experimental::filesystem::filesystem_error &e)
        {
            std::cout << e.what() << '\n';
        }
    }
    else
    {
        std::cout << "No recording specified.\n";
        return 0;
    }
}