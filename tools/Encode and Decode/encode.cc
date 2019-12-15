#include <stdint.h>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <complex>
#include <liquid/liquid.h>
#include <unordered_map>


typedef std::complex<float>			gr_complex;

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
            "\e[1m\e[36m"
};

inline int32_t wrap_index(int32_t i, int32_t n) {
    return ((i % n) + n) % n;
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


template <typename T>
inline void print_vector(std::ostream &out, const T *v, const std::string &prefix, const int size, const int element_len_bits)
{
    out << prefix << ": ";

    for (int i = 0; i < size; i++)
        out << to_bin(v[i], element_len_bits) << ", ";

    out << std::endl
        << std::flush;
}

inline uint32_t gray_decode(uint32_t x)
{
    for (uint32_t bit = 1u << 31; bit >= 2; bit >>= 1)
    {
        if (x & bit)
        {
            x ^= bit >> 1;
        }
    }

    return x;
}

inline uint32_t rotl(uint32_t bits, uint32_t count = 1u, const uint32_t size = 8u)
{
    const uint32_t len_mask = (1u << size) - 1u;

    count %= size;    // Limit bit rotate count to size
    bits &= len_mask; // Limit given bits to size

    return ((bits << count) & len_mask) | (bits >> (size - count));
}

static inline gr_complex
gr_expj(float phase)
{
  float	t_imag, t_real;
  sincosf(phase, &t_imag, &t_real);
  return gr_complex(t_real, t_imag);
}

typedef struct __attribute__((__packed__)) loratap_rssi
{
    uint8_t packet_rssi;  /* LoRa packet
RSSI, if snr >= 0 then dBm value is -139 + packet_rssi, otherwise dBm
value is -139 + packet_rssi * .25 */
    uint8_t max_rssi;     /* LoRa
receiver max RSSI (dBm value is -139 + rssi) */
    uint8_t current_rssi; /* LoRa
receiver current RSSI (dBm value is -139 + rssi) */
    uint8_t snr;          /* LoRa SNR (dB
value is (snr[two's complement])/4) */
} loratap_rssi_t;

typedef struct __attribute__((__packed__)) loratap_channel
{
    uint32_t frequency; /* LoRa
frequency (Hz) */
    uint8_t bandwidth;  /* Channel
bandwidth (KHz) in 125 KHz steps */
    uint8_t sf;         /* LoRa SF
(sf_t) [7, 8, 9, 10, 11, 12] */
} loratap_channel_t;

typedef struct __attribute__((__packed__)) loratap_header
{
    uint8_t lt_version; /* LoRatap
header version */
    uint8_t lt_padding;
    uint16_t lt_length; /* LoRatap
header length */
    loratap_channel_t channel;
    loratap_rssi_t rssi;
    uint8_t sync_word; /* LoRa radio
sync word [0x34 = LoRaWAN] */
} loratap_header_t;

typedef struct __attribute__((__packed__)) loraphy_header
{
    uint8_t length;
    uint8_t crc_msn : 4;
    uint8_t has_mac_crc : 1;
    uint8_t cr : 3;
    uint8_t crc_lsn : 4;
    uint8_t reserved : 4;
} loraphy_header_t;

typedef struct __attribute__((__packed__)) loraconf
{
    loratap_header_t tap;
    loraphy_header_t phy;
} loraconf_t;

std::unordered_map<uint8_t, std::vector<gr_complex>> d_chirps;
fec d_h48_fec;
uint32_t d_samples_per_second;
uint32_t d_bw;
bool d_explicit;
bool d_reduced_rate;
double d_dt;
uint8_t d_osr;
uint16_t d_num_preamble_symbols;
double d_chirp_phi0;
std::vector<gr_complex> d_sample_buffer;



bool parse_packet_conf(loraconf_t &conf, uint8_t *packet, uint32_t packet_len)
{
    if (packet_len <= sizeof(loraconf_t))
    {
        return false;
    }

    memcpy(&conf.tap, packet, sizeof(loratap_header_t));
    memcpy(&conf.phy, packet + sizeof(loratap_header_t), sizeof(loraphy_header_t));

    return true;
}
/**
 *  Whitening sequence
 */
const uint8_t prng_payload[] = {
    0xff, 0xff, 0x2d, 0xff, 0x78, 0xff, 0xe1, 0xff, 0x00, 0xff, 0xd2, 0x2d, 0x55, 0x78, 0x4b, 0xe1, 0x66, 0x00, 0x1e, 0xd2, 0xff, 0x55, 0x2d, 0x4b, 0x78, 0x66, 0xe1, 0x1e, 0xd2, 0xff, 0x87, 0x2d, 0xcc, 0x78, 0xaa, 0xe1, 0xb4, 0xd2, 0x99, 0x87, 0xe1, 0xcc, 0x00, 0xaa, 0x00, 0xb4, 0x00, 0x99, 0x00, 0xe1, 0xd2, 0x00, 0x55, 0x00, 0x99, 0x00, 0xe1, 0x00, 0xd2, 0xd2, 0x87, 0x55, 0x1e, 0x99, 0x2d, 0xe1, 0x78, 0xd2, 0xe1, 0x87, 0xd2, 0x1e, 0x55, 0x2d, 0x99, 0x78, 0x33, 0xe1, 0x55, 0xd2, 0x4b, 0x55, 0x66, 0x99, 0x1e, 0x33, 0x2d, 0x55, 0x78, 0x4b, 0xe1, 0x66, 0x00, 0x1e, 0x00, 0x2d, 0x00, 0x78, 0xd2, 0xe1, 0x87, 0x00, 0xcc, 0x00, 0x78, 0x00, 0x33, 0xd2, 0x55, 0x87, 0x99, 0xcc, 0x33, 0x78, 0x55, 0x33, 0x99, 0x55, 0x33, 0x99, 0x87, 0x33, 0xcc, 0x55, 0xaa, 0x99, 0x66, 0x33, 0x1e, 0x87, 0x2d, 0xcc, 0x78, 0xaa, 0x33, 0x66, 0x55, 0x1e, 0x99, 0x2d, 0xe1, 0x78, 0x00, 0x33, 0x00, 0x55, 0xd2, 0x99, 0x55, 0xe1, 0x4b, 0x00, 0xb4, 0x00, 0x4b, 0xd2, 0x66, 0x55, 0xcc, 0x4b, 0xaa, 0xb4, 0x66, 0x4b, 0xcc, 0x66, 0xaa, 0xcc, 0xb4, 0xaa, 0x4b, 0x66, 0x66, 0xcc, 0xcc, 0xaa, 0x78, 0xb4, 0x33, 0x4b, 0x55, 0x66, 0x4b, 0xcc, 0x66, 0x78, 0xcc, 0x33, 0x78, 0x55, 0xe1, 0x4b, 0x00, 0x66, 0xd2, 0xcc, 0x87, 0x78, 0x1e, 0xe1, 0xff, 0x00, 0xff, 0xd2, 0x2d, 0x87, 0xaa, 0x1e, 0x66, 0xff, 0xcc, 0xff, 0xaa, 0x2d, 0x66, 0xaa, 0x1e, 0x66, 0xff, 0xcc, 0x2d, 0xaa, 0xaa, 0x66, 0xb4, 0x1e, 0x4b, 0xff, 0x66, 0x2d, 0x1e, 0xaa, 0x2d, 0xb4, 0xaa, 0x4b, 0xb4, 0x66, 0x99, 0x1e, 0xe1, 0x2d, 0xd2, 0xaa, 0x55, 0xb4, 0x99, 0x99, 0xe1, 0xe1, 0x00, 0xd2, 0xd2, 0x55, 0x87, 0x99, 0xcc, 0xe1, 0xaa, 0x00, 0x66, 0xd2, 0xcc, 0x87, 0x78, 0xcc, 0xe1, 0xaa, 0xd2, 0x66, 0x87, 0xcc, 0x1e, 0x78, 0xff, 0xe1, 0x2d, 0xd2, 0x78, 0x87, 0x33, 0x1e, 0x87, 0xff, 0x1e, 0x2d, 0x2d, 0x78, 0x78, 0x33, 0x33, 0x87, 0x87, 0x1e, 0xcc, 0x2d, 0x78, 0x78, 0xe1, 0x33, 0xd2, 0x87, 0x55, 0xcc, 0x4b, 0x78, 0x66, 0xe1, 0xcc, 0xd2, 0xaa, 0x55, 0xb4, 0x4b, 0x99, 0x66, 0x33, 0xcc, 0x55, 0xaa, 0x99, 0xb4, 0xe1, 0x99, 0xd2, 0x33, 0x55, 0x55, 0x4b, 0x99, 0xb4, 0xe1, 0x99, 0xd2, 0x33, 0x55, 0x55, 0x4b, 0x4b, 0xb4, 0xb4, 0x99, 0x4b, 0x33, 0xb4, 0x55, 0x99, 0x4b, 0x33, 0xb4, 0x87, 0x4b, 0x1e, 0xb4, 0x2d, 0x99, 0xaa, 0x33, 0x66, 0xc7, 0x1e, 0x1e, 0x2d, 0x2d, 0xaa, 0xaa, 0x66, 0x66, 0xcc, 0x1e, 0x78, 0x2d, 0x33, 0xaa, 0x87, 0x66, 0x1e, 0xcc, 0xff, 0x78, 0x2d, 0x33, 0xaa, 0x87, 0x66, 0x1e, 0x1e, 0xff, 0xff, 0x2d, 0xff, 0xaa, 0xff, 0x66, 0x2d, 0x1e, 0xaa, 0xff, 0xb4, 0xff, 0x99, 0xff, 0x33, 0x2d, 0x87, 0xaa, 0xcc, 0xb4, 0x78, 0x99, 0x33, 0x33, 0x87, 0x87, 0xcc, 0xcc, 0xaa, 0x78, 0xb4, 0x33, 0x4b, 0x87, 0xb4, 0xcc, 0x99, 0xaa, 0xe1, 0xb4, 0xd2, 0x4b, 0x87, 0xb4, 0xcc, 0x99, 0x78, 0xe1, 0xe1, 0xd2, 0x00, 0x87, 0x00, 0xcc, 0xd2, 0x78, 0x87, 0xe1, 0x1e, 0x00, 0x2d, 0x00, 0xaa, 0xd2, 0xb4, 0x87, 0x4b, 0x1e, 0xb4, 0x2d, 0x4b, 0xaa, 0xb4, 0xb4, 0x4b, 0x4b, 0x66, 0xb4, 0x1e, 0x4b, 0xff, 0xb4, 0xff, 0x4b, 0x2d, 0x66, 0x78, 0x1e, 0x33, 0xff, 0x55, 0xff, 0x4b, 0x2d, 0xb4, 0x78, 0x99, 0x33, 0xe1, 0x55, 0x00, 0x4b, 0xd2, 0xb4, 0x55, 0x99, 0x99, 0xe1, 0x33, 0x00, 0x87, 0xd2, 0x1e, 0x55, 0xff, 0x99, 0xff, 0x33, 0xff, 0x87, 0xff, 0x1e, 0x00, 0x00, 0x00, 0x00, 0x87, 0xe1, 0xaa, 0xcc,
};

void transmit_chirp(bool up, uint8_t sf, uint16_t symbol, bool quarter = false)
{
    const uint32_t chips_per_symbol = pow(2, sf);
    const double symbols_per_second = (double)d_bw / chips_per_symbol;
    const uint32_t samples_per_symbol = d_samples_per_second / symbols_per_second;
    const double T = 0.5 * d_bw * symbols_per_second;
    const double f0 = -(d_bw / 2.0);
    double pre_dir = 2.0 * M_PI;
    double t;

    if (!up)
        pre_dir *= -1;

    std::vector<gr_complex> chirp(samples_per_symbol);
    double phase = 0;
    for (uint32_t i = 0u; i < samples_per_symbol; i++)
    {
        t = d_dt * ((i + (d_osr * symbol)) % samples_per_symbol);
        phase = d_chirp_phi0 + (pre_dir * t * (f0 + T * t));
        chirp[i] = gr_expj(phase);
    }

    // Add chirp to buffer
    if (quarter)
        d_sample_buffer.insert(d_sample_buffer.end(), chirp.begin(), chirp.begin() + chirp.size() / 4);
    else
        d_sample_buffer.insert(d_sample_buffer.end(), chirp.begin(), chirp.end());

    // Set phase
    d_chirp_phi0 = phase;
}

void shuffle(uint8_t *data, uint32_t data_len, const uint8_t *shuffle_pattern)
{
    for (uint32_t i = 0u; i < data_len; i++)
    {
        uint8_t result = 0u;

        for (uint32_t j = 0u; j < 8; j++)
        {
            result |= !!(data[i] & (1u << shuffle_pattern[j])) << j;
        }

        data[i] = result;
    }
}

uint8_t interleave_block(uint16_t *symbols, uint8_t *data, uint8_t sf, uint8_t cr, bool reduced_rate)
{
    if (reduced_rate)
        sf -= 2;

    // Determine symbols for this block
    for (uint8_t symbol = 0; symbol < 4 + cr; symbol++)
    {
        for (int8_t bit = sf - 1; bit >= 0; bit--)
        {
            symbols[symbol] |= ((data[bit] >> symbol) & 0x01) << bit;
        }

        int32_t to_rotate = sf - symbol;
        if (to_rotate < 0)
            to_rotate += sf;

        symbols[symbol] = rotl(symbols[symbol], to_rotate, sf);
    }

    // Rotate to interleave
    std::vector<uint16_t> symbols_v(symbols, symbols + (4 + cr));
    print_interleave_matrix(std::cout, symbols_v, sf);
    print_vector(std::cout, symbols, "Chips", 4 + cr, sf);

    // Determine bins
    std::cout << "Bins: ";
    for (uint8_t symbol = 0; symbol < 4 + cr; symbol++)
    {
        symbols[symbol] = gray_decode(symbols[symbol]);
        if (reduced_rate)
            symbols[symbol] <<= 2;
        std::cout << (int)symbols[symbol] << ", ";
    }
    std::cout << std::endl;

    return sf;
}

void nibble_swap(uint8_t* encoded, uint32_t length) {
        for(uint32_t i = 0; i+1 < length; i += 2) {
            uint8_t tmp = encoded[i];
            encoded[i] = encoded[i+1];
            encoded[i+1] = tmp;
        }
}

inline void whiten(uint8_t* input, const uint8_t* sequence, uint32_t length) {
            for(uint32_t i = 0; i < length; i++) {
                input[i] = input[i] ^ sequence[i];
            }
        }

void transmit_packet(loraconf_t& conf, uint8_t* packet, bool up = true) { // TODO: clean up
        uint32_t packet_length = conf.phy.length + sizeof(loraphy_header_t); //
        uint32_t num_bytes = packet_length*2;

        // uint32_t num_symbols = 8.0 * sizeof(loraphy_header_t) / ((4.0/8)) + 


        uint8_t encoded[num_bytes];
        uint32_t num_symbols = num_bytes * ((4.0+conf.phy.cr) / conf.tap.channel.sf) + 0.5;
        uint32_t encoded_offset = 0;
        uint32_t packet_offset = 0;

        std::cout <<"phy length: " <<(int) conf.phy.length << std::endl;
        std::cout <<"phy header len: " <<(int) sizeof(loraphy_header_t) << std::endl;
    
        std::cout <<"numbytes: " << (int) num_bytes << std::endl;
        std::cout <<"calc: " <<(float) ((4.0+conf.phy.cr) / conf.tap.channel.sf) << std::endl;
        std::cout <<"num_symbols: " <<(int) num_symbols << std::endl;
        std::cout <<"cr: " << (int) conf.phy.cr << std::endl;
        std::cout <<"sf: " << (int)conf.tap.channel.sf << std::endl;

        // Add preamble symbols to queue
        for(uint32_t i = 0; i < d_num_preamble_symbols; i++) {
            transmit_chirp(up, conf.tap.channel.sf, 0);
        }

        // Add sync words to queue
        uint8_t sync_word = 0x12;
        uint32_t sync_offset_1 = ((sync_word & 0xf0) >> 4) * pow(2, conf.tap.channel.sf) * d_osr / 32;
        uint32_t sync_offset_2 = (sync_word & 0x0f) * pow(2, conf.tap.channel.sf) * d_osr / 32;
        transmit_chirp(up, conf.tap.channel.sf, sync_offset_1);
        transmit_chirp(up, conf.tap.channel.sf, sync_offset_2);

        // Add SFD to queue
        transmit_chirp(!up, conf.tap.channel.sf, 0);
        transmit_chirp(!up, conf.tap.channel.sf, 0);
        transmit_chirp(!up, conf.tap.channel.sf, 0, true);

        // If explicit header, add one block (= SF codewords) to queue in reduced rate mode (and always 4/8)
        if(d_explicit) {
            fec_encode(d_h48_fec, 3, packet, encoded); // Header is always 4/8
            packet_offset = 3;
            encoded_offset = 5;
        }

        // Add remaining blocks to queue
        print_vector(std::cout, packet, "Packet", packet_length, 8);

        fec_encode(d_h48_fec, packet_length - packet_offset, packet+packet_offset, encoded+encoded_offset); // TODO: change to appropriate scheme
        nibble_swap(encoded+encoded_offset, num_bytes-encoded_offset);
        print_vector(std::cout, encoded, "Encoded", num_bytes, 8);

        whiten(encoded+encoded_offset, prng_payload, num_bytes-encoded_offset);
        print_vector(std::cout, encoded, "Whitened", num_bytes, 8);

        const uint8_t shuffle_pattern[] = {1, 2, 3, 5, 4, 0, 6, 7};
        shuffle(encoded, num_bytes, shuffle_pattern);
        print_vector(std::cout, encoded, "Shuffled", num_bytes, 8);

        // Interleaving
        uint16_t symbols[num_symbols];
        memset(symbols, 0x00, num_symbols * sizeof(uint16_t));
        uint32_t symbols_done = 0;
        uint32_t interleave_offset = 0;

        if(d_explicit) {
            interleave_offset += interleave_block(symbols, encoded, conf.tap.channel.sf, 4, true);
            symbols_done += 8;
        }

        while(interleave_offset+conf.tap.channel.sf <= num_bytes) { // TODO: needs to be exact number of bytes
            interleave_offset += interleave_block(symbols+symbols_done, encoded+interleave_offset, conf.tap.channel.sf, conf.phy.cr, d_reduced_rate);
            symbols_done += 4 + conf.phy.cr;
        }

        // Transmission
        for(uint32_t i = 0; i < num_symbols; i++) {
            transmit_chirp(up, conf.tap.channel.sf, symbols[i]);
        }
    }

int main()
{
    d_osr = 8;
    d_samples_per_second = 125000 * d_osr;
    d_num_preamble_symbols = 8;
    d_bw = 125000;
    d_explicit = true;
    d_reduced_rate = false;
    d_chirp_phi0 = -M_PI;

    d_dt = 1.0f / d_samples_per_second;
    d_sample_buffer.reserve(d_samples_per_second); // Allocate for one second of samples

    fec_scheme fs = LIQUID_FEC_HAMMING84;
    d_h48_fec = fec_create(fs, NULL);


    // Temporary         ve  pa      le              fq  bw  sf  pr  mr  cr  sn  sy  H1  H1  H1
    char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x17\x91\xa0\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x20\x21\x22\xb8\x73";
    
    //char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x03\x30\x30\x01\x02\x03\xda\xe2";
    // Header + Goodbye! + CRC,  (08 30 00 47 6f 6f 64 62 79 65 21 3f e4) i.e what the lora_receiver has decoded from the arduino raw uplink
    //                            le H2 H3 G  o  o  d  b   y  e  ! (CRC ) le is correctly 8 i.e Goodbye! has 8 characters
     //char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x09\x00\x00\x00\x00\x12\x03\x10\x00\x30\x03\x10\x64\x62\x79\x65\x21\x3f\xe4";

   // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x03\x30\x30\x41\x43\x4b\x03\x56\xeb"; // ACK

    //char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x03\x30\x30\x01\x02\x03\xda\xe2"; // 123


    
    loraconf_t conf;
    memset(&conf, 0, sizeof(loraconf_t));

    if(!parse_packet_conf(conf, (uint8_t*)test_pkt, sizeof(test_pkt))) {
    std::cerr << "Malformed LoRa packet received" << std::endl;
    exit(1);
    }

    transmit_packet(conf, (uint8_t*)(test_pkt + sizeof(loratap_header_t)), true);


    std::string pathname("output.bin");

std::ofstream textout(pathname.c_str(), std::ios::out | std::ios::binary);
textout.write((const char*)&d_sample_buffer[0], d_sample_buffer.size()*sizeof(gr_complex)); // numer of items in vector times the size of each item in vector





textout.close();
    return 0;
}
