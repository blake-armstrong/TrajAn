// Self-contained XTC (GROMACS trajectory) reader/writer.
// XTC uses XDR (big-endian) encoding with optional lossy coordinate compression.
//
// Frame layout (all fields big-endian):
//   Header (56 bytes):
//     int32   magic    = 1995
//     uint32  natoms
//     uint32  step
//     float   time     (ps)
//     float   box[3][3] (nm, row-major)
//     uint32  lsize    (natoms if uncompressed; total compressed-block bytes if > 9)
//
//   If lsize <= 9  →  uncompressed: natoms*3 floats (big-endian, nm)
//   If lsize >  9  →  compressed block (36-byte sub-header + bit-packed data):
//     float   precision
//     int32   minint[3]
//     int32   maxint[3]
//     int32   smallidx
//     int32   nbytes
//     uint8   data[align4(nbytes)]

#include <trajan/io/xtc.h>
#include <trajan/core/util.h>
#include <occ/core/linear_algebra.h>
#include <occ/crystal/unitcell.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace trajan::io {

using occ::crystal::UnitCell;
using trajan::core::Frame;

// ── Endian utilities ──────────────────────────────────────────────────────

static inline bool little_endian() {
  const uint16_t v = 1;
  return *reinterpret_cast<const uint8_t *>(&v) == 1;
}

template <typename T> static T bswap(T v) {
  union {
    T val;
    uint8_t b[sizeof(T)];
  } s, d;
  s.val = v;
  for (size_t i = 0; i < sizeof(T); ++i)
    d.b[i] = s.b[sizeof(T) - 1 - i];
  return d.val;
}

static int32_t rd_i32(std::ifstream &f) {
  int32_t v;
  f.read(reinterpret_cast<char *>(&v), 4);
  return little_endian() ? bswap(v) : v;
}
static uint32_t rd_u32(std::ifstream &f) {
  uint32_t v;
  f.read(reinterpret_cast<char *>(&v), 4);
  return little_endian() ? bswap(v) : v;
}
static float rd_f32(std::ifstream &f) {
  uint32_t v;
  f.read(reinterpret_cast<char *>(&v), 4);
  if (little_endian()) v = bswap(v);
  float r;
  std::memcpy(&r, &v, 4);
  return r;
}
static void wr_i32(std::ofstream &f, int32_t v) {
  if (little_endian()) v = bswap(v);
  f.write(reinterpret_cast<const char *>(&v), 4);
}
static void wr_u32(std::ofstream &f, uint32_t v) {
  if (little_endian()) v = bswap(v);
  f.write(reinterpret_cast<const char *>(&v), 4);
}
static void wr_f32(std::ofstream &f, float v) {
  uint32_t tmp;
  std::memcpy(&tmp, &v, 4);
  if (little_endian()) tmp = bswap(tmp);
  f.write(reinterpret_cast<const char *>(&tmp), 4);
}

static inline int align4(int n) { return (n + 3) & ~3; }

// ── Compression tables ────────────────────────────────────────────────────

static constexpr int32_t XTC_MAGIC = 1995;
static constexpr int FIRSTIDX = 9;
static const int magicints[] = {
    0,       0,       0,       0,       0,       0,       0,        0,
    0,       8,       10,      12,      16,      20,      25,       32,
    40,      50,      64,      80,      101,     128,     161,      203,
    256,     322,     406,     512,     645,     812,     1024,     1290,
    1625,    2048,    2580,    3250,    4096,    5060,    6501,     8192,
    10321,   13003,   16384,   20642,   26007,   32768,   41285,    52015,
    65536,   82570,   104031,  131072,  165140,  208063,  262144,   330280,
    416127,  524287,  660561,  832255,  1048576, 1321122, 1664510,  2097152,
    2642245, 3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085,
    16777216};
static const int LASTIDX = static_cast<int>(sizeof(magicints) / sizeof(int));

static int sizeofint(uint32_t n) {
  uint32_t num = 1;
  int bits = 0;
  while (n >= num && bits < 32) { bits++; num <<= 1; }
  return bits;
}

static int sizeofints(int count, const uint32_t sizes[]) {
  uint32_t bytes[32];
  int nbytes = 1, bits = 0;
  bytes[0] = 1;
  for (int i = 0; i < count; ++i) {
    uint32_t carry = 0;
    for (int j = 0; j < nbytes; ++j) {
      uint32_t tmp = bytes[j] * sizes[i] + carry;
      bytes[j] = tmp & 0xff;
      carry = tmp >> 8;
    }
    while (carry) { bytes[nbytes++] = carry & 0xff; carry >>= 8; }
  }
  int num = 1;
  --nbytes;
  while ((int)bytes[nbytes] >= num) { bits++; num *= 2; }
  return bits + nbytes * 8;
}

// ── Bit reader (for decompression) ───────────────────────────────────────
//
// The XTC bit stream is big-endian: the first bit of byte 0 is the MSB.
// Multi-byte integers are assembled little-endian (first byte = LSB)
// to match the mixed-radix encoding used by GROMACS xdrfile.

class BitReader {
public:
  explicit BitReader(const uint8_t *data) : m_data(data), m_pos(0) {}

  // Read nbits from the stream MSB-first; result bit 0 = last bit read.
  uint32_t read(int nbits) {
    uint32_t r = 0;
    for (int i = 0; i < nbits; ++i)
      r = (r << 1) | next_bit();
    return r;
  }

  // Read bitsize bits and decode as a mixed-radix triple using sizeint[3].
  // Bytes are assembled little-endian (LSB first) then mixed-radix factored.
  void unpack(int bitsize, const int32_t sizeint[3], int32_t out[3]) {
    uint64_t v = 0;
    int fullbytes = bitsize / 8;
    int partbits  = bitsize % 8;
    for (int i = 0; i < fullbytes; ++i)
      v |= static_cast<uint64_t>(read(8)) << (8 * i);
    if (partbits)
      v |= static_cast<uint64_t>(read(partbits)) << (8 * fullbytes);

    uint64_t sz  = static_cast<uint64_t>(sizeint[2]);
    uint64_t szy = sz * static_cast<uint64_t>(sizeint[1]);
    out[0] = static_cast<int32_t>(v / szy);
    uint64_t q = v - static_cast<uint64_t>(out[0]) * szy;
    out[1] = static_cast<int32_t>(q / sz);
    out[2] = static_cast<int32_t>(q - static_cast<uint64_t>(out[1]) * sz);
  }

private:
  int next_bit() {
    int byte_idx    = m_pos / 8;
    int bit_in_byte = 7 - (m_pos % 8); // MSB first within each byte
    ++m_pos;
    return (m_data[byte_idx] >> bit_in_byte) & 1;
  }

  const uint8_t *m_data;
  int m_pos;
};

// ── Decompression ─────────────────────────────────────────────────────────

static bool xtc_decompress(uint32_t natoms, float inv_p,
                           const int32_t minint[3], const int32_t maxint[3],
                           int32_t smallidx, const uint8_t *data,
                           float *out) {
  int smlim = smallidx - 1;
  if (smlim < FIRSTIDX) smlim = FIRSTIDX;

  uint32_t sizeint[3];
  int bitsizeint[3] = {0, 0, 0};
  int bitsize;
  bool large;

  for (int d = 0; d < 3; ++d) sizeint[d] = static_cast<uint32_t>(maxint[d] - minint[d] + 1);

  if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffffu) {
    for (int d = 0; d < 3; ++d) bitsizeint[d] = sizeofint(sizeint[d]);
    bitsize = 0;
    large   = true;
  } else {
    bitsize = sizeofints(3, sizeint);
    large   = false;
  }

  BitReader br(data);
  int smaller  = magicints[smlim] / 2;
  int smallnum = magicints[smallidx] / 2;
  int sizesmall[3] = {magicints[smallidx], magicints[smallidx], magicints[smallidx]};

  int32_t thiscoord[3], prevcoord[3] = {0, 0, 0}, tmp[3];
  int32_t isizeint[3] = {static_cast<int32_t>(sizeint[0]),
                         static_cast<int32_t>(sizeint[1]),
                         static_cast<int32_t>(sizeint[2])};

  float *lfp = out;
  const float *const end = out + natoms * 3;
  int i = 0;

  while (i < static_cast<int>(natoms)) {
    // ── Decode base atom ─────────────────────────────────────────────────
    if (!large) {
      br.unpack(bitsize, isizeint, thiscoord);
    } else {
      for (int d = 0; d < 3; ++d)
        thiscoord[d] = static_cast<int32_t>(br.read(bitsizeint[d]));
    }
    for (int d = 0; d < 3; ++d) {
      thiscoord[d] += minint[d];
      prevcoord[d]  = thiscoord[d];
    }
    ++i;

    // ── Flag + optional run-length block ─────────────────────────────────
    int flag = br.read(1);
    int is_smaller = 0, run = 0;
    if (flag) {
      run        = static_cast<int>(br.read(5));
      is_smaller = run % 3;
      run       -= is_smaller;
      is_smaller--;
    }

    // Debug first few and near-failure iterations
    if (i <= 5 || i >= 655)
      fprintf(stderr, "iter i=%d flag=%d run=%d is_smaller=%d smallidx=%d sizesmall=%d\n",
              i, flag, run, is_smaller, smallidx, sizesmall[0]);

    // Buffer overrun check (matches GROMACS xdrfile.c)
    if ((lfp - out) + run > static_cast<int>(natoms) * 3) {
      fprintf(stderr, "XTC overrun: lfp_off=%td run=%d natoms*3=%d\n",
              lfp - out, run, static_cast<int>(natoms) * 3);
      return false;
    }

    if (run > 0) {
      // Decode small (delta-encoded) atoms
      int32_t szs[3] = {sizesmall[0], sizesmall[1], sizesmall[2]};
      for (int k = 0; k < run; k += 3) {
        br.unpack(smallidx, szs, thiscoord);
        ++i;
        for (int d = 0; d < 3; ++d)
          thiscoord[d] += prevcoord[d] - smallnum;

        if (k == 0) {
          // Interchange first small atom with the base atom (water optimisation)
          std::memcpy(tmp,       thiscoord, 12);
          std::memcpy(thiscoord, prevcoord, 12);
          std::memcpy(prevcoord, tmp,       12);
          *lfp++ = prevcoord[0] * inv_p;
          *lfp++ = prevcoord[1] * inv_p;
          *lfp++ = prevcoord[2] * inv_p;
        } else {
          prevcoord[0] = thiscoord[0];
          prevcoord[1] = thiscoord[1];
          prevcoord[2] = thiscoord[2];
        }
        *lfp++ = thiscoord[0] * inv_p;
        *lfp++ = thiscoord[1] * inv_p;
        *lfp++ = thiscoord[2] * inv_p;
      }
    } else {
      *lfp++ = thiscoord[0] * inv_p;
      *lfp++ = thiscoord[1] * inv_p;
      *lfp++ = thiscoord[2] * inv_p;
    }

    smallidx += is_smaller;
    if (is_smaller < 0) {
      smallnum = smaller;
      smaller  = (smallidx > FIRSTIDX) ? magicints[smallidx - 1] / 2 : 0;
    } else if (is_smaller > 0) {
      smaller    = smallnum;
      smallnum   = magicints[smallidx] / 2;
    }
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
    if (smallidx < FIRSTIDX || smallidx >= LASTIDX) {
      fprintf(stderr, "XTC: smallidx=%d out of range [%d,%d) at i=%d is_smaller=%d\n",
              smallidx, FIRSTIDX, LASTIDX, i, is_smaller);
      return false;
    }
    if (sizesmall[0] == 0) { fprintf(stderr, "XTC: sizesmall==0 at i=%d smallidx=%d\n", i, smallidx); return false; }
  }
  if (lfp != end) {
    fprintf(stderr, "XTC: lfp mismatch: wrote %td floats, expected %u\n", lfp - out, natoms * 3);
    return false;
  }
  return true;
}

// ── Bit writer (for compression) ─────────────────────────────────────────

class BitWriter {
public:
  // Write nbits of value MSB-first into the stream.
  void write(uint32_t value, int nbits) {
    for (int i = nbits - 1; i >= 0; --i) {
      m_cur = (m_cur << 1) | ((value >> i) & 1);
      if (++m_bits == 8) { m_buf.push_back(m_cur); m_cur = 0; m_bits = 0; }
    }
  }

  // Encode a mixed-radix triple and write as bitsize-bit little-endian integer.
  void encode(int bitsize, const uint32_t sizes[3], const int32_t vals[3]) {
    // Mixed-radix encode: v = vals[0]*sizes[1]*sizes[2] + vals[1]*sizes[2] + vals[2]
    // Store as bytes LSB first, each byte written MSB-first into stream.
    uint8_t bytes[8] = {};
    int nbytes = 1;
    bytes[0] = static_cast<uint8_t>(vals[0]) & 0xff;
    bytes[1] = static_cast<uint8_t>(vals[0] >> 8) & 0xff;
    bytes[2] = static_cast<uint8_t>(vals[0] >> 16) & 0xff;
    bytes[3] = static_cast<uint8_t>(vals[0] >> 24) & 0xff;
    nbytes   = 4;

    for (int i = 1; i < 3; ++i) {
      uint32_t carry = static_cast<uint32_t>(vals[i]);
      for (int j = 0; j < nbytes; ++j) {
        uint32_t tmp = bytes[j] * sizes[i] + carry;
        bytes[j]  = tmp & 0xff;
        carry     = tmp >> 8;
      }
      while (carry) { bytes[nbytes++] = carry & 0xff; carry >>= 8; }
    }

    int bits_left = bitsize;
    for (int i = 0; i < nbytes && bits_left > 0; ++i) {
      int n = std::min(8, bits_left);
      write(bytes[i], n);
      bits_left -= n;
    }
  }

  std::vector<uint8_t> finish() {
    if (m_bits > 0) {
      m_buf.push_back(static_cast<uint8_t>(m_cur << (8 - m_bits)));
      m_cur  = 0;
      m_bits = 0;
    }
    return std::move(m_buf);
  }

private:
  std::vector<uint8_t> m_buf;
  uint8_t m_cur{0};
  int     m_bits{0};
};

// ── Compression ───────────────────────────────────────────────────────────

static void xtc_compress(uint32_t natoms, float precision, const float *coords,
                         std::vector<uint8_t> &out, int32_t &nbytes,
                         int32_t minint[3], int32_t maxint[3], int &smallidx) {
  // Quantise to integers
  std::vector<int32_t> ic(natoms * 3);
  for (int d = 0; d < 3; ++d) {
    minint[d] = std::numeric_limits<int32_t>::max();
    maxint[d] = std::numeric_limits<int32_t>::min();
  }
  for (uint32_t i = 0; i < natoms; ++i) {
    for (int d = 0; d < 3; ++d) {
      int32_t iv = static_cast<int32_t>(std::round(coords[i * 3 + d] * precision));
      ic[i * 3 + d] = iv;
      minint[d] = std::min(minint[d], iv);
      maxint[d] = std::max(maxint[d], iv);
    }
  }

  uint32_t sizeint[3];
  for (int d = 0; d < 3; ++d) sizeint[d] = maxint[d] - minint[d] + 1;

  bool large = (sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffffu;
  int bitsize   = 0;
  int bitsizeint[3] = {0, 0, 0};
  if (large) {
    for (int d = 0; d < 3; ++d) bitsizeint[d] = sizeofint(sizeint[d]);
  } else {
    bitsize = sizeofints(3, sizeint);
  }

  // Determine smallidx (smallest magicints value >= all sizeint components)
  smallidx = FIRSTIDX;
  while (smallidx < LASTIDX - 1 &&
         magicints[smallidx] < static_cast<int>(sizeint[0]) &&
         magicints[smallidx] < static_cast<int>(sizeint[1]) &&
         magicints[smallidx] < static_cast<int>(sizeint[2]))
    ++smallidx;

  BitWriter bw;
  for (uint32_t i = 0; i < natoms; ++i) {
    int32_t rel[3] = {
        ic[i * 3 + 0] - minint[0],
        ic[i * 3 + 1] - minint[1],
        ic[i * 3 + 2] - minint[2],
    };

    if (large) {
      for (int d = 0; d < 3; ++d)
        bw.write(static_cast<uint32_t>(rel[d]), bitsizeint[d]);
    } else {
      bw.encode(bitsize, sizeint, rel);
    }
    bw.write(0, 1); // flag = 0: no small-coord run follows
  }

  out    = bw.finish();
  nbytes = static_cast<int32_t>(out.size());
}

// ── XTCHandler ────────────────────────────────────────────────────────────

bool XTCHandler::_initialise() {
  if (m_mode == Mode::Read) {
    m_infile.open(this->file_path(), std::ios::binary);
    if (!m_infile.is_open()) return false;
    m_infile.seekg(0, std::ios::end);
    m_file_size = m_infile.tellg();
    m_infile.seekg(0, std::ios::beg);
    return true;
  } else {
    m_outfile.open(this->file_path(), std::ios::binary | std::ios::trunc);
    return m_outfile.is_open();
  }
}

void XTCHandler::_finalise() {
  if (m_infile.is_open())  m_infile.close();
  if (m_outfile.is_open()) { m_outfile.flush(); m_outfile.close(); }
}

bool XTCHandler::read_next_frame(Frame &frame) {
  if (!m_infile.good() || m_infile.tellg() >= m_file_size) return false;

  // ── Header ────────────────────────────────────────────────────────────
  int32_t magic = rd_i32(m_infile);
  if (!m_infile.good()) return false;
  if (magic != XTC_MAGIC)
    throw std::runtime_error("XTC: bad magic number in frame header");

  m_natoms = rd_u32(m_infile);
  /*uint32_t step =*/ rd_u32(m_infile);
  /*float time    =*/ rd_f32(m_infile);

  // Box matrix (row-major, nm)
  float box_nm[3][3];
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      box_nm[r][c] = rd_f32(m_infile);

  {
    // Columns of the direct matrix are the lattice vectors (transposed from row-major)
    occ::Mat3 m;
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c)
        m(c, r) = static_cast<double>(box_nm[r][c]) * 10.0; // nm → Å

    occ::Vec3 a = m.col(0), b = m.col(1), c = m.col(2);
    double na = a.norm(), nb = b.norm(), nc = c.norm();
    double ca = std::clamp(b.dot(c) / (nb * nc), -1.0, 1.0);
    double cb = std::clamp(a.dot(c) / (na * nc), -1.0, 1.0);
    double cg = std::clamp(a.dot(b) / (na * nb), -1.0, 1.0);
    double alpha = std::acos(ca), beta = std::acos(cb), gamma = std::acos(cg);
    if (trajan::util::unitcell_is_reasonable(na, nb, nc, alpha, beta, gamma))
      frame.set_unit_cell(occ::crystal::triclinic_cell(na, nb, nc, alpha, beta, gamma));
  }

  uint32_t lsize = rd_u32(m_infile);

  m_coords.resize(m_natoms * 3);

  if (lsize <= 9u) {
    // ── Uncompressed ────────────────────────────────────────────────────
    for (uint32_t i = 0; i < m_natoms * 3u; ++i)
      m_coords[i] = rd_f32(m_infile);
  } else {
    // ── Compressed ──────────────────────────────────────────────────────
    m_precision   = rd_f32(m_infile);
    int32_t minint[3], maxint[3];
    for (int d = 0; d < 3; ++d) minint[d] = rd_i32(m_infile);
    for (int d = 0; d < 3; ++d) maxint[d] = rd_i32(m_infile);
    int32_t smallidx = rd_i32(m_infile);
    int32_t nbytes   = rd_i32(m_infile);

    int buf_size = (nbytes + 7) & ~7; // pad to 8 bytes for bit reader
    m_buf.assign(buf_size, 0);
    m_infile.read(reinterpret_cast<char *>(m_buf.data()), align4(nbytes));

    float inv_p = 1.0f / m_precision;
    if (!xtc_decompress(m_natoms, inv_p, minint, maxint, smallidx,
                        m_buf.data(), m_coords.data()))
      throw std::runtime_error("XTC: coordinate decompression failed");
  }

  // Copy nm → Å into frame
  for (uint32_t i = 0; i < m_natoms; ++i) {
    occ::Vec3 pos = {
        static_cast<double>(m_coords[i * 3 + 0]) * 10.0,
        static_cast<double>(m_coords[i * 3 + 1]) * 10.0,
        static_cast<double>(m_coords[i * 3 + 2]) * 10.0,
    };
    frame.update_atom_position(i, pos);
  }

  if (!m_infile.good() && !m_infile.eof())
    throw std::runtime_error("XTC: read error in frame body");
  return true;
}

bool XTCHandler::write_next_frame(const Frame &frame) {
  uint32_t natoms = static_cast<uint32_t>(frame.num_atoms());

  // ── Header ────────────────────────────────────────────────────────────
  wr_i32(m_outfile, XTC_MAGIC);
  wr_u32(m_outfile, natoms);
  wr_u32(m_outfile, m_step++);
  wr_f32(m_outfile, 0.0f); // time (ps) — not tracked

  // Box matrix (row-major, nm)
  occ::Mat3 direct = occ::Mat3::Zero();
  if (frame.has_unit_cell())
    direct = frame.unit_cell().value().direct();
  // direct columns = lattice vectors; write as rows (transposed)
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      wr_f32(m_outfile, static_cast<float>(direct(c, r) / 10.0)); // Å → nm

  // ── Coordinates → nm ─────────────────────────────────────────────────
  const auto &atoms = frame.atoms();
  m_coords.resize(natoms * 3);
  for (uint32_t i = 0; i < natoms; ++i) {
    m_coords[i * 3 + 0] = static_cast<float>(atoms[i].x / 10.0);
    m_coords[i * 3 + 1] = static_cast<float>(atoms[i].y / 10.0);
    m_coords[i * 3 + 2] = static_cast<float>(atoms[i].z / 10.0);
  }

  if (natoms <= 9u) {
    // ── Uncompressed (small systems) ─────────────────────────────────────
    wr_u32(m_outfile, natoms);
    for (uint32_t i = 0; i < natoms * 3u; ++i)
      wr_f32(m_outfile, m_coords[i]);
  } else {
    // ── Compressed ───────────────────────────────────────────────────────
    std::vector<uint8_t> comp;
    int32_t minint[3], maxint[3], nbytes;
    int smallidx;
    xtc_compress(natoms, m_precision, m_coords.data(),
                 comp, nbytes, minint, maxint, smallidx);

    uint32_t lsize = static_cast<uint32_t>(36 + align4(nbytes)); // hdr2 + data
    wr_u32(m_outfile, lsize);

    wr_f32(m_outfile, m_precision);
    for (int d = 0; d < 3; ++d) wr_i32(m_outfile, minint[d]);
    for (int d = 0; d < 3; ++d) wr_i32(m_outfile, maxint[d]);
    wr_i32(m_outfile, static_cast<int32_t>(smallidx));
    wr_i32(m_outfile, nbytes);

    m_outfile.write(reinterpret_cast<const char *>(comp.data()), nbytes);
    // Pad to 4-byte alignment
    static const uint8_t zeros[4] = {};
    int pad = align4(nbytes) - nbytes;
    if (pad > 0) m_outfile.write(reinterpret_cast<const char *>(zeros), pad);
  }

  return m_outfile.good();
}

} // namespace trajan::io
