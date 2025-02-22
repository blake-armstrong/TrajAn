#include "trajan/core/frame.h"
#include "trajan/core/trajectory.h"
#include "trajan/io/file.h"
#include <fstream>
#include <stdexcept>

namespace trajan::io {

namespace core = trajan::core;

bool DCDHandler::read_frame(core::Frame &frame) { return true; }
}; // namespace trajan::io
