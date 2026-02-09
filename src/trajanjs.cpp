#include "js/core_bindings.h"
#include <emscripten/bind.h>
#include <emscripten/val.h>

using namespace emscripten;

EMSCRIPTEN_BINDINGS(trajan) {

  register_core_bindings();

  function("setLogFile", &occ::log::set_log_file);

  constant("version", std::string("0.0.0"));
}
