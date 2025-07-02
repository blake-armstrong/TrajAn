#pragma once
#include <nanobind/nanobind.h>

namespace nb = nanobind;

void register_trajectory_bindings(nb::module_ &m);
