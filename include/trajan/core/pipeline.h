#pragma once
#include <functional>
#include <string>
#include <trajan/core/frame.h>
#include <trajan/core/log.h>
#include <vector>

namespace trajan::core {

using FrameTransform = std::function<void(Frame &)>;

class Pipeline {
public:
  void add_transform(std::string subcommand, FrameTransform fn) {
    m_transforms.push_back({std::move(subcommand), std::move(fn)});
  }

  void apply(Frame &frame, const std::string &caller) const {
    for (const auto &[name, fn] : m_transforms) {
      trajan::log::set_subcommand_log_pattern(name);
      fn(frame);
    }
    trajan::log::set_subcommand_log_pattern(caller);
  }

  bool empty() const { return m_transforms.empty(); }
  size_t size() const { return m_transforms.size(); }

private:
  struct NamedTransform {
    std::string subcommand;
    FrameTransform fn;
  };
  std::vector<NamedTransform> m_transforms;
};

} // namespace trajan::core
