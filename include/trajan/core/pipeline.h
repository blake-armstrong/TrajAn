#pragma once
#include <functional>
#include <trajan/core/frame.h>
#include <vector>

namespace trajan::core {

using FrameTransform = std::function<void(Frame &)>;

class Pipeline {
public:
  void add_transform(FrameTransform fn) {
    m_transforms.push_back(std::move(fn));
  }

  void apply(Frame &frame) const {
    for (const auto &fn : m_transforms)
      fn(frame);
  }

  bool empty() const { return m_transforms.empty(); }
  size_t size() const { return m_transforms.size(); }

private:
  std::vector<FrameTransform> m_transforms;
};

} // namespace trajan::core
