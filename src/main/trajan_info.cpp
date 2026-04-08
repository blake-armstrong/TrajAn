#include <map>
#include <trajan/core/log.h>
#include <trajan/main/trajan_info.h>

namespace trajan::main {

namespace {

// Separator width
constexpr size_t SEP_WIDTH = 80;
constexpr char SEP_HEAVY = '=';
constexpr char SEP_LIGHT = '-';
constexpr char FILL_CHAR = '.';

std::string make_separator(char c, size_t width = SEP_WIDTH) {
  return std::string(width, c);
}

std::string make_fill(int count) {
  return count > 0 ? std::string(count, FILL_CHAR) : "";
}

// Prints:  "  |  label ............: value"
// `fill_width` is the total width of the label+fill column (before the colon).
void print_field(size_t fill_width, const std::string &label,
                 const std::string &value,
                 const std::string &indent = "  |  ") {
  const int fill_count =
      static_cast<int>(fill_width) - static_cast<int>(label.size()) - 1;
  trajan::log::info("{}{} {} : {}", indent, label, make_fill(fill_count),
                    value);
}

void print_section_header(const std::string &title) {
  trajan::log::info("");
  // trajan::log::info("{}", make_separator(SEP_HEAVY));
  const std::string prefix = "  +- ";
  const size_t remaining = SEP_WIDTH - prefix.size() - title.size() - 1;
  trajan::log::info("{}{} {}", prefix, title,
                    make_separator(SEP_LIGHT, remaining));
}

// Prints a field whose value is a list of strings, wrapping at SEP_WIDTH.
// Continuation lines use an empty label (same fill width) so the value column
// stays aligned.
void print_wrapped_field(size_t fill_width, const std::string &label,
                         std::vector<std::string> items,
                         const std::string &indent = "  |  ") {
  // prefix width: indent + fill_width chars + " : " (3)
  const size_t prefix_width = indent.size() + fill_width + 1 + 3;
  const size_t value_width =
      SEP_WIDTH > prefix_width ? SEP_WIDTH - prefix_width : 10;

  std::vector<std::string> lines;
  std::string current;
  for (size_t i = 0; i < items.size(); ++i) {
    const std::string sep = current.empty() ? "" : ", ";
    if (!current.empty() &&
        current.size() + sep.size() + items[i].size() > value_width) {
      lines.push_back(current);
      current = items[i];
    } else {
      current += sep + items[i];
    }
  }
  if (!current.empty())
    lines.push_back(current);
  if (lines.empty())
    lines.push_back("None");

  print_field(fill_width, label, lines[0], indent);
  for (size_t i = 1; i < lines.size(); ++i)
    print_field(fill_width, "", lines[i], indent);
}

void print_section_footer() {
  trajan::log::info("  +{}", make_separator(SEP_LIGHT, SEP_WIDTH - 3));
}

void print_subsection(const std::string &title) {
  trajan::log::info("  |");
  const std::string prefix = "  |  ";
  const size_t remaining = SEP_WIDTH - prefix.size() - title.size() - 1;
  trajan::log::info("{}{} {}", prefix, title,
                    make_separator(SEP_LIGHT, remaining));
}

struct UnitCellAccum {
  double a{0}, b{0}, c{0}, alpha{0}, beta{0}, gamma{0}, volume{0};
  size_t count{0};

  void add(const occ::crystal::UnitCell &uc) {
    a += uc.a(); b += uc.b(); c += uc.c();
    alpha += uc.alpha(); beta += uc.beta(); gamma += uc.gamma();
    volume += uc.volume();
    ++count;
  }
};

} // namespace

void run_info_subcommand(const InfoOpts &opts, Trajectory &traj) {

  UnitCellAccum uc_accum;

  // Capture snapshot data from first frame, accumulate unit cell over all.
  std::optional<occ::crystal::UnitCell> last_uc;
  std::map<std::string, size_t> atom_type_count;
  std::map<std::string, std::string> atom_type_elements;
  double mass = 0.0;
  size_t total_atoms = 0;
  size_t frame_count = 0;

  {
    trajan::log::Progress progress("info: reading frames");
    while (traj.next_frame()) {
      ++frame_count;
      progress.increment();

      const auto uc = traj.unit_cell();
      if (uc) {
        uc_accum.add(uc.value());
        last_uc = uc;
      }

      // Composition and topology only need one frame — use the first.
      if (frame_count == 1) {
        atom_type_count.clear();
        atom_type_elements.clear();
        mass = 0.0;
        total_atoms = traj.atoms().size();
        for (const auto &a : traj.atoms()) {
          atom_type_count[a.type]++;
          atom_type_elements[a.type] = a.element.symbol();
          mass += a.element.mass();
        }
      }
    }
    progress.finish(fmt::format("info: processed {} frame(s)", frame_count));
  }

  if (frame_count == 0) {
    trajan::log::warn("info: no frames available");
    return;
  }

  // -- Unit Cell (averaged) -------------------------------------------------
  print_section_header(frame_count > 1 ? "Unit Cell (averaged over all frames)"
                                       : "Unit Cell");
  if (uc_accum.count > 0) {
    const double n = static_cast<double>(uc_accum.count);
    const double a     = uc_accum.a / n;
    const double b     = uc_accum.b / n;
    const double c     = uc_accum.c / n;
    const double alpha = uc_accum.alpha / n;
    const double beta  = uc_accum.beta  / n;
    const double gamma = uc_accum.gamma / n;
    const double vol   = uc_accum.volume / n;

    // Rebuild an averaged unit cell for derived quantities.
    const occ::crystal::UnitCell avg_uc =
        occ::crystal::triclinic_cell(a, b, c, alpha, beta, gamma);

    constexpr size_t fw = 32;
    using occ::units::degrees;
    print_field(fw, "lengths (a, b, c)",
                fmt::format("{:10.4f}  {:10.4f}  {:10.4f}", a, b, c));
    print_field(fw, "angles (alpha, beta, gamma)",
                fmt::format("{:10.4f}  {:10.4f}  {:10.4f}",
                            degrees(alpha), degrees(beta), degrees(gamma)));
    print_field(fw, "cell type", fmt::format("{}", avg_uc.cell_type()));
    print_field(fw, "volume", fmt::format("{:.4f} A^3", vol));

    double area_bc = b * c * sin(alpha);
    double area_ac = a * c * sin(beta);
    double area_ab = a * b * sin(gamma);
    print_field(fw, "surface area",
                fmt::format("{:.4f} A^2", 2.0 * (area_bc + area_ac + area_ab)));

    std::array<std::string, 3> labels{"a", "b", "c"};

    print_subsection("Direct lattice vectors");
    const auto &dt = avg_uc.direct();
    for (size_t i = 0; i < 3; i++) {
      const auto v = dt.col(i);
      print_field(fw, fmt::format("({})", labels[i]),
                  fmt::format("{:10.4f}  {:10.4f}  {:10.4f}", v(0), v(1), v(2)));
    }

    print_subsection("Reciprocal lattice vectors");
    const auto &rt = avg_uc.reciprocal();
    for (size_t i = 0; i < 3; i++) {
      const auto v = rt.col(i);
      print_field(fw, fmt::format("({}*)", labels[i]),
                  fmt::format("{:10.4f}  {:10.4f}  {:10.4f}", v(0), v(1), v(2)));
    }
  } else {
    trajan::log::info("  |  No unit cell information available.");
  }
  print_section_footer();

  // -- Composition (from first frame) ---------------------------------------
  print_section_header("Composition");
  constexpr size_t fw = 32;
  print_field(fw, "total atoms", fmt::format("{}", total_atoms));
  print_field(fw, "total mass", fmt::format("{:.4f} g/mol", mass));
  if (uc_accum.count > 0) {
    const double avg_vol = uc_accum.volume / static_cast<double>(uc_accum.count);
    print_field(fw, "mass density",
                fmt::format("{:.4f} g/cm^3", mass / avg_vol * 1.6605388));
    print_field(fw, "number density",
                fmt::format("{:.4f} atom/A^3", total_atoms / avg_vol));
  }

  print_subsection("Atom types");
  trajan::log::info("  |  {:>8}  {:>4}  {:>10}", "type", "elem", "count");
  trajan::log::info("  |  {:>8}  {:>4}  {:>10}", "--------", "----",
                    "----------");
  for (const auto &[at, count] : atom_type_count) {
    trajan::log::info("  |  {:>8}  {:>4}  {:>10}", at, atom_type_elements[at],
                      count);
  }
  print_section_footer();

  // -- Topology (from last frame) -------------------------------------------
  print_section_header("Topology");
  core::Topology &top = traj.get_topology();
  constexpr size_t tfw = 32;
  print_field(tfw, "atoms",              fmt::format("{}", top.num_atoms()));
  print_field(tfw, "bonds",              fmt::format("{}", top.num_bonds()));
  print_field(tfw, "angles",             fmt::format("{}", top.num_angles()));
  print_field(tfw, "proper dihedrals",   fmt::format("{}", top.num_proper_dihedrals()));
  print_field(tfw, "improper dihedrals", fmt::format("{}", top.num_improper_dihedrals()));
  print_section_footer();

  // -- Molecules ------------------------------------------------------------
  auto &molecules = top.get_molecules();
  if (molecules.size() > 0) {
    print_section_header("Molecules");

    ankerl::unordered_dense::map<std::string, size_t> unique_molecules;
    ankerl::unordered_dense::map<std::string, size_t> unique_molecule_counts;

    for (const auto &m : molecules) {
      trajan::log::debug("  {}", m.repr());
      unique_molecules[m.type] = m.index;
      unique_molecule_counts[m.type]++;
    }

    const double avg_vol = uc_accum.count > 0
                               ? uc_accum.volume / static_cast<double>(uc_accum.count)
                               : 0.0;
    constexpr size_t mfw = 42;

    bool first = true;
    for (const auto &[t, i] : unique_molecules) {
      const auto &m = molecules[i];
      const size_t count = unique_molecule_counts[t];

      if (!first)
        trajan::log::info("  |");
      first = false;

      const std::string mol_prefix = "  |  +- Molecule type: ";
      const size_t mol_remaining = SEP_WIDTH - mol_prefix.size() - t.size() - 1;
      trajan::log::info("{}{} {}", mol_prefix, t,
                        make_separator(SEP_LIGHT, mol_remaining));

      auto mfield = [&](const std::string &label, const std::string &value) {
        print_field(mfw, label, value, "  |  |  ");
      };

      mfield("count", fmt::format("{}", count));
      mfield("atoms per molecule", fmt::format("{}", m.atoms().size()));

      auto atom_types_vec = m.atom_types();
      std::sort(atom_types_vec.begin(), atom_types_vec.end());
      print_wrapped_field(mfw, "atom types", atom_types_vec, "  |  |  ");

      auto elements_vec = m.element_symbols();
      std::sort(elements_vec.begin(), elements_vec.end());
      if (atom_types_vec != elements_vec)
        print_wrapped_field(mfw, "elements", elements_vec, "  |  |  ");

      mfield("molar mass", fmt::format("{:.4f} g/mol", m.molar_mass() * 1e3));

      if (avg_vol > 0.0) {
        mfield("concentration",
               fmt::format("{:.4f} mol/L", count / avg_vol * 1660.5388));
        mfield("mass fraction",
               fmt::format("{:.4f} %m/m", count * m.molar_mass() * 1e5 / mass));
      }

      const double T = 300;
      mfield(fmt::format("translational dG (@ {} K)", T),
             fmt::format("{:.4f} kJ/mol", m.translational_free_energy(T)));
      mfield(fmt::format("rotational dG (@ {} K)", T),
             fmt::format("{:.4f} kJ/mol", m.rotational_free_energy(T)));

      trajan::log::info("  |  +{}", make_separator(SEP_LIGHT, SEP_WIDTH - 7));
    }
    print_section_footer();
  }

  if (opts.detailed_top) {
    print_section_header("Detailed Topology");
    traj.frame().populate_angles(top);
    top.print_detailed();
    print_section_footer();
  }
};

CLI::App *add_info_subcommand(CLI::App &app, Trajectory &traj) {
  auto *info =
      app.add_subcommand("info", "Prints out a variety of information within "
                                 "the current Trajectory object.");
  auto opts = std::make_shared<InfoOpts>();
  info->add_flag("--detailed-topology", opts->detailed_top,
                 "Whether to output detailed topology info (if there is any "
                 "topology info).");
  info->add_option("--timings", opts->timings,
                   "Whether to output timing info.");

  info->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("info");
    run_info_subcommand(*opts, traj);
  });

  return info;
};

} // namespace trajan::main
