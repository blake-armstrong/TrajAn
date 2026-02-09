#include <CLI/CLI.hpp>
#include <occ/core/bondgraph.h>
#include <optional>
#include <stdexcept>
#include <trajan/core/neigh.h>
#include <trajan/core/trajectory.h>
#include <trajan/io/selection.h>
#include <trajan/main/trajan_topology.h>

namespace trajan::main {

using trajan::core::Atom;

BondCriteria
bond_criteria_validator(const std::string &input,
                        std::optional<std::vector<char>> restrictions) {
  // parse format: "sel1&sel2</>num"
  trajan::log::debug("Attempting to parse bond criteria: '{}'", input);
  size_t op_pos = input.find_first_of("<>");
  if (op_pos == std::string::npos) {
    throw std::invalid_argument("Missing comparison operator (< or >)");
  }
  const char op_char = input[op_pos];
  std::string left = input.substr(0, op_pos);
  std::string right = input.substr(op_pos + 1);

  size_t and_pos = left.find('&');
  if (and_pos == std::string::npos) {
    throw std::invalid_argument("Missing comma separator between selections");
  }
  std::string sel1_str = left.substr(0, and_pos);
  std::string sel2_str = left.substr(and_pos + 1);

  double threshold;
  try {
    threshold = std::stod(right);
  } catch (const std::exception &e) {
    throw std::invalid_argument(
        fmt::format("Invalid threshold value: {} ({})", right, e.what()));
  }

  BondCriteria bc;

  try {
    bc.sel1 = io::selection_validator(sel1_str, restrictions);
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(fmt::format(
        "Error parsing bond criteria input '{}': {} ", sel1_str, e.what()));
  }
  try {
    bc.sel2 = io::selection_validator(sel2_str, restrictions);
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(fmt::format(
        "Error parsing bond criteria input '{}': {} ", sel2_str, e.what()));
  }

  bc.threshold = threshold;
  bc.op = (op_char == '<') ? BondCriteria::ComparisonOp::LessThan
                           : BondCriteria::ComparisonOp::GreaterThan;
  return bc;
}

void run_topology_subcommand(const TopologyOpts &opts, Trajectory &traj) {
  core::TopologyUpdateSettings settings;
  for (const std::string &bc_str : opts.bc_raw_sel) {
    BondCriteria bond_criteria =
        bond_criteria_validator(bc_str, MOLECULE_RESTRICTIONS);
    core::BondCutoff bond_cutoff;
    auto atoms1 = traj.get_atoms(bond_criteria.sel1);
    bond_cutoff.atom_indices1.reserve(atoms1.size());
    for (const auto &a : atoms1) {
      bond_cutoff.atom_indices1.push_back(a.index);
    }
    std::sort(bond_cutoff.atom_indices1.begin(),
              bond_cutoff.atom_indices1.end());
    auto atoms2 = traj.get_atoms(bond_criteria.sel2);
    bond_cutoff.atom_indices2.reserve(atoms2.size());
    for (const auto &a : atoms2) {
      bond_cutoff.atom_indices2.push_back(a.index);
    }
    std::sort(bond_cutoff.atom_indices2.begin(),
              bond_cutoff.atom_indices2.end());
    bond_cutoff.op = (bond_criteria.op == BondCriteria::ComparisonOp::LessThan)
                         ? core::BondCutoff::ComparisonOp::LessThan
                         : core::BondCutoff::ComparisonOp::GreaterThan;
    bond_cutoff.threshold = bond_criteria.threshold;
    settings.bond_cutoffs.push_back(bond_cutoff);
  }
  std::vector<io::SelectionCriteria> nb_parsed_sel;
  for (const std::string &input : opts.nb_raw_sel) {
    auto nb_parsed_sels = io::selection_validator(input, MOLECULE_RESTRICTIONS);
    for (const auto &sel : nb_parsed_sels) {
      nb_parsed_sel.push_back(sel);
    }
  }
  auto nbatoms = traj.get_atoms(nb_parsed_sel);
  settings.no_bonds.reserve(nbatoms.size());
  for (const auto &a : nbatoms) {
    settings.no_bonds.push_back(a.index);
  }
  std::sort(settings.no_bonds.begin(), settings.no_bonds.end());
  settings.top_auto = opts.top_auto;
  settings.bond_tolerance = opts.bond_tolerance;
  settings.update_frequency = opts.update_frequency;
  settings.compute_topology = true;
  occ::core::set_bond_tolerance(opts.bond_tolerance);
  traj.set_topology_settings(settings);
}

CLI::App *add_topology_subcommand(CLI::App &app, Trajectory &traj) {
  CLI::App *top = app.add_subcommand("top", "Compute the molecular topology.");
  auto opts = std::make_shared<TopologyOpts>();
  top->add_flag("-a,--auto", opts->top_auto,
                "Generate topology automatically from atom covalent radii.");
  top->add_option("--no-bond", opts->nb_raw_sel,
                  "Selected atoms will not be allowed to bond to other atoms\n"
                  "First selection (prefix: i=atom indices, a=atom types)\n"
                  "Examples:\n"
                  "  i1,2,3-5    (atom indices 1,2,3,4,5)\n"
                  "  aC,N,O      (atom types C, N, O)\n"
                  "  j/m not allowed\n");
  // ->check(io::selection_validator(
  //     opts->nb_parsed_sel,
  //     std::make_optional<std::vector<char>>({'j', 'm'})));
  top->add_option(
         "--bond-criteria", opts->bc_raw_sel,
         "Selected pairs of atoms will only be allowed to form a bond if the "
         "distance criteria is less than or greater than the input "
         "threshold.\n Selection one and two should be separated by '&'.\n "
         "Examples:\n  i1,2,3-5&aC,N,O<1.4 (bonds only formed between atom "
         "indices 1,2,3,4,5 and atom types C,N,O when separated by less than "
         "1.4 Angstroms.)\n  i1,2&aN,H>2.5\n")
      ->type_size(0, -1);
  // ->each(bond_criteria_validator(
  //     opts->bond_criterias,
  //     std::make_optional<std::vector<char>>({'j', 'm'})));
  opts->bond_tolerance = occ::core::get_bond_tolerance();
  top->add_option("--bond-tolerance", opts->bond_tolerance,
                  "Bond tolerance in Angstroms to use when deciding if two "
                  "atoms are bonded. This number is added to the sum of the "
                  "covalent radii.\n");
  top->add_option(
      "--update-frequency", opts->update_frequency,
      "How often to recompute the topology. Default is set to 0 which means to "
      "only compute the topology once at the start of the trajectory analysis. "
      "This means the topology is assumed to unchage throughout. Setting this "
      "to anything other than 0 will slow down the trajectory analysis, but it "
      "can be useful when analysing non-classical or reactive-classical "
      "simulations.\n");
  top->callback([opts, &traj]() {
    trajan::log::set_subcommand_log_pattern("top");
    run_topology_subcommand(*opts, traj);
  });
  return top;
}

} // namespace trajan::main
