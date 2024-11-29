#include <unordered_map>
#include <complex>
#include <string>
#include <cmath>
#include <cstdint>
#include <vector>
#include <regex>
#include <thread>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

// Namespace for Pybind11
namespace py = pybind11;

// Type definition for the quantum state
using State = std::unordered_map<uint64_t, std::complex<double>>;

/**
 * Structure to represent an exponential Pauli term like "X(0)Y(1)Z(2)" and an angle θ.
 */
struct ExpPauliTerm {
    std::string pauli_string;
    double angle;

    // Constructors
    ExpPauliTerm() : pauli_string(""), angle(0.0) {}
    ExpPauliTerm(const std::string& p, double a) : pauli_string(p), angle(a) {}
};

/**
 * Parses a Pauli string like "X(0)X(2)" into a map from qubit indices to Pauli operators.
 *
 * @param pauli_string The Pauli string to parse, e.g., "X(0)Y(1)Z(2)".
 * @return A map from qubit indices to Pauli operators ('I', 'X', 'Y', 'Z').
 *
 * Example:
 *   auto pauli_map = parse_pauli_string("X(0)Y(1)Z(2)");
 *   // pauli_map[0] == 'X', pauli_map[1] == 'Y', pauli_map[2] == 'Z'
 */
std::unordered_map<int, char> parse_pauli_string(const std::string& pauli_string) {
    std::unordered_map<int, char> pauli_map;
    static const std::regex rgx("(I|X|Y|Z)\\((\\d+)\\)");
    std::smatch match;
    std::string s = pauli_string;

    if (pauli_string.empty()) {
        throw std::invalid_argument("Pauli string cannot be empty.");
    }

    while (std::regex_search(s, match, rgx)) {
        char op = match[1].str()[0];
        int qubit = std::stoi(match[2].str());
        if (qubit < 0) {
            throw std::invalid_argument("Qubit indices must be non-negative integers.");
        }
        pauli_map[qubit] = op;
        s = match.suffix().str();
    }

    if (pauli_map.empty() && pauli_string != "I") {
        throw std::invalid_argument("Invalid Pauli string format: " + pauli_string);
    }

    return pauli_map;
}

/**
 * Applies a Pauli operator to a basis state, returning the new basis state and the associated phase factor.
 *
 * @param pauli_map A map from qubit indices to Pauli operators ('I', 'X', 'Y', 'Z').
 * @param basis_state The basis state represented as an unsigned 64-bit integer.
 * @return A pair containing the new basis state and the phase factor (std::complex<double>).
 *
 * Example:
 *   auto pauli_map = parse_pauli_string("X(0)Y(1)");
 *   uint64_t basis_state = 0b01; // |01⟩
 *   auto [new_state, phase] = apply_Pk(pauli_map, basis_state);
 *   // new_state corresponds to the basis state after applying X on qubit 0 and Y on qubit 1
 */
std::pair<uint64_t, std::complex<double>> apply_Pk(
    const std::unordered_map<int, char>& pauli_map, uint64_t basis_state) {
    uint64_t new_basis_state = basis_state;
    std::complex<double> phase = 1.0;

    for (const auto& [qubit, p] : pauli_map) {
        if (qubit < 0 || qubit >= 64) {
            throw std::out_of_range("Qubit index out of valid range (0-63): " + std::to_string(qubit));
        }

        bool bit = (basis_state >> qubit) & 1ULL;

        switch (p) {
            case 'I':
                // Identity operator: do nothing
                break;
            case 'X':
                // Flip the qubit
                new_basis_state ^= (1ULL << qubit);
                break;
            case 'Y':
                // Flip the qubit and add phase factor
                new_basis_state ^= (1ULL << qubit);
                phase *= bit ? std::complex<double>(0, -1) : std::complex<double>(0, 1);
                break;
            case 'Z':
                // Add phase factor if the qubit is in state |1⟩
                if (bit) {
                    phase *= -1.0;
                }
                break;
            default:
                throw std::invalid_argument(std::string("Invalid Pauli operator: ") + p);
        }
    }
    return {new_basis_state, phase};
}

/**
 * Computes the inner product ⟨ψ|φ⟩ of two quantum states.
 *
 * @param psi The first quantum state (bra), represented as a State.
 * @param phi The second quantum state (ket), represented as a State.
 * @return The inner product ⟨ψ|φ⟩ as a std::complex<double>.
 */
std::complex<double> inner_product(const State& psi, const State& phi) {
    std::complex<double> result = 0.0;

    if (psi.empty() || phi.empty()) {
        throw std::invalid_argument("Quantum states cannot be empty.");
    }

    // Determine the smaller state for optimization
    const State* smaller_state = (psi.size() < phi.size()) ? &psi : &phi;
    const State* larger_state = (smaller_state == &psi) ? &phi : &psi;

    // Iterate over the basis states of the smaller state
    for (const auto& [basis_state, coeff_smaller] : *smaller_state) {
        auto it = larger_state->find(basis_state);
        if (it != larger_state->end()) {
            const std::complex<double>& coeff_larger = it->second;
            result += std::conj(coeff_smaller) * coeff_larger;
        }
    }

    return result;
}

/**
 * Computes the expectation value ⟨φ|H|ψ⟩ for a given Hamiltonian H represented as a sum of Pauli terms.
 *
 * @param phi The bra state ⟨φ|, represented as a State.
 * @param psi The ket state |ψ⟩, represented as a State.
 * @param H_terms A vector of pairs, each containing a Pauli string and its coefficient in H.
 * @return The expectation value ⟨φ|H|ψ⟩ as a std::complex<double>.
 */
std::complex<double> expectation_value(
    const State& phi, const State& psi,
    const std::vector<std::pair<std::string, double>>& H_terms) {
    if (psi.empty() || phi.empty()) {
        throw std::invalid_argument("Quantum states cannot be empty.");
    }

    if (H_terms.empty()) {
        throw std::invalid_argument("Hamiltonian terms cannot be empty.");
    }

    std::complex<double> result = 0.0;
    State psi_k;

    for (const auto& [pauli_string, coeff] : H_terms) {
        if (coeff == 0.0) {
            continue; // Skip zero coefficients
        }

        // Parse the Pauli string outside the loop
        const auto pauli_map = parse_pauli_string(pauli_string);
        psi_k.clear();

        // Apply Pk to ψ to obtain ψk
        for (const auto& [basis_state, amplitude] : psi) {
            auto [new_basis_state, phase] = apply_Pk(pauli_map, basis_state);
            psi_k[new_basis_state] += phase * amplitude;
        }
        // Compute ⟨φ|ψk⟩
        std::complex<double> inner_prod = inner_product(phi, psi_k);
        // Accumulate c_k * ⟨φ|ψk⟩
        result += coeff * inner_prod;
    }
    return result;
}

/**
 * Computes the expectation value ⟨φ|H|ψ⟩ in parallel.
 *
 * @param phi The bra state ⟨φ|, represented as a State.
 * @param psi The ket state |ψ⟩, represented as a State.
 * @param H_terms A vector of pairs, each containing a Pauli string and its coefficient in H.
 * @return The expectation value ⟨φ|H|ψ⟩ as a std::complex<double>.
 *
 * Note: This function parallelizes over the terms in H for performance on large Hamiltonians.
 */
std::complex<double> expectation_value_parallel(
    const State& phi, const State& psi,
    const std::vector<std::pair<std::string, double>>& H_terms) {
    if (psi.empty() || phi.empty()) {
        throw std::invalid_argument("Quantum states cannot be empty.");
    }

    if (H_terms.empty()) {
        throw std::invalid_argument("Hamiltonian terms cannot be empty.");
    }

    std::complex<double> result = 0.0;

    // Determine the number of threads to use
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;

    size_t total_terms = H_terms.size();
    size_t chunk_size = (total_terms + num_threads - 1) / num_threads;

    // Vector to hold partial results
    std::vector<std::complex<double>> partial_results(num_threads, 0.0);

    // Worker function for each thread
    auto worker = [&](size_t start_idx, size_t end_idx, unsigned int thread_id) {
        std::complex<double> local_result = 0.0;
        State psi_k;

        for (size_t idx = start_idx; idx < end_idx; ++idx) {
            const auto& [pauli_string, coeff] = H_terms[idx];
            if (coeff == 0.0) {
                continue; // Skip zero coefficients
            }

            // Parse the Pauli string once per term
            const auto pauli_map = parse_pauli_string(pauli_string);
            psi_k.clear();

            // Apply Pk to ψ to obtain ψk
            for (const auto& [basis_state, amplitude] : psi) {
                auto [new_basis_state, phase] = apply_Pk(pauli_map, basis_state);
                psi_k[new_basis_state] += phase * amplitude;
            }

            // Compute ⟨φ|ψk⟩
            std::complex<double> inner_prod = inner_product(phi, psi_k);

            // Accumulate c_k * ⟨φ|ψk⟩
            local_result += coeff * inner_prod;
        }
        partial_results[thread_id] = local_result;
    };

    // Launch threads
    std::vector<std::thread> threads;
    for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start_idx = t * chunk_size;
        size_t end_idx = std::min(start_idx + chunk_size, total_terms);
        if (start_idx >= end_idx) break;
        threads.emplace_back(worker, start_idx, end_idx, t);
    }

    // Join threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Combine partial results
    for (const auto& partial : partial_results) {
        result += partial;
    }

    return result;
}

/**
 * Applies the exponential of a Pauli operator e^{-iθP} to a quantum state |ψ⟩.
 *
 * @param pauli_string The Pauli operator P in string format, e.g., "X(0)Y(1)".
 * @param theta The angle θ in radians.
 * @param state The input quantum state |ψ⟩, represented as a State.
 * @return The new quantum state after applying e^{-iθP}|ψ⟩.
 */
State apply_exp_pauli(const std::string& pauli_string, double theta, const State& state) {
    if (state.empty()) {
        throw std::invalid_argument("Input quantum state cannot be empty.");
    }

    if (theta == 0.0) {
        return state; // No change if theta is zero
    }
    
    State new_state;
    const double cos_theta = std::cos(theta / 2.0);
    const double sin_theta = std::sin(theta / 2.0);
    const auto pauli_map = parse_pauli_string(pauli_string);

    // Compute P|ψ⟩ and update new_state
    for (const auto& [basis_state, amplitude] : state) {
        auto [new_basis_state, phase] = apply_Pk(pauli_map, basis_state);

        // Update new_state
        new_state[basis_state] += cos_theta * amplitude;
        new_state[new_basis_state] += (-std::complex<double>(0, sin_theta)) * phase * amplitude;
    }

    return new_state;
}

/**
 * Applies a sequence of exponential Pauli operators to an initial quantum state.
 *
 * @param U_terms A vector of ExpPauliTerm, each containing a Pauli string and an angle θ.
 * @param initial_state The initial quantum state, represented as a State.
 * @return The final quantum state after applying all exponential operators.
 */
State apply_U(const std::vector<ExpPauliTerm>& U_terms, const State& initial_state) {
    if (initial_state.empty()) {
        throw std::invalid_argument("Initial quantum state cannot be empty.");
    }

    State state = initial_state;
    for (const auto& term : U_terms) {
        if (term.pauli_string.empty()) {
            throw std::invalid_argument("Pauli string in U_terms cannot be empty.");
        }
        state = apply_exp_pauli(term.pauli_string, term.angle, state);
    }
    return state;
}

/**
 * Initializes the zero state |0⟩^n for a given number of qubits.
 *
 * @param num_qubits The number of qubits n.
 * @return The zero quantum state |0⟩^n, represented as a State.
 */
State initialize_zero_state(int num_qubits) {
    if (num_qubits <= 0 || num_qubits > 64) {
        throw std::invalid_argument("Number of qubits must be between 1 and 64.");
    }

    // Basis state |0...0⟩ corresponds to the basis state with index 0
    State state;
    state.emplace(0ULL, 1.0);
    return state;
}

/**
 * Converts a quantum state to a string representation for display or debugging.
 *
 * @param state The quantum state to convert, represented as a State.
 * @param num_qubits The number of qubits in the state.
 * @return A string representation of the quantum state.
 */
std::string state_to_string(const State& state, int num_qubits) {
    if (num_qubits <= 0 || num_qubits > 64) {
        throw std::invalid_argument("Number of qubits must be between 1 and 64.");
    }
    
    std::string s;
    for (const auto& [basis_state, coeff] : state) {
        s += "|";
        for (int i = num_qubits - 1; i >= 0; --i) {
            s += ((basis_state >> i) & 1ULL) ? "1" : "0";
        }
        s += "⟩: " + std::to_string(coeff.real()) +
             (coeff.imag() >= 0 ? "+" : "") + std::to_string(coeff.imag()) + "i\n";
    }
    return s;
}

// Bindings for Pybind11
PYBIND11_MODULE(spex, p) {
    p.doc() = "Expectation value computation module on sparse pauli states for tequila, implemented in C++ using Pybind11";

    // Expose ExpPauliTerm structure
    py::class_<ExpPauliTerm>(p, "ExpPauliTerm")
        .def(py::init<>())
        .def(py::init<const std::string&, double>(),
             py::arg("pauli_string"), py::arg("angle"))
        .def_readwrite("pauli_string", &ExpPauliTerm::pauli_string)
        .def_readwrite("angle", &ExpPauliTerm::angle);

 // Expose parse_pauli_string function
    p.def("parse_pauli_string", &parse_pauli_string,
          "Parse a Pauli string into a map from qubit indices to Pauli operators",
          py::arg("pauli_string"));


    // Expose apply_Pk function
    p.def("apply_Pk",
          [](const std::string& pauli_string, uint64_t basis_state) {
              const auto pauli_map = parse_pauli_string(pauli_string);
              return apply_Pk(pauli_map, basis_state);
          },
          "Apply a Pauli operator to a basis state",
          py::arg("pauli_string"), py::arg("basis_state"));

    // Expose inner_product function
    p.def("inner_product", &inner_product,
          "Compute the inner product of two quantum states",
          py::arg("psi"), py::arg("phi"));

    // Expose expectation_value function
    p.def("expectation_value", &expectation_value,
          "Compute the expectation value ⟨φ|H|ψ⟩",
          py::arg("phi"), py::arg("psi"), py::arg("H_terms"));

    // Expose expectation_value_parallel function
    p.def("expectation_value_parallel", &expectation_value_parallel,
          "Compute the expectation value ⟨φ|H|ψ⟩ (parallelized)",
          py::arg("phi"), py::arg("psi"), py::arg("H_terms"),
          py::call_guard<py::gil_scoped_release>());

    // Expose apply_exp_pauli function
    p.def("apply_exp_pauli", &apply_exp_pauli,
          "Apply e^{-iθP} to a quantum state",
          py::arg("pauli_string"), py::arg("theta"), py::arg("state"));

    // Expose apply_U function
    p.def("apply_U", &apply_U,
          "Apply a sequence of exponential Pauli operators to a quantum state",
          py::arg("U_terms"), py::arg("initial_state"));

    // Expose initialize_zero_state function
    p.def("initialize_zero_state", &initialize_zero_state,
          "Initialize the |0⟩ state",
          py::arg("num_qubits"));

    // Expose state_to_string function
    p.def("state_to_string", &state_to_string,
          "Convert a quantum state to a string representation",
          py::arg("state"), py::arg("num_qubits"));
}
