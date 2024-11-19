#include <iostream>
#include <unordered_map>
#include <complex>
#include <string>
#include <cmath>
#include <cstdint>
#include <vector>

// represents a quantum state, maps basis state to its complex amplitude
using State = std::unordered_map<uint64_t, std::complex<double>>;

// Function for applying pauli operator Pk to a given basis state
std::pair<uint64_t, std::complex<double>> apply_Pk(const std::string& pauli_string, uint64_t basis_state) {
    uint64_t new_basis_state = basis_state;
    std::complex<double> phase = 1.0;

    for (size_t i = 0; i < pauli_string.length(); ++i) {
        char p = pauli_string[i];
        bool bit = (basis_state >> i) & 1;

        switch (p) {
            case 'I':
                // identity operator
                break;
            case 'X':
                // flips the qubit
                new_basis_state ^= (1ULL << i);
                break;
            case 'Y':
                // flips the qubit and adds a phase factor
                new_basis_state ^= (1ULL << i);
                phase *= bit ? std::complex<double>(0, -1) : std::complex<double>(0, 1);
                break;
            case 'Z':
                // adds a phase factor if the qubit is in state |1⟩
                if (bit) {
                    phase *= -1.0;
                }
                break;
            default:
                std::cerr << "Invalid Pauli operator: " << p << std::endl;
                exit(1);
        }
    }
    return {new_basis_state, phase};
}

// Function for calculating the inner product of two quantum states
std::complex<double> inner_product(const State& psi, const State& phi) {
    std::complex<double> result = 0.0;

    // Pointers to avoid copying states
    const State* smaller_state;
    const State* larger_state;

    // Determine which state is smaller
    if (psi.size() < phi.size()) {
        smaller_state = &psi;
        larger_state = &phi;
    } else {
        smaller_state = &phi;
        larger_state = &psi;
    }

    // Iterate over the basis states and coefficients of the smaller state
    for (const auto& [basis_state, coeff_smaller] : *smaller_state) {
        // Try to find the matching basis state in the larger state
        auto it = larger_state->find(basis_state);
        if (it != larger_state->end()) {
            // it exists in smaller and larger state
            // Extract coefficient(second) from larger state
            std::complex<double> coeff_larger = it->second;
            if (smaller_state == &phi) {
                result += std::conj(coeff_smaller) * coeff_larger;
            } else {
                result += std::conj(coeff_larger) * coeff_smaller;
            }
        }
    }

    return result;
}


int main() {
    // List of tuples (Pauli-string, theta)
    std::vector<std::pair<std::string, double>> input = {
        {"ZYXZYX", 1.5708},
        {"XYZXYZ", 3.1416}
    };

    size_t num_sequences = input.size();
    if (num_sequences == 0) {
        std::cerr << "No pauli-strings entered." << std::endl;
        return 1;
    }

    // Check, if pauli-strings have the same length
    size_t num_qubits = input[0].first.length();
    for (size_t i = 1; i < num_sequences; ++i) {
        if (input[i].first.length() != num_qubits) {
            std::cerr << "Pauli-strings have to be the same length. May use identity operator 'I'." << std::endl;
            return 1;
        }
    }

    // Initial sate |0⟩
    State psi;
    psi[0] = 1.0;

    // Applying the gate sequences
    for (size_t seq = 0; seq < num_sequences; ++seq) {
        const std::string& pauli_string = input[seq].first;
        double theta = input[seq].second;

        State new_psi;
        double cos_theta = std::cos(theta / 2.0);
        double sin_theta = std::sin(theta / 2.0);

        for (const auto& [basis_state, coeff] : psi) {
            // cos(θ/2) * |l⟩
            new_psi[basis_state] += cos_theta * coeff;

            // -i sin(θ/2) * Pk |l⟩
            auto [new_basis_state, phase] = apply_Pk(pauli_string, basis_state);
            new_psi[new_basis_state] += (-std::complex<double>(0, 1) * sin_theta) * phase * coeff;
        }

        // update the state for the next sequence
        psi = std::move(new_psi);
    }

    // output
    for (const auto& [basis_state, coeff] : psi) {
        std::cout << "|";
        for (int i = num_qubits - 1; i >= 0; --i) {
            std::cout << ((basis_state >> i) & 1);
        }
        std::cout << "⟩: " << coeff << std::endl;
    }


    // Testing inner_product
    State phi;

    // Initialize psi and phi with basis states and coefficients
    psi[0b00] = std::complex<double>(1.0, 0.0); // |00⟩
    psi[0b11] = std::complex<double>(0.0, 1.0); // |11⟩

    phi[0b00] = std::complex<double>(0.5, 0.5); // |00⟩
    phi[0b01] = std::complex<double>(-0.5, 0.5); // |01⟩

    // Calculate the inner product
    std::complex<double> result = inner_product(psi, phi);

    std::cout << "Inner product ⟨φ|ψ⟩ = " << result << std::endl;

    return 0;

    // c++ template<>

    // pybind11
}