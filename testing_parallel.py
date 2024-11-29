import tequila as tq
import spex_tequila as spex
import numpy as np
import time

def term_to_string(term):
    """Convert a Tequila PauliString to a string like 'X(0)Y(1)Z(2)'."""
    if term == ():  # Identity operator
        return 'I'
    pauli_string = ''.join(f"{pauli}({qubit})" for qubit, pauli in term)
    return pauli_string

def convert_H_to_H_terms(H):
    """Convert a Tequila Hamiltonian to a list of (pauli_string, coeff)."""
    return [(term_to_string(term), coeff) for term, coeff in H.items()]

def state_vector_to_State(psi_wavefunction, num_qubits=None):
    """Convert a QubitWaveFunction to the C++ State format (dict of index to amplitude)."""
    State = {}
    # Determine number of qubits, if not given
    if num_qubits is None:
        max_bitstring_length = max(len(state.split('>')[0].split('|')[-1]) for state in str(psi_wavefunction).split('>') if state)
        num_qubits = max_bitstring_length
    
    for state in str(psi_wavefunction).split('>'):
        if '|' not in state:
            continue
        amplitude_str, bitstring = state.strip().split('|')
        amplitude = complex(amplitude_str.strip().strip('()'))
        # change bitstring to little endian
        bitstring = bitstring.strip()
        index = int(bitstring[::-1], 2)
        if abs(amplitude) > 1e-12:
            State[index] = amplitude
    return State


def expectation_value_tequila(n):
    R = 1.5
    geom = ""
    for k in range(2*n):
        geom += "h 0.0 0.0 {}\n".format(R*k)

    edges = [(2*i, 2*i+1) for i in range(n)]

    mol = tq.Molecule(geometry=geom, basis_set="sto-3g")
    U = mol.make_ansatz(name="SPA", edges=edges)
    U = U.map_variables({k:1.0 for k in U.extract_variables()})

    H = mol.make_hamiltonian()

    start_time = time.time()
    E = tq.ExpectationValue(H=H, U=U)
    E = tq.simulate(E)
    end_time = time.time()
    total_time = end_time - start_time

    print(f"Expectation Value: {E}")
    print(f"Total Time: {total_time:.6f} seconds")

def expectation_value_spex(n):
    R = 1.5
    geom = ""
    for k in range(2*n):
        geom += "h 0.0 0.0 {}\n".format(R*k)

    edges = [(2*i, 2*i+1) for i in range(n)]

    mol = tq.Molecule(geometry=geom, basis_set="sto-3g")
    U = mol.make_ansatz(name="SPA", edges=edges)
    U = U.map_variables({k:1.0 for k in U.extract_variables()})

    H = mol.make_hamiltonian()

    start_time_simulate = time.time()
    psi_vector = tq.simulate(U)
    end_time_simulate = time.time()
    total_time_simulate = end_time_simulate - start_time_simulate

    start_time_conversion = time.time()

    psi_state = state_vector_to_State(psi_vector, n**2)
    phi_state = psi_state

    H_terms = convert_H_to_H_terms(H)
    end_time_conversion = time.time()
    total_time_conversion = end_time_conversion - start_time_conversion

    start_time_spex = time.time()
    E_spex = spex.expectation_value(phi_state, psi_state, H_terms)
    end_time_spex = time.time()
    total_time_spex = end_time_spex - start_time_spex

    total_time = total_time_conversion + total_time_spex + total_time_simulate

    print(f"Expectation Value: {E_spex}")
    print(f"Simulate Time: {total_time_simulate:.6f} seconds")
    print(f"Conversion Time: {total_time_conversion:.6f} seconds")
    print(f"spex Computation Time: {total_time_spex:.6f} seconds")
    print(f"Total Time: {total_time:.6f} seconds")

def expectation_value_spex_parallel(n):
    R = 1.5
    geom = ""
    for k in range(2*n):
        geom += "h 0.0 0.0 {}\n".format(R*k)

    edges = [(2*i, 2*i+1) for i in range(n)]

    mol = tq.Molecule(geometry=geom, basis_set="sto-3g")
    U = mol.make_ansatz(name="SPA", edges=edges)
    U = U.map_variables({k:1.0 for k in U.extract_variables()})

    H = mol.make_hamiltonian()

    start_time_simulate = time.time()
    psi_vector = tq.simulate(U)
    end_time_simulate = time.time()
    total_time_simulate = end_time_simulate - start_time_simulate

    start_time_conversion = time.time()

    psi_state = state_vector_to_State(psi_vector, n**2)
    phi_state = psi_state

    H_terms = convert_H_to_H_terms(H)
    end_time_conversion = time.time()
    total_time_conversion = end_time_conversion - start_time_conversion

    start_time_spex = time.time()
    E_spex = spex.expectation_value_parallel(phi_state, psi_state, H_terms)
    end_time_spex = time.time()
    total_time_spex = end_time_spex - start_time_spex

    total_time = total_time_conversion + total_time_spex + total_time_simulate

    print(f"Expectation Value: {E_spex}")
    print(f"Simulate Time: {total_time_simulate:.6f} seconds")
    print(f"Conversion Time: {total_time_conversion:.6f} seconds")
    print(f"spex Parallel Computation Time: {total_time_spex:.6f} seconds")
    print(f"Total Time: {total_time:.6f} seconds")


def main():
    n = 4

    print("\n### Tequila ###")
    expectation_value_tequila(n)
    
    print("\n### spex ###")
    expectation_value_spex(n)
    
    print("\n### spex-parallel ###")
    expectation_value_spex_parallel(n)
    
    print("\n")


if __name__ == "__main__":
    main()
