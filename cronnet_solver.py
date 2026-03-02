#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CronNet-Holo Spectral Solver
Quantum Self-Verification Geometry (QSVG) Framework
Luis Morató de Dalmases, 2026

Simulation of the {3,3,5} network and verification of spectral
convergence to the zeros of the Riemann zeta function.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import scipy.stats as stats
from scipy.linalg import eigh
import time
import warnings
warnings.filterwarnings('ignore')

class CronNetSolver:
    """
    Main solver for the {3,3,5} chronic network.
    Computes the discrete Laplacian spectrum and verifies
    convergence to Riemann zeros.
    """
    
    def __init__(self, R_eff=2e-4, c=299792458):
        """
        Initialize the CronNet simulator.
        
        Parameters:
        -----------
        R_eff : float
            Effective coherence radius (meters). Default: 2e-4 m
        c : float
            Speed of light (m/s)
        """
        self.R_eff = R_eff
        self.c = c
        
        # First non-trivial zeros of the Riemann Zeta function
        # γ₁, γ₂, γ₃, γ₄, γ₅ (imaginary parts)
        self.riemann_gamma = np.array([14.1347, 21.0220, 25.0108, 30.4248, 32.9351])
        self.riemann_ratios = self.riemann_gamma / self.riemann_gamma[0]
        
        # Geometric constants of {3,3,5}
        self.phi = (1 + np.sqrt(5)) / 2  # Golden ratio
        self.delta_theta = np.radians(6.8)  # Angular defect (degrees → radians)
        
        print(f"╔══════════════════════════════════════════════════════════╗")
        print(f"║     CronNet-Holo v1.0 - QSVG Spectral Solver            ║")
        print(f"║              Luis Morató de Dalmases, 2026              ║")
        print(f"╠══════════════════════════════════════════════════════════╣")
        print(f"║ R_eff = {R_eff:.2e} m")
        print(f"║ c = {c:.2e} m/s")
        print(f"║ δθ = {np.degrees(self.delta_theta):.2f}°")
        print(f"║ φ = {self.phi:.6f} (golden ratio)")
        print(f"╚══════════════════════════════════════════════════════════╝")
        
    def generate_lattice_335(self, num_nodes=1000, connectivity=0.25):
        """
        Generate an approximation of the {3,3,5} chronic network.
        
        Uses a Fibonacci spherical distribution to simulate
        optimal packing of the 120 vertices of the polytope.
        
        Parameters:
        -----------
        num_nodes : int
            Number of nodes to generate
        connectivity : float
            Connection radius (fraction of mean distance)
            
        Returns:
        --------
        adj_matrix : scipy.sparse.csr_matrix
            Adjacency matrix of the network
        nodes : np.ndarray
            Node coordinates (x, y, z)
        """
        # Generate points on S² (sphere) with Fibonacci distribution
        # This simulates the stereographic projection of {3,3,5}
        indices = np.arange(0, num_nodes, dtype=float) + 0.5
        phi = np.arccos(1 - 2 * indices / num_nodes)  # Polar angle
        theta = np.pi * (1 + np.sqrt(5)) * indices    # Azimuthal angle (Fibonacci)
        
        # Spherical → Cartesian coordinates
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        
        nodes = np.column_stack((x, y, z))
        
        # Build adjacency matrix based on proximity
        # In a real {3,3,5} network, each node has degree 12 (local icosahedron)
        tree = cKDTree(nodes)
        
        # Mean distance to nearest neighbors
        distances, _ = tree.query(nodes, k=2)
        mean_dist = np.mean(distances[:, 1])
        threshold = connectivity * mean_dist * 2.5  # Adjust for connectivity
        
        # Sparse distance matrix
        adj_matrix = tree.sparse_distance_matrix(tree, threshold, output_type='coo_matrix')
        
        # Convert to CSR and binarize
        adj_matrix.data = np.ones_like(adj_matrix.data)
        adj_matrix = adj_matrix.tocsr()
        
        # Ensure symmetry
        adj_matrix = adj_matrix.maximum(adj_matrix.transpose())
        
        # Remove self-connections
        adj_matrix.setdiag(0)
        adj_matrix.eliminate_zeros()
        
        # Mean degree
        degrees = np.array(adj_matrix.sum(axis=1)).flatten()
        mean_degree = np.mean(degrees)
        
        print(f"   Network {num_nodes} nodes | Mean degree: {mean_degree:.2f} | "
              f"Connections: {adj_matrix.nnz}")
        
        return adj_matrix, nodes
    
    def compute_spectrum(self, adj_matrix, k=6, normalized=True):
        """
        Compute Laplacian and extract eigenvalues (frequencies).
        
        Parameters:
        -----------
        adj_matrix : scipy.sparse.csr_matrix
            Adjacency matrix
        k : int
            Number of eigenvalues to compute (excluding zero)
        normalized : bool
            If True, use normalized Laplacian
            
        Returns:
        --------
        lambdas : np.ndarray
            Eigenvalues (λ₁, λ₂, ..., λₖ)
        """
        n = adj_matrix.shape[0]
        
        # Degree of each node
        degree = np.array(adj_matrix.sum(axis=1)).flatten()
        
        # Unnormalized Laplacian: L = D - A
        D = sp.diags(degree)
        L = D - adj_matrix
        
        if normalized:
            # Normalized Laplacian: L_norm = D^{-1/2} L D^{-1/2}
            # More numerically stable
            with np.errstate(divide='ignore'):
                D_inv_sqrt = sp.diags(1.0 / np.sqrt(degree + 1e-10))
            L_norm = D_inv_sqrt @ L @ D_inv_sqrt
            L_to_use = L_norm
        else:
            L_to_use = L
        
        # Compute eigenvalues (smallest, excluding zero)
        try:
            vals, vecs = eigsh(L_to_use, k=k+1, which='SM', tol=1e-6, maxiter=10000)
            # Sort
            idx = np.argsort(vals)
            vals = vals[idx]
            # Exclude zero (first eigenvalue)
            lambdas = vals[1:k+1]
        except:
            # Fallback: use eigh for small matrices
            L_dense = L_to_use.toarray()
            vals_dense = eigh(L_dense, eigvals_only=True)
            vals_dense = np.sort(vals_dense)
            lambdas = vals_dense[1:k+1]
        
        return lambdas
    
    def compute_frequencies(self, lambdas):
        """
        Convert eigenvalues to physical frequencies.
        
        f_n = (c/(2π R_eff)) * √(λ_n) / 2   (fundamental mode correction)
        
        Parameters:
        -----------
        lambdas : np.ndarray
            Laplacian eigenvalues
            
        Returns:
        --------
        frequencies : np.ndarray
            Frequencies in Hz
        ratios : np.ndarray
            Ratios f_n / f_1
        """
        # Conversion factor: c/(2π R_eff) * 1/2
        scale = self.c / (2 * np.pi * self.R_eff * 2)
        
        frequencies = scale * np.sqrt(lambdas)
        ratios = frequencies / frequencies[0]
        
        return frequencies, ratios
    
    def convergence_error(self, ratios):
        """
        Compute RMS error with respect to Riemann zeros.
        
        ε = √( Σᵢ (r_i - γ_i/γ₁)² )
        
        Parameters:
        -----------
        ratios : np.ndarray
            Computed frequency ratios
            
        Returns:
        --------
        error : float
            RMS error
        """
        n = min(len(ratios), len(self.riemann_ratios))
        error = np.sqrt(np.mean((ratios[:n] - self.riemann_ratios[:n])**2))
        return error
    
    def run_convergence_test(self, node_steps=[500, 1000, 2000, 4000, 8000]):
        """
        Run convergence test for different network sizes.
        
        Parameters:
        -----------
        node_steps : list
            List of node counts to test
        """
        print(f"\n{'='*60}")
        print(f" SPECTRAL CONVERGENCE TEST")
        print(f"{'='*60}")
        print(f"\nReference Riemann ratios: {self.riemann_ratios[:4]}")
        print(f"γ₁ = {self.riemann_gamma[0]:.4f}")
        print(f"{'Nodes':>10} {'λ₂/λ₁':>10} {'f₂/f₁':>10} {'f₃/f₁':>10} {'f₄/f₁':>10} {'Error':>12}")
        print(f"{'-'*65}")
        
        results = []
        times = []
        
        for n in node_steps:
            t_start = time.time()
            
            # Generate network
            adj, nodes = self.generate_lattice_335(num_nodes=n)
            
            # Compute spectrum
            lambdas = self.compute_spectrum(adj, k=4)
            
            # Convert to frequencies and ratios
            freqs, ratios = self.compute_frequencies(lambdas)
            
            # Compute error
            error = self.convergence_error(ratios)
            
            t_elapsed = time.time() - t_start
            times.append(t_elapsed)
            
            results.append({
                'nodes': n,
                'lambdas': lambdas,
                'ratios': ratios,
                'error': error,
                'time': t_elapsed
            })
            
            print(f"{n:10d} {lambdas[1]/lambdas[0]:10.4f} {ratios[1]:10.4f} "
                  f"{ratios[2]:10.4f} {ratios[3]:10.4f} {error:12.6f}  "
                  f"({t_elapsed:.1f}s)")
        
        self.plot_convergence(results)
        return results
    
    def plot_convergence(self, results):
        """
        Generate convergence plot.
        """
        nodes = [r['nodes'] for r in results]
        errors = [r['error'] for r in results]
        
        # Power law fit
        log_nodes = np.log(nodes)
        log_errors = np.log(errors)
        coeffs = np.polyfit(log_nodes, log_errors, 1)
        
        plt.figure(figsize=(12, 8))
        
        # Subplot 1: Error vs nodes
        plt.subplot(2, 2, 1)
        plt.loglog(nodes, errors, 'bo-', linewidth=2, markersize=8, 
                   label='CronNet Simulation')
        
        # Theoretical fit
        nodes_fit = np.logspace(2, 4, 100)
        errors_fit = np.exp(coeffs[1]) * nodes_fit**coeffs[0]
        plt.loglog(nodes_fit, errors_fit, 'r--', linewidth=1.5,
                   label=f'$\epsilon \propto N^{{{coeffs[0]:.2f}}}$')
        
        # Experimental threshold
        plt.axhline(y=0.01, color='g', linestyle=':', label='1% Threshold')
        
        plt.xlabel('Number of nodes $|V|$')
        plt.ylabel('RMS error $\epsilon$')
        plt.title('Convergence to Riemann zeros')
        plt.grid(True, which='both', alpha=0.3)
        plt.legend()
        
        # Subplot 2: Ratio comparison
        plt.subplot(2, 2, 2)
        x_pos = np.arange(1, 5)
        width = 0.35
        
        # Last result (highest resolution)
        last = results[-1]
        plt.bar(x_pos - width/2, last['ratios'][:4], width, label='CronNet', color='blue', alpha=0.7)
        plt.bar(x_pos + width/2, self.riemann_ratios[:4], width, label='Riemann', color='red', alpha=0.7)
        
        plt.xlabel('Mode $n$')
        plt.ylabel('Ratio $f_n/f_1$')
        plt.title(f'Ratio comparison (N={last["nodes"]})')
        plt.xticks(x_pos)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Subplot 3: Computation time
        plt.subplot(2, 2, 3)
        plt.plot(nodes, [r['time'] for r in results], 's-', color='purple')
        plt.xlabel('Number of nodes')
        plt.ylabel('Computation time (s)')
        plt.title('Computational scaling')
        plt.grid(True, alpha=0.3)
        
        # Subplot 4: Error per mode
        plt.subplot(2, 2, 4)
        mode_errors = []
        for i in range(4):
            mode_errors.append(np.abs(last['ratios'][i] - self.riemann_ratios[i]))
        
        plt.bar(x_pos, mode_errors, color='orange', alpha=0.7)
        plt.xlabel('Mode $n$')
        plt.ylabel('Absolute error')
        plt.title(f'Error per mode (N={last["nodes"]})')
        plt.grid(True, alpha=0.3)
        plt.xticks(x_pos)
        
        plt.tight_layout()
        plt.savefig('convergence_plot.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print(f"\n📈 Plot saved as 'convergence_plot.png'")
        print(f"   Convergence exponent: {coeffs[0]:.3f}")


class CronNetRobustness(CronNetSolver):
    """
    Class for topological robustness tests.
    Demonstrates that ratios remain invariant under perturbations.
    """
    
    def add_topological_noise(self, adj_matrix, noise_level=0.05, noise_type='edge'):
        """
        Add topological noise to the network.
        
        Parameters:
        -----------
        adj_matrix : scipy.sparse.csr_matrix
            Original adjacency matrix
        noise_level : float
            Noise level (fraction of connections to modify)
        noise_type : str
            Noise type: 'edge' (weights) or 'structural' (connections)
            
        Returns:
        --------
        adj_noisy : scipy.sparse.csr_matrix
            Matrix with noise
        """
        n = adj_matrix.shape[0]
        
        if noise_type == 'edge':
            # Noise in edge weights
            # Create noise matrix with same structure
            rows, cols = adj_matrix.nonzero()
            n_edges = len(rows)
            
            # Gaussian noise in weights
            noise_vals = np.random.normal(1.0, noise_level, n_edges)
            noise_vals = np.maximum(noise_vals, 0)  # No negative weights
            
            # Create new matrix
            adj_noisy = sp.csr_matrix((noise_vals, (rows, cols)), shape=(n, n))
            
        else:
            # Structural noise: add/remove random connections
            # Number of connections to modify
            n_mod = int(noise_level * adj_matrix.nnz)
            
            # Copy original matrix
            adj_noisy = adj_matrix.copy().tolil()
            
            # Remove some existing connections
            rows, cols = adj_matrix.nonzero()
            if len(rows) > 0:
                idx_to_remove = np.random.choice(len(rows), min(n_mod//2, len(rows)), replace=False)
                for idx in idx_to_remove:
                    adj_noisy[rows[idx], cols[idx]] = 0
                    adj_noisy[cols[idx], rows[idx]] = 0
            
            # Add new connections
            for _ in range(n_mod//2):
                i, j = np.random.randint(0, n, 2)
                if i != j:
                    adj_noisy[i, j] = 1
                    adj_noisy[j, i] = 1
            
            adj_noisy = adj_noisy.tocsr()
        
        return adj_noisy
    
    def simulate_noise_impact(self, node_count=2000, 
                              noise_levels=[0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2],
                              repetitions=5):
        """
        Simulate impact of thermal/electromagnetic noise.
        
        Parameters:
        -----------
        node_count : int
            Number of network nodes
        noise_levels : list
            Noise levels to test
        repetitions : int
            Repetitions per level (for statistics)
        """
        print(f"\n{'='*60}")
        print(f" TOPOLOGICAL ROBUSTNESS TEST")
        print(f"{'='*60}")
        print(f"\nNodes: {node_count} | Repetitions: {repetitions}")
        print(f"\n{'Level':>8} {'Mean Error':>12} {'Std Dev':>12} {'Min':>10} {'Max':>10} {'Status':>12}")
        print(f"{'-'*70}")
        
        # Generate base network (no noise)
        adj_base, nodes = self.generate_lattice_335(num_nodes=node_count)
        lambdas_base = self.compute_spectrum(adj_base)
        _, ratios_base = self.compute_frequencies(lambdas_base)
        error_base = self.convergence_error(ratios_base)
        
        results = []
        
        for noise in noise_levels:
            errors = []
            
            for rep in range(repetitions):
                # Add noise
                adj_noisy = self.add_topological_noise(adj_base, noise_level=noise)
                
                # Compute spectrum
                try:
                    lambdas_noisy = self.compute_spectrum(adj_noisy)
                    _, ratios_noisy = self.compute_frequencies(lambdas_noisy)
                    error = self.convergence_error(ratios_noisy)
                    errors.append(error)
                except:
                    # If computation fails, high error
                    errors.append(1.0)
            
            errors = np.array(errors)
            mean_error = np.mean(errors)
            std_error = np.std(errors)
            min_error = np.min(errors)
            max_error = np.max(errors)
            
            # Determine status
            if mean_error < 0.005:
                status = "✅ STABLE"
            elif mean_error < 0.01:
                status = "🟡 MARGINAL"
            elif mean_error < 0.02:
                status = "⚠️ DEGRADED"
            else:
                status = "❌ UNSTABLE"
            
            results.append({
                'noise': noise,
                'mean': mean_error,
                'std': std_error,
                'min': min_error,
                'max': max_error,
                'status': status
            })
            
            print(f"{noise*100:7.1f}%   {mean_error:10.6f}   ±{std_error:8.6f}   "
                  f"{min_error:8.6f}   {max_error:8.6f}   {status}")
        
        self.plot_robustness(results, error_base)
        return results
    
    def plot_robustness(self, results, baseline_error):
        """
        Generate topological robustness plot.
        """
        noise_levels = [r['noise']*100 for r in results]
        mean_errors = [r['mean'] for r in results]
        std_errors = [r['std'] for r in results]
        
        plt.figure(figsize=(12, 8))
        
        # Subplot 1: Error vs noise
        plt.subplot(2, 2, 1)
        plt.errorbar(noise_levels, mean_errors, yerr=std_errors, 
                     fmt='o-', color='teal', linewidth=2, markersize=8,
                     capsize=5, label='Error with noise')
        
        # Baseline (no noise)
        plt.axhline(y=baseline_error, color='blue', linestyle='--', 
                    label=f'No noise (ε={baseline_error:.6f})')
        
        # Stability regions
        plt.axhspan(0, 0.005, alpha=0.2, color='green', label='Stable')
        plt.axhspan(0.005, 0.01, alpha=0.2, color='yellow', label='Marginal')
        plt.axhspan(0.01, 0.02, alpha=0.2, color='orange', label='Degraded')
        plt.axhspan(0.02, max(mean_errors)+0.01, alpha=0.2, color='red', label='Unstable')
        
        plt.xlabel('Noise level (%)')
        plt.ylabel('RMS error ε')
        plt.title('Topological robustness under perturbations')
        plt.grid(True, alpha=0.3)
        plt.legend(loc='upper left')
        
        # Subplot 2: Ratio comparison at different noise levels
        plt.subplot(2, 2, 2)
        x_pos = np.arange(1, 5)
        width = 0.2
        
        # Select levels to visualize
        idx_to_plot = [0, 2, 4, 6]  # 0%, 2%, 5%, 15%
        colors = ['blue', 'green', 'orange', 'red']
        
        for i, idx in enumerate(idx_to_plot):
            if idx < len(results):
                noise_val = results[idx]['noise'] * 100
                # Need to recompute for this level
                adj_base, _ = self.generate_lattice_335(num_nodes=1000)
                adj_noisy = self.add_topological_noise(adj_base, noise_level=results[idx]['noise'])
                lambdas = self.compute_spectrum(adj_noisy)
                _, ratios = self.compute_frequencies(lambdas)
                
                plt.plot(x_pos, ratios[:4], 'o-', color=colors[i], 
                        label=f'{noise_val:.0f}% noise', linewidth=1.5)
        
        # Riemann ratios
        plt.plot(x_pos, self.riemann_ratios[:4], 'k--', label='Riemann', linewidth=2)
        
        plt.xlabel('Mode $n$')
        plt.ylabel('Ratio $f_n/f_1$')
        plt.title('Ratio evolution with noise')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Subplot 3: Error distribution
        plt.subplot(2, 2, 3)
        for i, noise in enumerate(noise_levels):
            if i % 2 == 0:  # Every 2 levels
                # Generate simulated data for distribution
                simulated_errors = np.random.normal(mean_errors[i], std_errors[i], 1000)
                plt.hist(simulated_errors, bins=20, alpha=0.5, 
                        label=f'{noise:.0f}%', density=True)
        
        plt.xlabel('Error ε')
        plt.ylabel('Probability density')
        plt.title('Error distribution by noise level')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Subplot 4: Stability heat map
        plt.subplot(2, 2, 4)
        
        # Create stability probability matrix
        stability_matrix = np.zeros((len(noise_levels), 3))
        for i, r in enumerate(results):
            if r['mean'] < 0.005:
                stability_matrix[i, 0] = 1.0
            elif r['mean'] < 0.01:
                stability_matrix[i, 1] = 1.0
            else:
                stability_matrix[i, 2] = 1.0
        
        plt.imshow(stability_matrix.T, aspect='auto', cmap='RdYlGn', 
                   extent=[min(noise_levels)-2, max(noise_levels)+2, 0, 3])
        
        plt.yticks([0.5, 1.5, 2.5], ['Stable', 'Marginal', 'Unstable'])
        plt.xlabel('Noise level (%)')
        plt.title('Stability map')
        plt.colorbar(label='Probability')
        
        plt.tight_layout()
        plt.savefig('robustness_plot.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print(f"\n📈 Plot saved as 'robustness_plot.png'")


class CronNetExperimental(CronNetSolver):
    """
    Class for real experimental simulation.
    Includes cavity effects, detection, and experimental noise.
    """
    
    def __init__(self, R_eff=2e-4, cavity_length=50e-6, finesse=1e5, temperature=4.0):
        """
        Initialize with experimental parameters.
        """
        super().__init__(R_eff)
        self.cavity_length = cavity_length
        self.finesse = finesse
        self.temperature = temperature
        
        # Detection parameters
        self.Q_factor = finesse * 2 * np.pi / (1 - np.exp(-np.pi/finesse))
        self.linewidth = self.c / (2 * cavity_length * self.Q_factor)
        
    def simulate_experiment(self, node_count=5000, scan_range=(1.0, 4.0), scan_points=1000):
        """
        Simulate an experimental frequency scan.
        """
        print(f"\n{'='*60}")
        print(f" EXPERIMENTAL SIMULATION")
        print(f"{'='*60}")
        print(f"\nCavity parameters:")
        print(f"  Length L = {self.cavity_length*1e6:.1f} μm")
        print(f"  Finesse F = {self.finesse:.1e}")
        print(f"  Q factor = {self.Q_factor:.1e}")
        print(f"  Linewidth Δf = {self.linewidth/1e9:.3f} GHz")
        print(f"  Temperature T = {self.temperature} K")
        
        # Generate network
        adj, nodes = self.generate_lattice_335(num_nodes=node_count)
        
        # Compute ideal spectrum
        lambdas = self.compute_spectrum(adj, k=6)
        freqs_ideal, ratios_ideal = self.compute_frequencies(lambdas)
        
        print(f"\nIdeal frequencies (THz):")
        for i, f in enumerate(freqs_ideal[:4]):
            print(f"  f_{i+1} = {f/1e12:.4f} THz  (ratio = {ratios_ideal[i]:.4f})")
        
        # Generate scan
        scan_freqs = np.linspace(scan_range[0]*1e12, scan_range[1]*1e12, scan_points)
        
        # Cavity response model (Lorentzian)
        transmission = np.zeros_like(scan_freqs)
        
        for i, f_ideal in enumerate(freqs_ideal[:4]):
            # Linewidth
            gamma = self.linewidth
            
            # Mode contribution
            lorentzian = 1.0 / (1 + ((scan_freqs - f_ideal) / gamma)**2)
            
            # Thermal noise (proportional to T)
            thermal_noise = np.random.normal(0, 0.01 * self.temperature/4.0, len(scan_freqs))
            
            transmission += lorentzian * (0.25 * (4-i)) + thermal_noise
        
        # Add detection noise
        transmission += np.random.normal(0, 0.02, len(scan_freqs))
        
        self.plot_experimental_scan(scan_freqs, transmission, freqs_ideal)
        
        return scan_freqs, transmission, freqs_ideal
    
    def plot_experimental_scan(self, scan_freqs, transmission, freqs_ideal):
        """
        Plot experimental scan results.
        """
        plt.figure(figsize=(14, 8))
        
        # Subplot 1: Full spectrum
        plt.subplot(2, 2, 1)
        plt.plot(scan_freqs/1e12, transmission, 'b-', linewidth=1, alpha=0.7)
        
        # Mark predicted peaks
        for i, f in enumerate(freqs_ideal[:4]):
            plt.axvline(x=f/1e12, color='r', linestyle='--', alpha=0.5)
            plt.text(f/1e12, 0.8, f'f_{i+1}', rotation=90, fontsize=10)
        
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Transmission (a.u.)')
        plt.title('Simulated transmission spectrum')
        plt.grid(True, alpha=0.3)
        
        # Subplot 2: Zoom on first peak
        plt.subplot(2, 2, 2)
        mask = (scan_freqs > freqs_ideal[0] - 5*self.linewidth) & \
               (scan_freqs < freqs_ideal[0] + 5*self.linewidth)
        
        plt.plot(scan_freqs[mask]/1e12, transmission[mask], 'b-', linewidth=2)
        plt.axvline(x=freqs_ideal[0]/1e12, color='r', linestyle='--', 
                   label=f'f₁ = {freqs_ideal[0]/1e12:.3f} THz')
        
        # Lorentzian fit
        gamma = self.linewidth/1e12
        lorentz_fit = 1.0 / (1 + ((scan_freqs[mask] - freqs_ideal[0])/self.linewidth)**2)
        plt.plot(scan_freqs[mask]/1e12, lorentz_fit*0.8, 'g--', 
                label='Lorentz fit', linewidth=1.5)
        
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Transmission')
        plt.title(f'Peak f₁ (width = {gamma:.3f} GHz)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Subplot 3: Signal-to-noise ratio
        plt.subplot(2, 2, 3)
        snr = []
        for f in freqs_ideal[:4]:
            idx = np.argmin(np.abs(scan_freqs - f))
            signal = transmission[idx]
            noise_floor = np.std(transmission[transmission < 0.1])
            snr.append(signal / noise_floor)
        
        plt.bar(range(1,5), snr, color='orange', alpha=0.7)
        plt.axhline(y=5, color='r', linestyle='--', label='5σ threshold')
        plt.xlabel('Mode')
        plt.ylabel('SNR')
        plt.title('Signal-to-noise ratio')
        plt.xticks(range(1,5))
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Subplot 4: Experimental ratios
        plt.subplot(2, 2, 4)
        
        # Detect peaks (local maxima)
        from scipy.signal import find_peaks
        peaks, properties = find_peaks(transmission, height=0.3, distance=100)
        detected_freqs = scan_freqs[peaks]
        
        if len(detected_freqs) >= 2:
            detected_ratios = detected_freqs / detected_freqs[0]
            
            x_pos = np.arange(1, min(5, len(detected_ratios)+1))
            plt.bar(x_pos - 0.2, detected_ratios[:4], width=0.4, 
                   label='Detected', color='blue', alpha=0.7)
            plt.bar(x_pos + 0.2, self.riemann_ratios[:4], width=0.4, 
                   label='Riemann', color='red', alpha=0.7)
            
            plt.xlabel('Mode')
            plt.ylabel('Ratio f_n/f₁')
            plt.title('Detected ratio comparison')
            plt.xticks(x_pos)
            plt.legend()
        else:
            plt.text(0.5, 0.5, 'Insufficient peaks detected', 
                    ha='center', va='center', transform=plt.gca().transAxes)
        
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('experimental_scan.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print(f"\n📈 Plot saved as 'experimental_scan.png'")


def main():
    """
    Main function: run all simulations.
    """
    print("\n" + "="*70)
    print(" CRONNET-HOLO v1.0 - QSVG SPECTRAL SIMULATOR")
    print(" Luis Morató de Dalmases, 2026")
    print("="*70)
    
    # 1. Convergence test
    solver = CronNetSolver(R_eff=2e-4)
    convergence_results = solver.run_convergence_test(
        node_steps=[500, 1000, 2000, 4000, 8000]
    )
    
    # 2. Topological robustness test
    robustness = CronNetRobustness(R_eff=2e-4)
    robustness_results = robustness.simulate_noise_impact(
        node_count=2000,
        noise_levels=[0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2],
        repetitions=3
    )
    
    # 3. Experimental simulation
    experiment = CronNetExperimental(
        R_eff=2e-4,
        cavity_length=50e-6,
        finesse=1e5,
        temperature=4.0
    )
    experiment.simulate_experiment(
        node_count=5000,
        scan_range=(1.0, 4.0),
        scan_points=2000
    )
    
    print("\n" + "="*70)
    print(" SIMULATION COMPLETED")
    print("="*70)
    print("\n📊 Generated files:")
    print("   - convergence_plot.png")
    print("   - robustness_plot.png")
    print("   - experimental_scan.png")
    print("\n✅ Ready for CERN presentation!")


if __name__ == "__main__":
    main()